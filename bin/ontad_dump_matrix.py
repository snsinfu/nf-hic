#!/usr/bin/env python3
"""Dump a single-chromosome contact matrix for OnTAD.

Writes an N*N space-separated text matrix that OnTAD can read. Only the
diagonal band |i-j| <= --diagonal is fetched from the (m)cooler; values outside
the band are written as --pad-value (default 0).

Missing pixels (NaN, e.g. masked/unbalanced bins) can be left as-is
(--nan-mode zero, the default) or replaced by the per-diagonal P(s)
expectation estimated from the finite values of the same chromosome
(--nan-mode expected).
"""

import argparse
import logging
import signal
import sys
import typing

import cooler
import numpy as np


LOG = logging.getLogger(__name__)
LOG_LEVEL = logging.INFO
LOG_FORMAT = "[%(asctime)s] %(message)s"
LOG_DATE_FORMAT = "%F %T"


def main(
    *,
    cool: str,
    output: str | None,
    region: str,
    bin_size: int | None,
    nan_mode: str | None,
    diagonal: int | None,
    no_balance: bool,
    pad_value: str,
) -> None:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    if bin_size is not None:
        cool += f"::/resolutions/{bin_size}"

    clr = cooler.Cooler(cool)

    chrom = region.split(":", 1)[0]
    if chrom not in set(clr.chromnames):
        LOG.error(
            "Region '%s' is not present in %s (available: %s)",
            region,
            cool,
            ", ".join(clr.chromnames),
        )
        sys.exit(1)

    region_start, region_end = clr.extent(region)
    region_size = region_end - region_start
    if region_size == 0:
        LOG.error("Region '%s' resolves to 0 bins in %s", region, cool)
        sys.exit(1)

    clr_matrix = clr.matrix(balance=(not no_balance))
    if not no_balance and "weight" not in clr.bins().columns:
        LOG.error(
            "%s has no stored balancing weights; use --no-balance or balance the cooler first",
            cool,
        )
        sys.exit(1)

    LOG.info("Region (%s): %d bins", region, region_size)

    # Default is to load the full matrix.
    if diagonal is None:
        diagonal = region_size - 1

    expectation: np.ndarray | None = None
    if nan_mode == "zero":
        expectation = np.zeros(diagonal + 1)
    elif nan_mode == "expected":
        LOG.info("Estimating contact expectation")
        expectation = _estimate_expectation(
            clr_matrix,
            region_start,
            region_end,
            diagonal,
        )
        if not np.all(np.isfinite(expectation)):
            LOG.warning(
                "Some diagonals of '%s' have no finite pixels; their NaNs are left "
                "as 'nan' in the output (OnTAD reads them as 0)",
                region,
            )

    LOG.info("Writing matrix out")

    pad_size = len(pad_value) + 1
    pad_seq = (pad_value + " ") * region_size
    pad_seq_left = pad_seq[pad_size:]
    pad_seq_right = pad_seq[pad_size - 1 : -1]

    with _open_or_stdout(output, "wt") as output_file:
        for index, center, band in scan_diagonal_bands(
            clr_matrix,
            region_start,
            region_end,
            band_size=diagonal,
            expectation=expectation,
        ):
            n_pads_left = index - center
            n_pads_right = region_size - len(band) - n_pads_left

            output_file.write(pad_seq_left[: n_pads_left * pad_size])
            output_file.write(" ".join(_format(value) for value in band))
            output_file.write(pad_seq_right[: n_pads_right * pad_size])
            output_file.write("\n")

    LOG.info("Done")


def _format(value: float) -> str:
    # NaNs are emitted as 'nan' on purpose: the explicit --nan-mode zero fills
    # them, while --nan-mode expected only fills diagonals that have data.
    # Whatever remains is a real gap and is kept for audit. OnTAD's C parser
    # reads 'nan' as 0.
    return "{:g}".format(value)


def _open_or_stdout(filename: str | None, mode: str) -> typing.IO:
    if filename is None:
        return open(sys.stdout.fileno(), mode, closefd=False)
    return open(filename, mode)


class DiagonalBand(typing.NamedTuple):
    position: int
    center: int
    band: np.ndarray


def scan_diagonal_bands(
    matrix,
    start: int,
    end: int,
    *,
    band_size: int,
    chunk_size: int | None = None,
    expectation: np.ndarray | None = None,
) -> typing.Iterator[DiagonalBand]:
    band_width = band_size * 2 + 1

    if chunk_size is None:
        chunk_size = band_width * 2
    if chunk_size < band_width:
        chunk_size = band_width

    step_size = chunk_size - band_width + 1

    for chunk_start in range(start, end, step_size):
        chunk_end = min(chunk_start + chunk_size, end)
        chunk_slice = slice(chunk_start, chunk_end)
        chunk = matrix[chunk_slice, chunk_slice]

        iter_start = band_size
        iter_end = chunk_size - band_size

        if chunk_start == start:
            iter_start = 0
        if chunk_end == end:
            iter_end = len(chunk)

        chunk_offset = chunk_start - start

        for i in range(iter_start, iter_end):
            band_start = max(i - band_size, 0)
            band_end = min(i + band_size + 1, len(chunk))
            band = chunk[i, band_start:band_end]
            center = i - band_start
            position = chunk_offset + i
            if expectation is not None:
                _nan_to_expectation(band, center, expectation)
            yield DiagonalBand(position, center, band)

        # step_size is smaller than chunk_size (i.e., we are iterating over
        # overlapping intervals) so the last iteration may double count
        # trailing rows.
        if chunk_end == end:
            break


def _nan_to_expectation(
    band: np.ndarray,
    center: int,
    expectation: np.ndarray,
) -> None:
    for i, value in enumerate(band):
        if np.isfinite(value):
            continue
        k = abs(i - center)
        if k < len(expectation) and np.isfinite(expectation[k]):
            band[i] = expectation[k]
        # else: leave the original NaN in place for audit


def _estimate_expectation(
    matrix,
    start: int,
    end: int,
    max_separation: int,
) -> np.ndarray:
    sums = np.zeros(max_separation + 1)
    counts = np.zeros(max_separation + 1)

    for _, center, band in scan_diagonal_bands(matrix, start, end, band_size=max_separation):
        for j, value in enumerate(band):
            k = abs(j - center)
            if np.isfinite(value):
                sums[k] += value
                counts[k] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return sums / counts


def parse_args() -> dict:
    parser = argparse.ArgumentParser(
        description="Dump a single-chromosome contact matrix as a space-separated text file",
    )
    arg = parser.add_argument
    arg("--output", "-o", metavar="FILE", default=None, help="Output file to write the matrix to")
    arg("--region", metavar="...", required=True, help="Genome region to slice ('chr1', 'chr2:0-100000', etc.)")
    arg("--bin-size", "-b", metavar="N", type=int, default=None, help="Use this resolution when the input is an mcool")
    arg("--diagonal", "-d", metavar="N", type=int, default=None, help="Only load the diagonal band of given size")
    arg("--pad-value", metavar="0", default="0", help="Output this value outside of the diagonal band")
    arg("--nan-mode", choices=["zero", "expected"], default=None, help="Fill NaNs by zeros or expected P(s) values")
    arg("--no-balance", action="store_true", default=False, help="Do not balance matrix using stored weights")
    arg("cool", metavar="COOL", help="Path to the cool dataset")
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
