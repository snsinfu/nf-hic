#!/usr/bin/env python3
"""Assemble per-chromosome OnTAD outputs into one TSV.

The bin table is read directly from the (m)cooler, so no separate dump_bins
step is needed.

OnTAD output rows are headerless, 1-based, closed intra-chromosomal bin
intervals: start_bin_id, end_bin_id, tad_level, tad_mean, tad_score. The
output keeps those column names but reports start_bin_id/end_bin_id as the
global bin ids of the source cool (so they can be joined against cooler bins).
`start`/`end` are the bin midpoints; full bin intervals are not repeated
because they can be recovered from the bin ids.
"""

import argparse
import json
import logging
import os
import sys

import cooler
import pandas as pd


LOG = logging.getLogger(__name__)
LOG_LEVEL = logging.INFO
LOG_FORMAT = "[%(asctime)s] %(message)s"
LOG_DATE_FORMAT = "%F %T"

ONTAD_COLUMNS = ["start_bin_id", "end_bin_id", "tad_level", "tad_mean", "tad_score"]
OUTPUT_COLUMNS = [
    "chrom",
    "start",
    "end",
    "start_bin_id",
    "end_bin_id",
    "tad_level",
    "tad_mean",
    "tad_score",
]


def main(
    *,
    cool: str,
    bin_size: int | None,
    output: str,
    prefix: str | None,
    metadata: str | None,
    tads: list[str],
) -> None:
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    if bin_size is not None:
        cool += f"::/resolutions/{bin_size}"

    clr = cooler.Cooler(cool)
    bins = clr.bins()[:][["chrom", "start", "end"]].reset_index(drop=True)
    bins["bin_index"] = bins.index

    frames = []
    for tad_path in tads:
        if not os.path.exists(tad_path) or os.path.getsize(tad_path) == 0:
            LOG.warning("Skipping empty OnTAD output: %s", tad_path)
            continue

        chrom = _chrom_from_name(os.path.basename(tad_path), prefix)
        section = bins[bins["chrom"] == chrom]
        if section.empty:
            LOG.error("Chromosome '%s' from %s is not in %s", chrom, tad_path, cool)
            sys.exit(1)

        section_tads = pd.read_csv(
            tad_path,
            sep="\t",
            header=None,
            names=ONTAD_COLUMNS,
        )
        if section_tads.empty:
            continue

        start_ids = section_tads["start_bin_id"].astype(int).to_numpy()
        end_ids = section_tads["end_bin_id"].astype(int).to_numpy()
        if start_ids.min() < 1 or end_ids.max() > len(section):
            LOG.error(
                "OnTAD bin ids out of range for '%s' in %s (1..%d)",
                chrom,
                tad_path,
                len(section),
            )
            sys.exit(1)

        opening = section.iloc[start_ids - 1]
        closing = section.iloc[end_ids - 1]
        frames.append(
            pd.DataFrame(
                {
                    "chrom": chrom,
                    "start": ((opening["start"] + opening["end"]) // 2).to_numpy(),
                    "end": ((closing["start"] + closing["end"]) // 2).to_numpy(),
                    "start_bin_id": opening["bin_index"].to_numpy(),
                    "end_bin_id": closing["bin_index"].to_numpy(),
                    "tad_level": section_tads["tad_level"].to_numpy(),
                    "tad_mean": section_tads["tad_mean"].to_numpy(),
                    "tad_score": section_tads["tad_score"].to_numpy(),
                }
            )
        )

    table = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=OUTPUT_COLUMNS)
    table = table[OUTPUT_COLUMNS]
    # deterministic genomic order (per-chromosome files arrive in arbitrary order)
    table = table.sort_values(["start_bin_id", "end_bin_id", "tad_level"]).reset_index(drop=True)

    record: dict = {}
    if prefix:
        record["prefix"] = prefix
    if bin_size is not None:
        record["resolution"] = bin_size
    if metadata:
        try:
            extra = json.loads(metadata)
        except json.JSONDecodeError as exc:
            LOG.error("--metadata is not valid JSON: %s", exc)
            sys.exit(1)
        if not isinstance(extra, dict):
            LOG.error("--metadata must be a JSON object")
            sys.exit(1)
        record.update(extra)

    with open(output, "wt") as output_file:
        if record:
            output_file.write("# " + json.dumps(record, sort_keys=True) + "\n")
        table.to_csv(output_file, sep="\t", index=False)


def _chrom_from_name(name: str, prefix: str | None) -> str:
    base = name[:-4] if name.endswith(".tad") else name
    if prefix and base.startswith(prefix + "."):
        return base[len(prefix) + 1 :]
    return base.rsplit(".", 1)[-1]


def parse_args() -> dict:
    parser = argparse.ArgumentParser(description="Assemble OnTAD outputs into a single TSV")
    arg = parser.add_argument
    arg("--cool", metavar="COOL", required=True, help="Path to the cool/mcool dataset")
    arg("--bin-size", "-b", metavar="N", type=int, default=None, help="Resolution when the input is an mcool")
    arg("--output", "-o", metavar="FILE", required=True, help="Output TSV file")
    arg("--prefix", metavar="STR", default=None, help="Prefix stripped from each .tad filename to recover the chromosome")
    arg("--metadata", metavar="JSON", default=None, help="JSON object merged into the single '# {...}' metadata line")
    arg("tads", metavar="FILE", nargs="*", help="OnTAD .tad files")
    return vars(parser.parse_args())


if __name__ == "__main__":
    main(**parse_args())
