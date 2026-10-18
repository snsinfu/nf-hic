#!/usr/bin/env python3
"""Post-process a CALDER2 output folder.

Three jobs:

1. **Chromosome-name remap.** CALDER2 strips a leading lowercase ``chr`` from
   every chromosome (``scripts/calder`` ``sanitize_chroms``) and then
   unconditionally prepends ``chr`` again when writing its BED/TSV outputs
   (``R/CALDER_hierarchy_v2.R``, ``R/build_comp_table_opt.R``,
   ``R/post_processing_nested_domains.R``). The transform is therefore
   ``orig -> "chr" + sub("^chr", "", orig)``. For assemblies whose chromosome
   names have no ``chr`` prefix (S. pombe ``I``/``II``/``III``, Ensembl-style
   ``1``/``2``/``X``/``MT``, ...) the output would be inconsistent with the
   mcool and with the rest of nf-hic (e.g. OnTAD). This script maps the names
   back to the original names present in the cooler.

2. **Diagnostics + intermediate control.** CALDER2 can exit 0 while producing
   nothing or only part of its output (``try()`` around
   ``build_comp_table_opt``/``CALDER_sub_domains`` in ``R/CALDER_main.R``
   swallows errors and the CLI ignores ``system()`` failures). With
   ``--keep_intermediate`` CALDER retains ``intermediate_data/``; this script
   copies a compact ``logs/`` subset (per-chrom logs, warnings, the cooler
   chrom table) so a silent failure leaves evidence. The full
   ``intermediate_data/`` (pixel dumps + ``.Rds``) is dropped unless
   ``--save-intermediates`` is given.

3. **Coordinate normalization.** CALDER2 writes 1-based, fully-closed intervals
   (``pos_start=(bin-1)*size+1``, ``pos_end=bin*size``; nested/domain boundaries
   use ``end_pos*size`` as a point). BED/TSV consumers expect 0-based half-open,
   so the start columns are decremented. The end coordinate is identical under
   both conventions, so it is left unchanged.
"""

import argparse
import logging
import os
import shutil
import sys

import cooler

logging.basicConfig(level=logging.INFO, format="[calder2_postprocess] %(levelname)s: %(message)s")
LOG = logging.getLogger("calder2_postprocess")

# Final (always published) tabular outputs whose first column is the chromosome.
FINAL_CHROM_FILES = [
    "sub_compartments/all_sub_compartments.bed",
    "sub_compartments/all_sub_compartments.tsv",
    "sub_compartments/cor_with_ref.txt",
    "sub_compartments/cor_with_ref.ALL.txt",
    "sub_domains/all_nested_boundaries.bed",
]

# Intermediate per-chromosome interval tables (remapped only when kept).
INTERMEDIATE_CHROM_SUFFIXES = (
    "_sub_compartments.bed",
    "_domain_hierachy.tsv",
    "_domain_boundaries.bed",
    "_nested_boundaries.bed",
)

# CALDER2 writes 1-based, fully-closed coordinates (pos_start=(bin-1)*size+1,
# pos_end=bin*size). BED/TSV consumers expect 0-based half-open. The end
# coordinate is the same under both conventions, so only the "start" columns
# are decremented. Columns are 1-based.
START_COLUMNS_BED_RANGE = (2, 7)   # all_sub_compartments.bed, chr*_sub_compartments.bed (start + thickStart)
START_COLUMNS_BED_POINT = (2, 6)   # all_nested_boundaries.bed, chr*_domain_boundaries.bed, chr*_nested_boundaries.bed
START_COLUMNS_TSV_RANGE = (2,)     # all_sub_compartments.tsv, chr*_domain_hierachy.tsv


def coordinate_start_columns(fname):
    """1-based start columns to make 0-based, chosen by output file name."""
    if fname.endswith(".tsv"):
        return START_COLUMNS_TSV_RANGE
    if fname.endswith(("_domain_boundaries.bed", "_nested_boundaries.bed")):
        return START_COLUMNS_BED_POINT
    if fname.endswith("_sub_compartments.bed"):
        return START_COLUMNS_BED_RANGE
    return ()


def calder_name(orig):
    """Reproduce CALDER2's chromosome-name transform: orig -> chr + strip(orig)."""
    sanitized = orig[3:] if orig.startswith("chr") else orig
    return "chr" + sanitized


def build_name_map(chroms):
    """Map CALDER names back to original cooler names.

    Warn (and leave ambiguous names unmapped) when two distinct original names
    collapse to the same CALDER name, e.g. an mcool with both ``1`` and ``chr1``.
    """
    name_map = {}
    ambiguous = {}
    for orig in chroms:
        key = calder_name(orig)
        if key in name_map and name_map[key] != orig:
            ambiguous.setdefault(key, [name_map[key]]).append(orig)
            name_map.pop(key, None)
        elif key not in ambiguous:
            name_map[key] = orig
    for key, names in ambiguous.items():
        LOG.warning(
            "Ambiguous chromosome name '%s' maps to multiple cooler names %s; "
            "leaving it unchanged in the output.", key, names,
        )
    return name_map


def _as_int(value):
    """Return int(value) for an integer-valued number, else None (headers, text)."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number.is_integer() else None


def rewrite_table(path, name_map, start_columns=()):
    """Rewrite *path* in place: remap CALDER chrom names and convert 1-based
    start columns to 0-based half-open. Non-integer cells (header rows) and
    missing cells are left unchanged."""
    if not os.path.isfile(path):
        return False
    changed = False
    out_lines = []
    with open(path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] in name_map:
                fields[0] = name_map[fields[0]]
                changed = True
            for col in start_columns:
                i = col - 1
                if i < len(fields):
                    number = _as_int(fields[i])
                    if number is not None:
                        fields[i] = str(number - 1)
                        changed = True
            out_lines.append("\t".join(fields))
    if not changed:
        return False
    tmp = path + ".tmp"
    with open(tmp, "w") as handle:
        handle.write("\n".join(out_lines) + ("\n" if out_lines else ""))
    os.replace(tmp, path)
    return True


def remap_intermediate_files(intermediate_dir, name_map):
    for root, _dirs, files in os.walk(intermediate_dir):
        for fname in files:
            if fname.endswith(INTERMEDIATE_CHROM_SUFFIXES):
                rewrite_table(os.path.join(root, fname), name_map,
                              coordinate_start_columns(fname))


def remap_log_filename(fname, name_map):
    """Rename a ``(WARNING_)chr<name>...`` log file to use the original name."""
    prefix = ""
    stem = fname
    if stem.startswith("WARNING_"):
        prefix, stem = "WARNING_", stem[len("WARNING_"):]
    for key in sorted(name_map, key=len, reverse=True):
        if stem == key or stem.startswith(key + "_") or stem.startswith(key + "."):
            return prefix + name_map[key] + stem[len(key):]
    return fname


def copy_diagnostics(intermediate_dir, logs_dir):
    """Copy the small diagnostic files from intermediate_data into logs/."""
    os.makedirs(logs_dir, exist_ok=True)

    chrom_table = os.path.join(intermediate_dir, "chroms.txt")
    if os.path.isfile(chrom_table):
        shutil.copy2(chrom_table, os.path.join(logs_dir, "chroms.txt"))

    for root, _dirs, files in os.walk(intermediate_dir):
        for fname in files:
            if not (fname.endswith("_log.txt") or fname.startswith("WARNING_")):
                continue
            rel = os.path.relpath(root, intermediate_dir)
            dest_dir = logs_dir if rel == "." else os.path.join(logs_dir, rel)
            os.makedirs(dest_dir, exist_ok=True)
            shutil.copy2(os.path.join(root, fname), os.path.join(dest_dir, fname))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cool", required=True, help="Cooler/mcool used for the analysis")
    parser.add_argument("--resolution", default="", help="Resolution group inside the mcool (optional)")
    parser.add_argument("--outdir", required=True, help="CALDER output folder (${prefix})")
    parser.add_argument(
        "--save-intermediates",
        action="store_true",
        help="Keep and publish the full intermediate_data/ folder (large)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = args.outdir
    intermediate_dir = os.path.join(outdir, "intermediate_data")
    logs_dir = os.path.join(outdir, "logs")

    cool_path = args.cool
    if args.resolution:
        cool_path = f"{cool_path}::/resolutions/{args.resolution}"
    try:
        chroms = list(cooler.Cooler(cool_path).chromnames)
    except Exception as exc:  # noqa: BLE001 - diagnostics first, never hide the failure
        LOG.error("Could not read chromosome names from %s: %s", cool_path, exc)
        chroms = []
    name_map = build_name_map(chroms)

    os.makedirs(logs_dir, exist_ok=True)

    for rel in FINAL_CHROM_FILES:
        path = os.path.join(outdir, rel)
        if rewrite_table(path, name_map, coordinate_start_columns(os.path.basename(rel))):
            LOG.info("Rewrote %s (chrom names / 0-based coordinates)", rel)

    if os.path.isdir(intermediate_dir):
        copy_diagnostics(intermediate_dir, logs_dir)
        for root, _dirs, files in os.walk(logs_dir):
            for fname in files:
                new = remap_log_filename(fname, name_map)
                if new != fname:
                    os.replace(os.path.join(root, fname), os.path.join(root, new))
        if args.save_intermediates:
            remap_intermediate_files(intermediate_dir, name_map)
            LOG.info("Keeping full intermediate_data/ (--save-intermediates)")
        else:
            shutil.rmtree(intermediate_dir, ignore_errors=True)
    else:
        LOG.warning(
            "No intermediate_data/ found under %s; CALDER produced no output. "
            "See logs/ (empty if CALDER failed before creating it).", outdir,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
