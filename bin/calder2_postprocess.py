#!/usr/bin/env python3
"""Post-process a CALDER2 output folder: publish diagnostics, drop intermediates.

CALDER2 can exit 0 while producing nothing or only part of its output
(``try()`` around ``build_comp_table_opt``/``CALDER_sub_domains`` in
``R/CALDER_main.R`` swallows errors and the CLI ignores ``system()``
failures). With ``--keep_intermediate`` CALDER retains ``intermediate_data/``;
this script copies a compact ``logs/`` subset (per-chrom logs, warnings, the
cooler chrom table) so a silent failure leaves evidence. The full
``intermediate_data/`` (pixel dumps + ``.Rds``) is dropped unless
``--save-intermediates`` is given.
"""

import argparse
import logging
import os
import shutil
import sys

logging.basicConfig(level=logging.INFO, format="[calder2_postprocess] %(levelname)s: %(message)s")
LOG = logging.getLogger("calder2_postprocess")


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
    # --cool/--resolution identify the cooler; consumed by the chromosome-name
    # remap step. Accepted here so the caller has one stable interface.
    parser.add_argument("--cool", default="", help="Cooler/mcool used for the analysis")
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

    os.makedirs(logs_dir, exist_ok=True)

    if os.path.isdir(intermediate_dir):
        copy_diagnostics(intermediate_dir, logs_dir)
        if args.save_intermediates:
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
