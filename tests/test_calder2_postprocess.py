#!/usr/bin/env python3
"""Regression test: CALDER2 outputs are rewritten to 0-based half-open BED.

Run: python3 tests/test_calder2_postprocess.py
"""
import importlib.util
import os
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.normpath(os.path.join(HERE, "..", "bin", "calder2_postprocess.py"))


def _load():
    # The module imports `cooler`; inject a stub so the test needs no cooler.
    fake = types.ModuleType("cooler")

    class Cooler:
        def __init__(self, path):
            self.chromnames = ["chr1"]

    fake.Cooler = Cooler
    sys.modules["cooler"] = fake
    spec = importlib.util.spec_from_file_location("calder2_postprocess", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class NormalizeCoordinatesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.outdir = self.tmp.name
        os.makedirs(os.path.join(self.outdir, "sub_compartments"))
        os.makedirs(os.path.join(self.outdir, "sub_domains"))
        self._write("sub_compartments/all_sub_compartments.bed",
                    "chr1\t10001\t20000\tA.1\t1\t.\t10001\t20000\t#FF0000\n")
        self._write("sub_compartments/all_sub_compartments.tsv",
                    "chr\tpos_start\tpos_end\tcomp_name\tcomp_rank\tcontinous_rank\n"
                    "chr1\t10001\t20000\tA.1\t1\t1\n")
        self._write("sub_domains/all_nested_boundaries.bed",
                    "chr1\t20000\t20000\t\t.\t20000\t20000\t#000000\n")

    def _write(self, rel, text):
        with open(os.path.join(self.outdir, rel), "w") as handle:
            handle.write(text)

    def _run(self):
        mod = _load()
        argv = sys.argv
        sys.argv = ["calder2_postprocess.py", "--cool", "fake.cool", "--outdir", self.outdir]
        try:
            mod.main()
        finally:
            sys.argv = argv

    def test_coordinates_are_zero_based_half_open(self):
        self._run()
        with open(os.path.join(self.outdir, "sub_compartments/all_sub_compartments.bed")) as fh:
            bed = fh.readline().rstrip("\n").split("\t")
        self.assertEqual(bed[1], "10000")   # start: 0-based
        self.assertEqual(bed[2], "20000")   # end: unchanged
        self.assertEqual(bed[6], "10000")   # thickStart
        self.assertEqual(bed[7], "20000")   # thickEnd

        with open(os.path.join(self.outdir, "sub_compartments/all_sub_compartments.tsv")) as fh:
            lines = [line.rstrip("\n").split("\t") for line in fh]
        self.assertEqual(lines[0][1], "pos_start")  # header untouched
        self.assertEqual(lines[1][1], "10000")
        self.assertEqual(lines[1][2], "20000")

        with open(os.path.join(self.outdir, "sub_domains/all_nested_boundaries.bed")) as fh:
            point = fh.readline().rstrip("\n").split("\t")
        self.assertEqual(point[1], "19999")  # point p -> [p-1, p)
        self.assertEqual(point[2], "20000")


if __name__ == "__main__":
    unittest.main()
