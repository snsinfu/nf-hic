#!/usr/bin/env python3
"""Regression test: OnTAD start/end are opening/closing bin midpoints.

Run: python3 tests/test_ontad_assemble.py
"""
import importlib.util
import os
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.normpath(os.path.join(HERE, "..", "bin", "ontad_assemble.py"))


class FakeArray(list):
    def __sub__(self, value):
        return FakeArray(item - value for item in self)

    def min(self):
        return min(self)

    def max(self):
        return max(self)


class FakeColumn:
    def __init__(self, values):
        self.values = list(values)

    def astype(self, dtype):
        if dtype in (int, "int", "int64"):
            return FakeColumn(int(value) for value in self.values)
        return FakeColumn(self.values)

    def to_numpy(self):
        return FakeArray(self.values)

    def __add__(self, other):
        if isinstance(other, FakeColumn):
            return FakeColumn(a + b for a, b in zip(self.values, other.values))
        return FakeColumn(value + other for value in self.values)

    def __floordiv__(self, other):
        return FakeColumn(value // other for value in self.values)

    def __eq__(self, other):
        return [value == other for value in self.values]


class FakeFrame:
    def __init__(self, data):
        items = list(data.items())
        length = next(
            (len(value) for _, value in items if isinstance(value, (FakeColumn, list, tuple, FakeArray))),
            0,
        )
        self._data = {}
        for name, value in items:
            if isinstance(value, FakeColumn):
                self._data[name] = list(value.values)
            elif isinstance(value, (list, tuple, FakeArray)):
                self._data[name] = list(value)
            else:
                self._data[name] = [value] * length
        self.index = list(range(length))

    def __len__(self):
        return len(self.index)

    @property
    def empty(self):
        return len(self) == 0

    def __getitem__(self, key):
        if isinstance(key, list):
            if key and all(isinstance(value, bool) for value in key):
                return FakeFrame({
                    name: [value for value, keep in zip(values, key) if keep]
                    for name, values in self._data.items()
                })
            return FakeFrame({name: self._data[name] for name in key})
        return FakeColumn(self._data[key])

    def __setitem__(self, key, value):
        if isinstance(value, FakeColumn):
            self._data[key] = list(value.values)
        elif isinstance(value, (list, tuple, FakeArray)):
            self._data[key] = list(value)
        else:
            self._data[key] = [value] * len(self)

    @property
    def iloc(self):
        return _Iloc(self)

    def reset_index(self, drop=True):
        return FakeFrame(self._data)

    def sort_values(self, by):
        columns = [by] if isinstance(by, str) else list(by)
        order = sorted(
            range(len(self)),
            key=lambda index: tuple(self._data[column][index] for column in columns),
        )
        return FakeFrame({
            name: [values[index] for index in order]
            for name, values in self._data.items()
        })

    def to_csv(self, output, sep="\t", index=False):
        columns = list(self._data)
        output.write(sep.join(columns) + "\n")
        for row in range(len(self)):
            output.write(sep.join(str(self._data[column][row]) for column in columns) + "\n")


class _Iloc:
    def __init__(self, frame):
        self.frame = frame

    def __getitem__(self, positions):
        if isinstance(positions, int):
            positions = [positions]
        positions = list(positions)
        return FakeFrame({
            name: [values[position] for position in positions]
            for name, values in self.frame._data.items()
        })


def _fake_pandas():
    module = types.ModuleType("pandas")

    def read_csv(path, sep="\t", header=None, names=None):
        with open(path) as handle:
            rows = [line.rstrip("\n").split(sep) for line in handle if line.strip()]
        return FakeFrame({
            name: [row[index] for row in rows]
            for index, name in enumerate(names)
        })

    def concat(frames, ignore_index=True):
        columns = list(frames[0]._data)
        data = {column: [] for column in columns}
        for frame in frames:
            for column in columns:
                data[column].extend(frame._data[column])
        return FakeFrame(data)

    module.DataFrame = FakeFrame
    module.read_csv = read_csv
    module.concat = concat
    return module


def _fake_cooler():
    module = types.ModuleType("cooler")

    class Bins:
        def __getitem__(self, key):
            return FakeFrame({
                "chrom": ["chr1", "chr1", "chr1"],
                "start": [0, 10000, 20000],
                "end": [10000, 20000, 30000],
            })

    class Cooler:
        def __init__(self, path):
            pass

        def bins(self):
            return Bins()

    module.Cooler = Cooler
    return module


def _load():
    sys.modules["pandas"] = _fake_pandas()
    sys.modules["cooler"] = _fake_cooler()
    spec = importlib.util.spec_from_file_location("ontad_assemble", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MidpointCoordinatesTest(unittest.TestCase):
    def test_start_and_end_are_opening_and_closing_bin_midpoints(self):
        with tempfile.TemporaryDirectory() as tmp:
            tad_path = os.path.join(tmp, "chr1.tad")
            with open(tad_path, "w") as handle:
                handle.write("1\t2\t0\t1.0\t1.0\n")
            output_path = os.path.join(tmp, "out.tsv")

            module = _load()
            module.main(
                cool="fake.cool",
                bin_size=10000,
                output=output_path,
                prefix=None,
                metadata=None,
                tads=[tad_path],
            )

            with open(output_path) as handle:
                lines = [
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if not line.startswith("#")
                ]
            record = dict(zip(lines[0], lines[1]))
            self.assertEqual(record["start"], "5000")
            self.assertEqual(record["end"], "15000")
            self.assertEqual(record["start_bin_id"], "0")
            self.assertEqual(record["end_bin_id"], "1")


if __name__ == "__main__":
    unittest.main()
