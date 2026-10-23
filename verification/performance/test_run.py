#!/usr/bin/env python3
"""Reject corrupted or numerically changed outputs in the performance harness."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

spec = importlib.util.spec_from_file_location("performance_run", Path(__file__).with_name("run.py"))
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)


def fixture(root):
    (root / "ledger").mkdir()
    (root / "detectors").mkdir()
    (root / "run.log").write_text("Configured stop reached: all valid particles passed\nElapsed wall time [s]: 1.25\n")
    ledger = np.zeros(2, dtype=[("step", "u8"), ("removed_macroparticles", "u8"), ("time_box_s", "f8")])
    ledger["step"] = [0, 1]
    ledger["time_box_s"] = [0., 1e-15]
    particles = np.zeros(2, dtype=[("source_id", "u8"), ("position_m", "f8", (3,))])
    particles["source_id"] = [9, 3]
    particles["position_m"][:, 2] = 3e-6
    fields = np.zeros((2, 2, 2), dtype=[("electric_V_per_m", "f8", (3,))])
    for product, values in zip(run.PRODUCTS, (ledger, particles, fields)):
        filename, group_name, data_name, commit = product
        with h5py.File(root / filename, "w") as file:
            group = file.create_group(group_name)
            group.create_dataset(data_name, data=values)
            group.create_dataset("complete", data=np.uint8(1))
            group.create_dataset(commit, data=np.uint64(len(values)))
            if group_name == "field_plane":
                group.create_dataset("time_s", data=[0., 1e-15])


class PerformanceChecks(unittest.TestCase):
    def test_same_rank_bitwise_rejects_even_tolerated_change(self):
        first = {"value": np.array([1.0])}
        second = {"value": np.array([np.nextafter(1., 2.)])}
        self.assertTrue(run.compare(first, second, 1e-9, 1e-12)["pass"])
        self.assertFalse(run.compare(first, second, 1e-9, 1e-12, True)["pass"])

    def test_zero_field_and_integer_changes_fail(self):
        self.assertFalse(run.compare({"E": np.zeros(1)}, {"E": np.array([1e-99])}, 1e-9, 1e-12)["pass"])
        self.assertFalse(run.compare({"n": np.array([3])}, {"n": np.array([4])}, 1e-9, 1e-12)["pass"])

    def test_complete_output_and_source_id_sort(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture(root)
            summary, arrays = run.inspect_run(root, 2)
            self.assertEqual(summary["steps"], 1)
            np.testing.assert_array_equal(arrays["particle_plane.source_id"], [3, 9])

    def test_incomplete_output_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture(root)
            with h5py.File(root / "detectors/field.h5", "r+") as file:
                file["field_plane/complete"][()] = 0
            with self.assertRaisesRegex(ValueError, "completion"):
                run.inspect_run(root, 2)

    def test_missing_particle_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture(root)
            with self.assertRaisesRegex(ValueError, "endpoint count"):
                run.inspect_run(root, 3)

    def test_nonfinite_output_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture(root)
            with h5py.File(root / "ledger/energy-ledger.h5", "r+") as file:
                values = file["energy_ledger/records"][:]
                values["time_box_s"][1] = np.nan
                file["energy_ledger/records"][:] = values
            with self.assertRaisesRegex(ValueError, "non-finite"):
                run.inspect_run(root, 2)


if __name__ == "__main__":
    unittest.main()
