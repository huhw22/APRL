"""Standalone checks for the streaming Elegant HDF5 v4 cropper.

Run with: python3 -m unittest verification/test_crop_elegant_hdf5.py
Requires numpy and h5py, the same Python dependencies as APRL prepare_case.py.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


TOOL = Path(__file__).resolve().parents[1] / "utilities" / "crop_elegant_hdf5.py"
LIGHT_SPEED_M_S = 299_792_458.0


class CropElegantHdf5Tests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "full.h5"
        self.output = self.root / "center.h5"
        self.dtype = np.dtype([
            ("plane_position_m", "<f8", (2,)),
            ("arrival_time_offset_s", "<f8"),
            ("proper_velocity", "<f8", (3,)),
            ("source_id", "<u8"),
            ("macro_weight", "<f8"),
        ])

    def write_source(self, z_um, weights=None, version=4, beta=0.8):
        z_um = np.asarray(z_um, dtype=float)
        records = np.zeros(len(z_um), dtype=self.dtype)
        records["plane_position_m"][:, 0] = np.arange(len(z_um)) * 1e-6
        records["plane_position_m"][:, 1] = np.arange(len(z_um)) * -1e-6
        records["proper_velocity"][:, 2] = beta / np.sqrt(1 - beta * beta)
        records["arrival_time_offset_s"] = -z_um * 1e-6 / (beta * LIGHT_SPEED_M_S)
        records["source_id"] = np.arange(101, 101 + len(z_um))
        records["macro_weight"] = 1.0 if weights is None else weights
        with h5py.File(self.source, "w") as file:
            group = file.create_group("particles")
            group.attrs["format_version"] = version
            group.attrs["arrival_time_reference"] = "macro-weighted-mean-of-selected-sdds-page"
            group.attrs["source_arrival_time_reference_s"] = 1.23e-9
            group.create_dataset("records", data=records)
        return records

    def run_tool(self, *extra):
        return subprocess.run(
            [sys.executable, str(TOOL), "--input", str(self.source),
             "--output", str(self.output), "--block-records", "2", *extra],
            capture_output=True, text=True, check=False)

    def test_window_preserves_complete_records_and_reference(self):
        source = self.write_source([-40, -30, -10, 0, 10, 30, 40],
                                   [1, 2, 3, 4, 5, 6, 7])
        result = self.run_tool("--width-um", "60")
        self.assertEqual(result.returncode, 0, result.stderr)
        with h5py.File(self.output, "r") as file:
            group = file["particles"]
            np.testing.assert_array_equal(group["records"][:], source[1:6])
            self.assertEqual(group.attrs["format_version"], 4)
            self.assertEqual(group.attrs["source_arrival_time_reference_s"], 1.23e-9)
            self.assertEqual(group.attrs["crop_source_records"], 7)
            self.assertEqual(group.attrs["crop_selected_records"], 5)
            self.assertEqual(group.attrs["crop_source_weight_sum"], 28)
            self.assertEqual(group.attrs["crop_selected_weight_sum"], 20)

    def test_center_shift(self):
        source = self.write_source([-25, 0, 25])
        # An offset of +10 um/c moves the centre toward negative z.
        center_fs = 10e-6 / (0.8 * LIGHT_SPEED_M_S) * 1e15
        result = self.run_tool("--width-um", "30", "--center-time-offset-fs", str(center_fs))
        self.assertEqual(result.returncode, 0, result.stderr)
        with h5py.File(self.output, "r") as file:
            np.testing.assert_array_equal(file["particles/records"][:], source[0:2])

    def test_uses_each_particle_longitudinal_speed(self):
        source = self.write_source([-31, -29, 0])
        slower_beta = 0.6
        source["proper_velocity"][0, 2] = slower_beta / np.sqrt(1 - slower_beta**2)
        with h5py.File(self.source, "r+") as file:
            file["particles/records"][:] = source
        result = self.run_tool()
        self.assertEqual(result.returncode, 0, result.stderr)
        with h5py.File(self.output, "r") as file:
            np.testing.assert_array_equal(file["particles/records"][:], source)

    def test_rejects_overwrite_and_preserves_source(self):
        source = self.write_source([-10, 0, 10])
        self.assertEqual(self.run_tool().returncode, 0)
        result = self.run_tool()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Refusing to overwrite", result.stderr)
        with h5py.File(self.source, "r") as file:
            np.testing.assert_array_equal(file["particles/records"][:], source)

    def test_invalid_version_or_empty_window_leaves_no_output(self):
        self.write_source([0], version=3)
        self.assertNotEqual(self.run_tool().returncode, 0)
        self.assertFalse(self.output.exists())
        self.write_source([100])
        self.assertNotEqual(self.run_tool().returncode, 0)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.root.glob("*.partial")), [])


if __name__ == "__main__":
    unittest.main()
