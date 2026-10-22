#!/usr/bin/env python3
"""Small, self-contained regression tests for compare_particle_planes.py."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SCRIPT = Path(__file__).with_name("compare_particle_planes.py")
RECORD_DTYPE = np.dtype([
    ("source_id", "<u8"),
    ("position_m", "<f8", (3,)),
    ("proper_velocity", "<f8", (3,)),
    ("crossing_direction", "<i4"),
])


def write_detector(path, events, complete=True, plane_z=1.00001):
    records = np.zeros(len(events), dtype=RECORD_DTYPE)
    for index, (source_id, direction, x_shift_um, y_shift_um,
                xp_shift_urad, yp_shift_urad) in enumerate(events):
        records["source_id"][index] = source_id
        records["crossing_direction"][index] = direction
        records["position_m"][index] = (
            (source_id + x_shift_um) * 1e-6,
            (-source_id / 2 + y_shift_um) * 1e-6, plane_z,
        )
        records["proper_velocity"][index] = (
            (2 * source_id + xp_shift_urad) * 1e-6,
            (-source_id + yp_shift_urad) * 1e-6, 1.0,
        )
    with h5py.File(path, "w") as output:
        group = output.create_group("particle_plane")
        group.attrs["plane_z_m"] = plane_z
        group.create_dataset("records", data=records)
        group.create_dataset("committed_records", data=np.uint64(len(records)))
        group.create_dataset("complete", data=np.uint8(complete))


class CompareParticlePlanesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="compare-planes-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.h5"
        self.variant = self.root / "variant.h5"
        self.result = self.root / "result.json"

    def command(self, *extra):
        return subprocess.run([
            sys.executable, str(SCRIPT),
            "--reference", str(self.reference),
            "--variant", str(self.variant),
            "--scratch-dir", str(self.root),
            "--json-out", str(self.result),
            "--chunk", "2", *extra,
        ], capture_output=True, text=True)

    def write_cases(self, complete=True):
        zero = (0, 0, 0, 0)
        shifted = (0.2, -0.1, 1.0, -2.0)
        write_detector(self.reference, [
            (3, 1, *zero), (1, 1, *zero), (2, 1, *zero),
            (3, 1, 100, 0, 0, 0),  # The second downstream crossing is ignored.
            (2, -1, *zero), (5, 1, *zero),
        ])
        write_detector(self.variant, [
            (2, 1, *shifted), (3, 1, *shifted), (1, 1, *shifted),
            (4, 1, *shifted), (3, 1, 100, 0, 0, 0),
            (1, -1, *zero),
        ], complete=complete)

    def test_matches_by_source_id_across_chunks_and_counts_records(self):
        self.write_cases()
        result = self.command()
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(self.result.read_text(encoding="utf-8"))
        self.assertEqual(report["matching"]["matched_particles"], 3)
        self.assertEqual(report["matching"]["reference_only_particles"], 1)
        self.assertEqual(report["matching"]["variant_only_particles"], 1)
        self.assertFalse(report["matching"]["same_unique_population"])
        self.assertEqual(report["reference"]["duplicate_downstream_records"], 1)
        self.assertEqual(report["variant"]["duplicate_downstream_records"], 1)
        self.assertEqual(report["reference"]["upstream_records"], 1)
        self.assertEqual(report["variant"]["upstream_records"], 1)
        self.assertEqual(report["reference"]["unique_downstream_particles"], 4)
        self.assertEqual(report["variant"]["unique_downstream_particles"], 4)
        self.assertAlmostEqual(report["metrics"]["reference"]["x_um"]["mean"], 2)
        self.assertAlmostEqual(report["metrics"]["variant"]["x_um"]["mean"], 2.2)
        for name, expected in (("dx_um", 0.2), ("dy_um", -0.1),
                               ("dxp_urad", 1.0), ("dyp_urad", -2.0)):
            value = report["metrics"]["difference"][name]
            self.assertAlmostEqual(value["mean"], expected)
            self.assertAlmostEqual(value["rms"], abs(expected))
            self.assertAlmostEqual(value["abs_p99"], abs(expected))
        self.assertAlmostEqual(
            report["metrics"]["difference"]["dxy_um"]["mean"],
            np.hypot(0.2, 0.1),
        )

    def test_rejects_incomplete_detector_by_default(self):
        self.write_cases(complete=False)
        result = self.command()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("incomplete detector", result.stderr)
        self.assertFalse(self.result.exists())
        allowed = self.command("--allow-incomplete")
        self.assertEqual(allowed.returncode, 0, allowed.stderr)
        report = json.loads(self.result.read_text(encoding="utf-8"))
        self.assertFalse(report["variant"]["complete"])


if __name__ == "__main__":
    unittest.main()
