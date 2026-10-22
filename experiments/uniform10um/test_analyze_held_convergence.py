#!/usr/bin/env python3
"""Tiny synthetic-HDF5 checks for held-start field convergence analysis."""

import importlib.util
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SPEC = importlib.util.spec_from_file_location(
    "analyze_held_convergence", Path(__file__).with_name("analyze_held_convergence.py"))
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def make_run(root, name, times, x_first, complete=True, reverse_flux=False):
    run = root / name
    (run / "detectors").mkdir(parents=True)
    (run / "ledger").mkdir()
    nx = 6 if x_first == -2.5e-6 else 4
    ny = 4
    dtype = np.dtype([("electric_V_per_m", "<f8", (3,)),
                      ("magnetic_T", "<f8", (3,))])
    fields = np.zeros((len(times), ny, nx), dtype=dtype)
    for i, time in enumerate(times):
        fields["electric_V_per_m"][i, ..., 0] = 1.0 + time
        sign = -1.0 if reverse_flux else 1.0
        fields["magnetic_T"][i, ..., 1] = sign * (1.0 + time) / analysis.C
    with h5py.File(run / "detectors" / "near-field.h5", "w") as file:
        group = file.create_group("field_plane")
        group.create_dataset("fields", data=fields)
        group.create_dataset("time_s", data=np.asarray(times, dtype=np.float64))
        group.create_dataset("committed_samples", data=np.uint64(len(times)))
        group.create_dataset("complete", data=np.uint8(complete))
        group.attrs.update({"plane_z_m": 1.0, "x_first_m": x_first,
                            "y_first_m": -1.5e-6, "dx_m": 1e-6,
                            "dy_m": 1e-6, "nx": nx, "ny": ny})
    ledger_dtype = np.dtype([("time_box_s", "<f8"),
                             ("field_energy_interior_J", "<f8"),
                             ("outward_field_energy_J", "<f8")])
    records = np.zeros(3, dtype=ledger_dtype)
    records["time_box_s"] = [0.0, 1.0, 2.0]
    records["field_energy_interior_J"] = [1.0, 1.5, 2.0]
    records["outward_field_energy_J"] = [0.0, 0.2, 0.4]
    with h5py.File(run / "ledger" / "energy-ledger.h5", "w") as file:
        group = file.create_group("energy_ledger")
        group.create_dataset("records", data=records)
        group.create_dataset("committed_records", data=np.uint64(3))
        group.create_dataset("complete", data=np.uint8(1))
    return run


class HeldConvergenceTest(unittest.TestCase):
    def test_time_interpolation_and_wide_aperture(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = make_run(root, "ref", [0.0, 1.0, 2.0], -1.5e-6)
            wide = make_run(root, "wide", [0.0, 0.5, 1.5, 2.0], -2.5e-6)
            report, trace = analysis.analyze(("ref", ref), [("wide", wide)],
                                             roi_half_width_um=1.0)
            pair = report["comparisons"]["wide"]
            self.assertAlmostEqual(pair["relative_E_spacetime_L2"], 0.0, places=14)
            self.assertAlmostEqual(pair["relative_B_spacetime_L2"], 0.0, places=14)
            self.assertAlmostEqual(pair["matched_forward_energy_difference_J"],
                                   0.0, places=25)
            self.assertEqual(report["fixed_roi"]["nx"], 2)
            self.assertEqual(report["common_lab_window"]["reference_samples"], 3)
            self.assertEqual(sum(row["kind"] == "lab_plane" for row in trace), 6)
            self.assertTrue(pair["tail_unresolved"])
            self.assertEqual(pair["transmitted_energy_status"],
                             "inconclusive_tail_unresolved")
            self.assertAlmostEqual(pair["common_box_field_energy_difference_J"], 0.0)

    def test_misaligned_grid_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = make_run(root, "ref", [0.0, 1.0], -1.5e-6)
            shifted = make_run(root, "shifted", [0.0, 1.0], -1.4e-6)
            with self.assertRaisesRegex(ValueError, "cell centres do not align"):
                analysis.analyze(("ref", ref), [("shifted", shifted)],
                                 roi_half_width_um=1.0)

    def test_cli_json_and_csv(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = make_run(root, "ref", [0.0, 1.0, 2.0], -1.5e-6)
            wide = make_run(root, "wide", [0.0, 0.5, 1.5, 2.0], -2.5e-6)
            json_out = root / "report.json"
            csv_out = root / "trace.csv"
            subprocess.run([
                sys.executable, str(Path(__file__).with_name("analyze_held_convergence.py")),
                "--reference", f"ref={ref}", "--variant", f"wide={wide}",
                "--roi-half-width-um", "1", "--json-out", str(json_out),
                "--csv-out", str(csv_out)], check=True, capture_output=True, text=True)
            report = json.loads(json_out.read_text())
            with csv_out.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(report["reference"], "ref")
            self.assertEqual(len(rows), 12)
            self.assertEqual({row["kind"] for row in rows},
                             {"lab_plane", "boosted_ledger"})

    def test_reverse_flux_and_incomplete_file(self):
        e = np.zeros((1, 1, 3))
        b = np.zeros_like(e)
        e[..., 0] = 2.0
        b[..., 1] = -3.0
        signed, forward = analysis._power(e, b, 1.0)
        self.assertLess(signed, 0)
        self.assertEqual(forward, 0.0)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ref = make_run(root, "ref", [0.0, 1.0], -1.5e-6)
            bad = make_run(root, "bad", [0.0, 1.0], -1.5e-6,
                           complete=False)
            with self.assertRaisesRegex(ValueError, "incomplete"):
                analysis.analyze(("ref", ref), [("bad", bad)],
                                 roi_half_width_um=1.0)


if __name__ == "__main__":
    unittest.main()
