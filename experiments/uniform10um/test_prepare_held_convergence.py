#!/usr/bin/env python3
"""Fast, stdlib-only checks for the held-field sensitivity cards."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_held_convergence import VARIANTS, variant_geometry


SCRIPT = Path(__file__).with_name("prepare_held_convergence.py")


class HeldConvergenceCardsTest(unittest.TestCase):
    def test_geometry_preserves_intended_physical_lengths(self):
        base = variant_geometry("base")
        wide = variant_geometry("wide")
        thick = variant_geometry("thick-cpml")
        long_guard = variant_geometry("long-guard")
        fine = variant_geometry("fine-z")
        self.assertEqual(base["physical_half_aperture_um"], 210)
        self.assertEqual(wide["physical_half_aperture_um"], 270)
        self.assertEqual(thick["physical_half_aperture_um"], 210)
        self.assertEqual(long_guard["longitudinal_guard_cells"], 768)
        self.assertEqual(base["z_cpml_cells"] * base["dz_um"],
                         fine["z_cpml_cells"] * fine["dz_um"])
        self.assertEqual(base["longitudinal_guard_cells"] * base["dz_um"],
                         fine["longitudinal_guard_cells"] * fine["dz_um"])

    def test_cards_are_held_and_reuse_one_particle_file(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            input_file = root / "particles.h5"
            input_file.touch()
            suite = root / "suite"
            command = [sys.executable, str(SCRIPT), "--input", str(input_file),
                       "--output-dir", str(suite), "--electrons", "10409720",
                       "--variants", "base,wide,thick-cpml,long-guard,fine-z"]
            subprocess.run(command, check=True, capture_output=True, text=True)
            plan = json.loads((suite / "held-convergence-plan.json").read_text())
            self.assertEqual(set(plan["variants"]), set(VARIANTS))
            for name in VARIANTS:
                contents = (suite / name / "input.yaml").read_text()
                self.assertIn("hold_transverse_until_activation: true", contents)
                self.assertIn("self_field_feedback: true", contents)
                self.assertIn("initial_self_field:\n  enabled: true", contents)
                self.assertIn("stop: {mode: after-last-element}", contents)
                self.assertIn("z: 1000010.0", contents)
                self.assertIn(str(input_file.resolve()), contents)
                self.assertIn("after-1m", contents)
            self.assertIn("cells: [128, 128, auto]",
                          (suite / "wide" / "input.yaml").read_text())
            self.assertIn("cells: [20, 20, 12]",
                          (suite / "thick-cpml" / "input.yaml").read_text())
            self.assertIn("upper_guard_cells: 768",
                          (suite / "long-guard" / "input.yaml").read_text())
            fine = (suite / "fine-z" / "input.yaml").read_text()
            self.assertIn("cell_size: [5, 5, 2.5]", fine)
            self.assertIn("cells: [10, 10, 24]", fine)
            self.assertIn("lower_guard_cells: 768", fine)
            self.assertIn("alignment_cells: 512", fine)
            rerun = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn("refusing to replace", rerun.stderr)


if __name__ == "__main__":
    unittest.main()
