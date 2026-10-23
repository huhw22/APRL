#!/usr/bin/env python3
"""Place a shortened retirement detector at the existing causal-guard floor."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import yaml

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--workspace", type=Path, required=True)
p.add_argument("--output-dir", type=Path, required=True)
p.add_argument("--retirement-length", type=float, default=0.3)
args = p.parse_args()
w = args.workspace.resolve()
out = args.output_dir.resolve()
out.mkdir(parents=True, exist_ok=True)
original = w / "results/elegant_4m_analysis_20261024/original"
historical = original / "runs/20261021T070618Z-861230"
plan = json.loads((historical / "input-plan.json").read_text())
card = yaml.safe_load((historical / "run.yaml").read_text())
if not math.isfinite(args.retirement_length) or args.retirement_length <= 0:
    raise ValueError("Retirement length must be positive")
mesh = card["mesh"]
entrance = float(plan["ballistic_reference_z_m"])
causal_gap = float(plan["maximum_particle_gamma"]) * math.hypot(
    mesh["cells"][0] * mesh["cell_size"][0],
    mesh["cells"][1] * mesh["cell_size"][1])
lab_cell = mesh["boost_gamma"] * mesh["cell_size"][2]
# One lab-equivalent cell is mandatory; one extra cell avoids roundoff at
# the strict geometric inequality. No old 7.5 m detector position is retained.
detector = entrance + args.retirement_length + causal_gap + 2 * lab_cell
# Let the established time-window calculation determine the stop. Do not keep
# the earlier arbitrary 50 m lower bound after moving the detector upstream.
stop_floor = detector + lab_cell
policy = dict(retirement_entrance_z_m=entrance,
              retirement_length_m=args.retirement_length,
              retirement_exit_z_m=entrance + args.retirement_length,
              detector_z_m=detector, causal_gap_m=causal_gap,
              lab_equivalent_z_cell_m=lab_cell,
              required_clearance_cells=1, configured_clearance_cells=2,
              stop_floor_z_m=stop_floor,
              protected_minimum_photon_energy_eV=50.0, protected_cycles=2,
              retained_tail_margin_fs=100.0, retained_dynamics_allowance_fs=20.0,
              rationale="Shorten the taper and move the detector to the existing causal guard plus one numerical cell; keep the previous tail-coverage policy.")
policy_file = out / "geometry-policy.json"
if policy_file.exists() or (out / "geometry.json").exists():
    raise ValueError("Refusing to overwrite a geometry plan")
policy_file.write_text(json.dumps(policy, indent=2) + "\n")
command = [
    sys.executable, str(Path(__file__).resolve().with_name("geometry_preflight.py")),
    "--input-file", str(original / "particles.h5"),
    "--historical-plan", str(historical / "input-plan.json"),
    "--output", str(out / "geometry.json"),
    "--retirement-entrance-z", repr(entrance),
    "--retirement-length", repr(args.retirement_length),
    "--detector-z", repr(detector), "--stop-z", repr(stop_floor),
    "--tail-margin-fs", "100", "--dynamics-delay-fs", "20",
]
print(json.dumps(policy, indent=2), flush=True)
subprocess.run(command, check=True)
