#!/usr/bin/env python3
"""Size a compact retirement run from observed exit events and geometric delays."""
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
              quiet_tail_fs=0.0, arbitrary_dynamics_allowance_fs=0.0,
              lower_guard_base_cells=64, upper_guard_base_cells=128,
              guard_alignment_cells=1, mesh_alignment_cells=256,
              endpoint_bracketing_samples=1,
              rationale="Keep original causal/frequency protection and left anchoring. Determine the last required signal from historical evolved exit events, retirement flight, and direct travel to the full physical detector aperture; add one sampling interval, no quiet tail.")
policy_file = out / "geometry-policy.json"
if policy_file.exists() or (out / "geometry.json").exists():
    raise ValueError("Refusing to overwrite a geometry plan")
policy_file.write_text(json.dumps(policy, indent=2) + "\n")
arrival_report = out / "arrival-projection.json"
boundary = card["boundary"]["cells"]
subprocess.run([
    sys.executable, str(Path(__file__).resolve().with_name("observed_arrival.py")),
    "--workspace", str(w), "--output", str(arrival_report),
    "--retirement-exit-z", repr(entrance + args.retirement_length),
    "--detector-z", repr(detector),
    "--physical-half-x", repr((mesh["cells"][0] / 2 - boundary[0]) * mesh["cell_size"][0]),
    "--physical-half-y", repr((mesh["cells"][1] / 2 - boundary[1]) * mesh["cell_size"][1]),
    "--dx", repr(mesh["cell_size"][0]), "--dy", repr(mesh["cell_size"][1]),
], check=True)
command = [
    sys.executable, str(Path(__file__).resolve().with_name("geometry_preflight.py")),
    "--input-file", str(original / "particles.h5"),
    "--historical-plan", str(historical / "input-plan.json"),
    "--output", str(out / "geometry.json"),
    "--retirement-entrance-z", repr(entrance),
    "--retirement-length", repr(args.retirement_length),
    "--detector-z", repr(detector), "--stop-z", repr(stop_floor),
    "--tail-margin-fs", "0", "--arrival-report", str(arrival_report),
    "--minimum-lower-guard-cells", "64", "--upper-guard-cells", "128",
    "--guard-alignment-cells", "1", "--mesh-alignment-cells", "256",
    "--stop-rounding-m", "0", "--endpoint-samples", "1",
]
print(json.dumps(policy, indent=2), flush=True)
subprocess.run(command, check=True)
