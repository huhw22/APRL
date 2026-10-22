#!/usr/bin/env python3
"""Write paired APRL cards for the independent 10 um / 50 A HDF5-v4 bunch.

Cards are written only; this utility never launches APRL.  A single input
file and identical mesh are shared by all selected comparisons.
"""

import argparse
import json
from pathlib import Path

PROFILES = {
    # profile: (boost gamma, nz, dz [um], duration [fs], stop z [um],
    #           field-plane rhythm [fs], longitudinal CPML cells)
    "g2-short": (2.0, 256, 1.0, 180.0, 30.0, 5.0, 12),
    # Same 256 um box and 12 um z-CPML thickness as g2-short, at half dz.
    "g2-short-fine": (2.0, 512, 0.5, 180.0, 30.0, 5.0, 24),
    "g2": (2.0, 256, 1.0, 300.0, 100.0, 5.0, 12),
    "g100": (100.0, 2048, 2.0, 7000.0, 200000.0, 25.0, 12),
    # The observation plane is exactly 1 m beyond the initial center at 10 um.
    # The whole 10 um lab bunch takes about 34 ps / 2040 Maxwell steps to pass.
    "g929-1m": (929.130549, "auto", 5.0, 50000.0, 1000010.0, 1.0, 12),
}

CASES = {
    "ballistic-field": (False, True, True),
    "held-field": (True, True, True),
    "ballistic-zero-field": (False, False, True),
    "held-zero-field": (True, False, True),
    "ballistic-no-feedback": (False, True, False),
    "held-no-feedback": (True, True, False),
}


def yaml_bool(value):
    return "true" if value else "false"


def card(profile, name, input_file, electrons, output_root,
         particle_plane, *, xy_cells=104, xy_cpml_cells=10,
         dz_um=None, z_cpml_cells=None, longitudinal_guard_cells=384,
         longitudinal_alignment_cells=256):
    gamma, nz, profile_dz, duration, stop_z, field_rhythm, profile_cpml_z = PROFILES[profile]
    dz = profile_dz if dz_um is None else dz_um
    cpml_z = profile_cpml_z if z_cpml_cells is None else z_cpml_cells
    held, initial_field, feedback = CASES[name]
    one_meter = profile == "g929-1m"
    run_dir = (output_root / "runs" / (profile + "-" + name)).resolve()
    input_file = input_file.resolve()
    plane_z = stop_z if one_meter else 20.0
    particle_name = "after-1m" if one_meter else "after-20um"
    particle_entry = (
        f"    - {{name: {particle_name}, z: {plane_z}, "
        "buffer_records: 65536, compression: 0}\n"
    ) if particle_plane else ""
    particle_plane_yaml = particle_entry or "    []\n"
    longitudinal_yaml = f"""  longitudinal:
    mode: observation-bounded
    lower_guard_cells: {longitudinal_guard_cells}
    upper_guard_cells: {longitudinal_guard_cells}
    self_field_padding_rms: 4.0
    ballistic_time_safety_factor: 1.1
    alignment_cells: {longitudinal_alignment_cells}
""" if one_meter else ""
    stop_yaml = ("stop: {mode: after-last-element}" if one_meter else
                 f"stop: {{mode: reference-center-z, z: {stop_z}}}")
    runtime_yaml = (
        "runtime: {mode: throughput, resource_monitor: {enabled: true, "
        "progress_interval_steps: 50, calibration_steps: 1}}" if one_meter
        else "runtime: {mode: throughput, resource_monitor: {enabled: false}}"
    )
    return f"""# Independent synthetic 10 um, 50 A laboratory-plane event bunch.
# Profile {profile}; case {name}. No Elegant data and no magnetic elements.
units: {{length: micrometer, time: femtosecond}}
{runtime_yaml}
output:
  overwrite: false
  manifest: {json.dumps(str(run_dir / 'run-manifest.yaml'))}
mesh:
  field_solver: cowan-z
  current_filter: none
  ballistic_outside_transverse: true
  self_field_feedback: {yaml_bool(feedback)}
  cells: [{xy_cells}, {xy_cells}, {nz}]
  cell_size: [5, 5, {dz}]
  center: [0, 0, 0]
  duration: {duration}
  boost_gamma: {gamma}
  particle_steps_per_undulator_period: 1
{longitudinal_yaml}boundary:
  type: cpml
  cells: [{xy_cpml_cells}, {xy_cpml_cells}, {cpml_z}]
  polynomial_order: 3
  target_reflection: 1.0e-8
  kappa_max: 8
  alpha_fraction: 0
initial_self_field:
  enabled: {yaml_bool(initial_field)}
  model: relativistic-poisson
  include_transverse_cpml_particles: {yaml_bool(initial_field)}
  relative_tolerance: 1.0e-7
  maximum_mean_transverse_beta: 0.01
  maximum_iterations: 10000
beam:
  reference:
    input_plane_z: 0
    initial_center_z: 10
    frozen_lab_plane_injection: true
    activation_plane_z: 0
    hold_transverse_until_activation: {yaml_bool(held)}
  input:
    type: hdf5
    file: {json.dumps(str(input_file))}
    electrons: {electrons:.17g}
    position_offset: [0, 0, 0]
sources: {{incident_waves: [], magnetic_elements: []}}
detectors:
  enabled: true
  directory: {json.dumps(str(run_dir / 'detectors'))}
  field_planes:
    - {{name: near-field, z: {plane_z}, rhythm: {field_rhythm}, buffer_samples: 4, compression: 0}}
  particle_planes:
{particle_plane_yaml}{stop_yaml}
trajectory: {{enabled: false}}
energy_ledger:
  enabled: true
  directory: {json.dumps(str(run_dir / 'ledger'))}
  filename: energy-ledger.h5
  sample_interval_steps: 1
  buffer_records: 256
  compression: 0
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="existing synthetic HDF5-v4 file")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--profile", choices=PROFILES, default="g2-short")
    parser.add_argument("--electrons", type=float, default=10409720.0,
                        help="represented physical electrons from generator metadata")
    parser.add_argument("--cases", default=",".join(CASES),
                        help="comma-separated names; default writes six comparisons")
    parser.add_argument("--particle-plane", action="store_true",
                        help="also write all particles at the profile's observation plane (roughly 112 bytes/row/case)")
    parser.add_argument("--overwrite-cards", action="store_true")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error("input file does not exist")
    if args.electrons <= 0:
        parser.error("--electrons must be positive")
    names = [part.strip() for part in args.cases.split(",") if part.strip()]
    if not names or len(set(names)) != len(names) or any(
            name not in CASES for name in names):
        parser.error("--cases must be unique names from: " + ", ".join(CASES))
    cards_dir = args.output_dir.resolve() / "cards"
    targets = [cards_dir / (args.profile + "-" + name + ".yaml")
               for name in names]
    if not args.overwrite_cards:
        existing = [path for path in targets if path.exists()]
        if existing:
            parser.error("refusing to replace existing card: " + str(existing[0]))
    cards_dir.mkdir(parents=True, exist_ok=True)
    for name, target in zip(names, targets):
        target.write_text(card(args.profile, name, args.input,
                               args.electrons, args.output_dir.resolve(),
                               args.particle_plane), encoding="utf-8")
        print(target)



if __name__ == "__main__":
    main()
