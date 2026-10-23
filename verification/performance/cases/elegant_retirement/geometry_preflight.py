#!/usr/bin/env python3
"""Stream the real v4 input on a Slurm node and size a held/retirement run.

This is a kinematic planning check, not an evolved-orbit or radiation test.
The extra dynamics delay is an explicit allowance. The production manifest,
retirement completion counts, and saved field tail still need checking later.
No particle subset or modified particle file is produced.
"""

import argparse
import json
import math
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

C = 299792458.0
H_EV_S = 4.135667696e-15


def rounded_up(value, multiple):
    return int(math.ceil(value / multiple) * multiple)


def kinematics(records, gamma_boost, beta_boost):
    u = records["proper_velocity"]
    gamma = np.sqrt(1.0 + np.sum(u * u, axis=1))
    transverse_mass2 = 1.0 + u[:, 0] ** 2 + u[:, 1] ** 2
    rapidity = gamma_boost * (1.0 + beta_boost)
    plus = (gamma + u[:, 2]) / rapidity
    minus = transverse_mass2 / plus
    beta_box = (plus - minus) / (plus + minus)
    # Same longitudinal-only laboratory momentum as heldLongitudinalBetaZ().
    held_plus = (np.hypot(1.0, u[:, 2]) + u[:, 2]) / rapidity
    beta_held = (held_plus * held_plus - 1.0) / (held_plus * held_plus + 1.0)
    beta_lab = u[:, 2] / gamma
    # Stable 1/beta_lab - 1 for ultra-relativistic positive u_z.
    delay_factor = transverse_mass2 / (u[:, 2] * (gamma + u[:, 2]))
    return gamma, beta_lab, beta_box, beta_held, delay_factor


def analyze(args):
    input_path = Path(args.input_file).resolve()
    stat_before = input_path.stat()
    plan = json.loads(Path(args.historical_plan).read_text())
    gamma_b = args.boost_gamma
    beta_b = math.sqrt((1.0 - 1.0 / gamma_b) * (1.0 + 1.0 / gamma_b))
    one_minus_beta = 1.0 / (gamma_b * gamma_b * (1.0 + beta_b))
    z_input = float(plan["input_plane_z_m"])
    z_reference = float(plan["initial_center_z_m"])
    common_time = float(plan["common_time_after_reference_s"])
    z_exit = args.retirement_entrance_z + args.retirement_length
    flight_to_exit = z_exit - z_input
    flight_to_detector = args.detector_z - z_input
    if min(flight_to_exit, flight_to_detector) <= 0:
        raise ValueError("Retirement exit and detector must follow the input plane")
    minima = {key: math.inf for key in ("arrival", "gamma", "beta_lab", "beta_box", "beta_held")}
    maxima = {key: -math.inf for key in minima}
    t_origin = math.inf
    latest_exit_photon_retarded = -math.inf
    outside_physical = outside_full = 0
    weight_sum = np.longdouble(0)
    beta_weight_sum = np.longdouble(0)
    beta_time_weight_sum = np.longdouble(0)
    records_seen = 0
    half_full = np.array([args.cells_x * args.dx, args.cells_y * args.dy]) * 0.5
    half_physical = half_full - args.transverse_cpml_cells * np.array([args.dx, args.dy])
    if np.any(half_physical <= 0):
        raise ValueError("Transverse CPML leaves no physical aperture")
    with h5py.File(input_path, "r") as handle:
        group = handle["particles"]
        if int(group.attrs.get("format_version", -1)) != 4:
            raise ValueError("This planner requires fixed-lab-plane HDF5 v4 input")
        dataset = group["records"]
        required = {"plane_position_m", "arrival_time_offset_s", "proper_velocity", "macro_weight"}
        if not required.issubset(dataset.dtype.names or ()):
            raise ValueError("Unexpected v4 input schema")
        count = len(dataset)
        if count == 0 or count != int(plan["records"]):
            raise ValueError("Input count does not match the historical plan")
        print(f"geometry pass 1/2: {count} records", flush=True)
        for begin in range(0, count, args.chunk_rows):
            rows = dataset[begin:begin + args.chunk_rows]
            for key in required:
                if not np.all(np.isfinite(rows[key])):
                    raise ValueError(f"Nonfinite {key} at input chunk {begin}")
            if np.any(rows["macro_weight"] <= 0) or np.any(rows["proper_velocity"][:, 2] <= 0):
                raise ValueError("Input requires positive weights and forward laboratory velocity")
            gamma, beta_lab, beta_box, beta_held, delay_factor = kinematics(rows, gamma_b, beta_b)
            arrival = rows["arrival_time_offset_s"] - common_time
            alpha = beta_b + beta_held
            if np.any(alpha <= 0):
                raise ValueError("Held worldlines do not advance in laboratory z")
            # simulation.cpp:1848-1862; input plane == activation plane.
            pretrigger = 2.0 * args.dz / (gamma_b * C * alpha)
            t_origin = min(t_origin, float(np.min(arrival - pretrigger)))
            latest_exit_photon_retarded = max(latest_exit_photon_retarded,
                float(np.max(arrival + flight_to_exit / C * delay_factor)))
            for key, values in (("arrival", arrival), ("gamma", gamma), ("beta_lab", beta_lab),
                                ("beta_box", beta_box), ("beta_held", beta_held)):
                minima[key] = min(minima[key], float(np.min(values)))
                maxima[key] = max(maxima[key], float(np.max(values)))
            xy = rows["plane_position_m"]
            outside_physical += int(np.count_nonzero(np.any(np.abs(xy) > half_physical, axis=1)))
            outside_full += int(np.count_nonzero(np.any(np.abs(xy) > half_full, axis=1)))
            weights = rows["macro_weight"].astype(np.longdouble)
            wide_beta = beta_lab.astype(np.longdouble)
            weight_sum += np.sum(weights, dtype=np.longdouble)
            beta_weight_sum += np.sum(weights * wide_beta, dtype=np.longdouble)
            beta_time_weight_sum += np.sum(weights * wide_beta * rows["arrival_time_offset_s"], dtype=np.longdouble)
            records_seen += len(rows)

        recomputed_common_time = float(((np.longdouble(z_reference) - np.longdouble(z_input)) / C *
            weight_sum + beta_time_weight_sum) / beta_weight_sum)
        common_time_error = recomputed_common_time - common_time
        if abs(common_time_error) > 1e-20:
            raise ValueError(f"Historical common-time mismatch: {common_time_error:.6g} s")
        # The stop worldline is defined by the lab event (t=0, z_reference).
        # Eliminating the held frame origin gives this exact fixed-plane relation.
        target_retarded = latest_exit_photon_retarded + (
            args.dynamics_delay_fs + args.tail_margin_fs) * 1e-15
        minimum_stop_z = z_input + gamma_b**2 * (z_reference - z_input) + (
            gamma_b**2 * beta_b * C * target_retarded) + (
            beta_b / (1.0 + beta_b) * flight_to_detector)
        stop_z = max(args.stop_z, math.ceil(minimum_stop_z / args.stop_rounding_m) * args.stop_rounding_m)
        reference_at_box_zero = z_input + gamma_b**2 * (
            z_reference - z_input + beta_b * C * t_origin)
        stop_box = (stop_z - reference_at_box_zero) / (gamma_b * beta_b * C)
        dt_box = args.dz / C  # Existing cowan-z solver, unchanged mesh spacing.
        step_count = math.ceil(stop_box / dt_box)
        rounded_stop_box = step_count * dt_box
        detector_first_lab = t_origin + beta_b * flight_to_detector / C
        detector_stop_lab = t_origin + rounded_stop_box / gamma_b + beta_b * flight_to_detector / C
        # Work with retarded time to avoid subtracting two ~25 ns timestamps.
        detector_first_retarded = t_origin - one_minus_beta * flight_to_detector / C
        detector_stop_retarded = detector_first_retarded + rounded_stop_box / gamma_b
        detector_z_at_start = flight_to_detector / gamma_b
        detector_z_at_stop = detector_z_at_start - beta_b * C * rounded_stop_box
        ballistic_time = stop_box * args.ballistic_time_safety_factor

        initial_min = ballistic_min = worldtube_min = math.inf
        initial_max = ballistic_max = worldtube_max = -math.inf
        retirement_latest_box = -math.inf
        projected_retirement_entry_outside = projected_retirement_exit_outside = 0
        print("geometry pass 2/2: exact held extents and swept detector window", flush=True)
        for begin in range(0, count, args.chunk_rows):
            rows = dataset[begin:begin + args.chunk_rows]
            _, _, beta_box, beta_held, delay_factor = kinematics(rows, gamma_b, beta_b)
            arrival = rows["arrival_time_offset_s"] - common_time
            activation_box = gamma_b * (arrival - t_origin)
            # eb_bunch.cpp:390-406, before automatic-layout translation.
            z_initial = -C * (beta_b + beta_held) * activation_box
            finish = z_initial + beta_box * C * ballistic_time
            initial_min = min(initial_min, float(np.min(z_initial)))
            initial_max = max(initial_max, float(np.max(z_initial)))
            ballistic_min = min(ballistic_min, float(np.min(z_initial)), float(np.min(finish)))
            ballistic_max = max(ballistic_max, float(np.max(z_initial)), float(np.max(finish)))
            z_activation = -beta_b * C * activation_box
            exit_box = gamma_b * (arrival - t_origin + flight_to_exit / C * (delay_factor + one_minus_beta))
            z_exit_box = flight_to_exit / (gamma_b * (1.0 + beta_b)) - (
                gamma_b * beta_b * C * (arrival - t_origin + flight_to_exit / C * delay_factor))
            retirement_latest_box = max(retirement_latest_box, float(np.max(exit_box)))
            worldtube_min = min(worldtube_min, float(np.min(z_initial)), float(np.min(z_activation)), float(np.min(z_exit_box)))
            worldtube_max = max(worldtube_max, float(np.max(z_initial)), float(np.max(z_activation)), float(np.max(z_exit_box)))
            slopes = rows["proper_velocity"][:, :2] / rows["proper_velocity"][:, 2, None]
            entry_xy = rows["plane_position_m"] + slopes * (args.retirement_entrance_z - z_input)
            exit_xy = rows["plane_position_m"] + slopes * flight_to_exit
            projected_retirement_entry_outside += int(np.count_nonzero(np.any(np.abs(entry_xy) > half_physical, axis=1)))
            projected_retirement_exit_outside += int(np.count_nonzero(np.any(np.abs(exit_xy) > half_physical, axis=1)))

    radiation_upper = initial_max + (stop_z - (z_input + gamma_b * initial_max)) / (
        gamma_b * (1.0 + beta_b))
    content_lower = min(initial_min, ballistic_min)
    content_upper = max(initial_max, ballistic_max, radiation_upper)
    # No self_field_padding_rms; deriveLongitudinalLayout() uses these extrema.
    required_lower_guard = max(args.minimum_lower_guard_cells, math.ceil(
        max(0.0, (content_lower - min(detector_z_at_stop, worldtube_min)) / args.dz)) + args.window_guard_cells)
    lower_guard = rounded_up(required_lower_guard, args.guard_alignment_cells)
    required_upper_guard = max(args.upper_guard_cells, math.ceil(
        max(0.0, (max(detector_z_at_start, worldtube_max) - content_upper) / args.dz)) + args.window_guard_cells)
    upper_guard = rounded_up(required_upper_guard, args.guard_alignment_cells)
    required_lower = content_lower - lower_guard * args.dz
    required_upper = content_upper + upper_guard * args.dz
    minimum_physical_cells = math.ceil((required_upper - required_lower) / args.dz)
    cells_z = rounded_up(minimum_physical_cells + 2 * args.z_cpml_cells, args.mesh_alignment_cells)
    extent_z = cells_z * args.dz
    outer_lower = -extent_z * 0.5
    physical_lower = outer_lower + args.z_cpml_cells * args.dz
    physical_upper = -outer_lower - args.z_cpml_cells * args.dz
    translation = physical_lower - required_lower
    plane_end_margin_cells = (detector_z_at_stop + translation - physical_lower) / args.dz
    plane_start_margin_cells = (physical_upper - detector_z_at_start - translation) / args.dz
    max_gamma = maxima["gamma"]
    delay_factor_gamma = 1.0 / max_gamma**2 / (
        math.sqrt(1.0 - 1.0 / max_gamma**2) * (1.0 + math.sqrt(1.0 - 1.0 / max_gamma**2)))
    protected_cycles = args.retirement_length * delay_factor_gamma * args.minimum_photon_energy_ev / (C * H_EV_S)
    required_taper_length = args.protection_cycles * C * H_EV_S / (
        args.minimum_photon_energy_ev * delay_factor_gamma)
    causal_distance = max_gamma * math.hypot(args.cells_x * args.dx, args.cells_y * args.dy)
    minimum_detector_z = z_exit + causal_distance + gamma_b * args.dz
    tail_margin_actual_fs = (detector_stop_retarded - latest_exit_photon_retarded) * 1e15 - args.dynamics_delay_fs
    latest_retirement_with_allowance = retirement_latest_box + gamma_b * args.dynamics_delay_fs * 1e-15
    stat_after = input_path.stat()
    checks = {
        "input_unchanged": (stat_before.st_size, stat_before.st_mtime_ns) == (stat_after.st_size, stat_after.st_mtime_ns),
        "count_matches_plan": records_seen == int(plan["records"]),
        "initial_all_inside_physical_aperture": outside_physical == 0,
        "detector_causal_clearance": args.detector_z >= minimum_detector_z,
        "protected_cycles_sufficient": protected_cycles >= args.protection_cycles,
        "stop_after_detector": stop_z > args.detector_z + gamma_b * args.dz,
        "full_tail_margin": tail_margin_actual_fs >= args.tail_margin_fs - 1e-6,
        "retirement_complete_before_stop_kinematic_with_allowance": rounded_stop_box > latest_retirement_with_allowance,
        "detector_in_physical_mesh_entire_run": min(plane_end_margin_cells, plane_start_margin_cells) >= args.window_guard_cells - 1e-6,
        "held_ballistic_worldtube_in_physical_mesh": worldtube_min + translation >= physical_lower and worldtube_max + translation <= physical_upper,
        "hard_duration_sufficient": args.duration > rounded_stop_box,
    }
    return {
        "format": "aprl-held-retirement-geometry-v1", "pass": all(checks.values()), "checks": checks,
        "input_file": str(input_path), "input_size_bytes": stat_before.st_size,
        "input_mtime_ns": stat_before.st_mtime_ns, "historical_plan": str(Path(args.historical_plan).resolve()),
        "records": records_seen, "input_weight_sum": float(weight_sum),
        "initial_outside_physical_count": outside_physical, "initial_outside_full_mesh_count": outside_full,
        "retirement_entrance_z_m": args.retirement_entrance_z, "retirement_length_m": args.retirement_length,
        "retirement_exit_z_m": z_exit, "detector_z_m": args.detector_z, "stop_z_m": stop_z,
        "candidate_stop_z_m": args.stop_z, "minimum_stop_z_m": minimum_stop_z,
        "expected_stop_box_s": stop_box, "step_rounded_stop_box_s": rounded_stop_box,
        "expected_maxwell_steps": step_count, "mesh_duration_s": args.duration,
        "recommended_lower_guard_cells": lower_guard, "minimum_lower_guard_cells": required_lower_guard,
        "recommended_upper_guard_cells": upper_guard, "estimated_cells_z": cells_z,
        "recommended_configuration": {
            "boost_gamma": gamma_b, "cells": [args.cells_x, args.cells_y, "auto"],
            "cell_size_m": [args.dx, args.dy, args.dz], "self_field_padding_rms": 0.0,
            "lower_guard_cells": lower_guard, "upper_guard_cells": upper_guard,
            "alignment_cells": args.mesh_alignment_cells, "ballistic_time_safety_factor": args.ballistic_time_safety_factor,
            "frozen_lab_plane_injection": True, "hold_transverse_until_activation": True,
            "activation_plane_z_m": z_input, "ballistic_outside_transverse": False,
            "include_transverse_cpml_particles": False,
        },
        "input_kinematics": {"minima": minima, "maxima": maxima, "common_time_s": common_time,
            "recomputed_common_time_s": recomputed_common_time, "common_time_error_s": common_time_error,
            "input_plane_z_m": z_input, "snapshot_center_z_m": z_reference},
        "held_geometry": {"frame_origin_lab_time_s": t_origin, "frame_origin_lab_z_m": z_input,
            "initial_z_untranslated_m": [initial_min, initial_max],
            "planner_ballistic_z_untranslated_m": [ballistic_min, ballistic_max],
            "held_then_ballistic_worldtube_through_retirement_m": [worldtube_min, worldtube_max],
            "latest_ballistic_retirement_exit_box_s": retirement_latest_box,
            "latest_retirement_exit_box_s_with_allowance": latest_retirement_with_allowance,
            "field_only_box_time_after_retirement_s": rounded_stop_box - latest_retirement_with_allowance,
            "projected_retirement_entry_outside_physical_count": projected_retirement_entry_outside,
            "projected_retirement_exit_outside_physical_count": projected_retirement_exit_outside},
        "detector_window": {"first_lab_time_s": detector_first_lab, "last_lab_time_s": detector_stop_lab,
            "duration_fs": (detector_stop_retarded - detector_first_retarded) * 1e15,
            "first_retarded_time_relative_input_plane_s": detector_first_retarded,
            "last_retarded_time_relative_input_plane_s": detector_stop_retarded,
            "latest_ballistic_retirement_photon_retarded_s": latest_exit_photon_retarded,
            "pre_head_photon_margin_fs": (minima["arrival"] - detector_first_retarded) * 1e15,
            "post_retirement_photon_tail_margin_fs_after_dynamics_allowance": tail_margin_actual_fs,
            "tail_margin_requested_fs": args.tail_margin_fs, "dynamics_delay_allowance_fs": args.dynamics_delay_fs,
            "z_untranslated_at_start_stop_m": [detector_z_at_start, detector_z_at_stop],
            "physical_margin_cells_start_stop": [plane_start_margin_cells, plane_end_margin_cells]},
        "mesh_geometry": {"extent_z_m": extent_z, "physical_bounds_z_m": [physical_lower, physical_upper],
            "translation_z_m": translation, "content_z_untranslated_m": [content_lower, content_upper],
            "radiation_upper_untranslated_m": radiation_upper, "guard_alignment_cells": args.guard_alignment_cells,
            "window_guard_cells": args.window_guard_cells, "z_cpml_cells": args.z_cpml_cells},
        "retirement_protection": {"minimum_photon_energy_eV": args.minimum_photon_energy_ev,
            "requested_cycles": args.protection_cycles, "observed_cycles": protected_cycles,
            "minimum_taper_length_m": required_taper_length, "causal_distance_m": causal_distance,
            "minimum_detector_z_m": minimum_detector_z},
        "assumptions": [
            "Fixed input and activation lab plane are identical; held pre-activation x/y exactly equal input x/y.",
            "The automatic layout is reproduced for cowan-z, zero self-field RMS padding, and the declared guards.",
            "The tail design uses the latest free-drift retirement-exit photon plus an explicit dynamics-delay allowance and the requested quiet tail.",
            "Kinematic bounds do not guarantee evolved magnetic/self-field trajectories; check final retired/completed counts and saved field tail.",
            "Projected transverse retirement counts use initial straight-line velocities and are diagnostic estimates, not evolved-particle losses.",
            "No zero-radiation comparison is scheduled or required by this geometry planner.",
        ],
    }


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-file")
    p.add_argument("--historical-plan")
    p.add_argument("--output")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--chunk-rows", type=int, default=262144)
    p.add_argument("--boost-gamma", type=float, default=929.13054928585427)
    p.add_argument("--dx", type=float, default=5e-6)
    p.add_argument("--dy", type=float, default=5e-6)
    p.add_argument("--dz", type=float, default=1.2543262415359034e-6)
    p.add_argument("--cells-x", type=int, default=160)
    p.add_argument("--cells-y", type=int, default=160)
    p.add_argument("--transverse-cpml-cells", type=int, default=10)
    p.add_argument("--z-cpml-cells", type=int, default=20)
    p.add_argument("--detector-z", type=float, default=7.5)
    p.add_argument("--retirement-entrance-z", type=float, default=4.03)
    p.add_argument("--retirement-length", type=float, default=2.0)
    p.add_argument("--stop-z", type=float, default=50.0, help="Candidate lower bound; extended if tail requires it")
    p.add_argument("--stop-rounding-m", type=float, default=0.1)
    p.add_argument("--tail-margin-fs", type=float, default=100.0)
    p.add_argument("--dynamics-delay-fs", type=float, default=20.0)
    p.add_argument("--window-guard-cells", type=int, default=64)
    p.add_argument("--minimum-lower-guard-cells", type=int, default=64)
    p.add_argument("--upper-guard-cells", type=int, default=128)
    p.add_argument("--guard-alignment-cells", type=int, default=256)
    p.add_argument("--mesh-alignment-cells", type=int, default=256)
    p.add_argument("--ballistic-time-safety-factor", type=float, default=1.1)
    p.add_argument("--duration", type=float, default=5e-10)
    p.add_argument("--minimum-photon-energy-ev", type=float, default=50.0)
    p.add_argument("--protection-cycles", type=int, default=2)
    return p


def self_test(args):
    with tempfile.TemporaryDirectory(prefix="aprl-retirement-geometry-fixture-") as directory:
        root = Path(directory)
        dtype = np.dtype([("plane_position_m", "f8", (2,)), ("arrival_time_offset_s", "f8"),
                          ("proper_velocity", "f8", (3,)), ("macro_weight", "f8")])
        rows = np.zeros(2, dtype=dtype)
        rows["proper_velocity"][:, 2] = math.sqrt(1174.0**2 - 1.0)
        rows["arrival_time_offset_s"] = [-20e-15, 20e-15]
        rows["macro_weight"] = 1.0
        input_file = root / "input.h5"
        with h5py.File(input_file, "w") as handle:
            group = handle.create_group("particles")
            group.attrs["format_version"] = 4
            group.create_dataset("records", data=rows)
        beta = math.sqrt(1.0 - 1.0 / 1174.0**2)
        common = 100e-15
        plan = {"records": 2, "input_plane_z_m": -0.05,
                "initial_center_z_m": -0.05 + beta * C * common,
                "common_time_after_reference_s": common}
        plan_file = root / "plan.json"
        plan_file.write_text(json.dumps(plan))
        args.input_file, args.historical_plan = str(input_file), str(plan_file)
        args.stop_z = 7.6  # Force the minimum-time calculation to extend it.
        result = analyze(args)
        assert result["pass"], result["checks"]
        assert result["stop_z_m"] > args.stop_z
        # Independent direct lab-to-boost transform of the final reference event.
        g = args.boost_gamma
        b = math.sqrt(1.0 - 1.0 / g**2)
        t0 = result["held_geometry"]["frame_origin_lab_time_s"]
        reference_box_z = g * (plan["initial_center_z_m"] + 0.05) + g * b * C * t0
        endpoint_reference_z = -0.05 + g * reference_box_z + g * b * C * result["expected_stop_box_s"]
        assert math.isclose(endpoint_reference_z, result["stop_z_m"], abs_tol=1e-10)
        assert result["detector_window"]["post_retirement_photon_tail_margin_fs_after_dynamics_allowance"] >= 100.0
        assert result["recommended_lower_guard_cells"] % args.guard_alignment_cells == 0
        print("geometry self-test PASS", flush=True)


def main():
    p = parser()
    args = p.parse_args()
    for name in ("chunk_rows", "boost_gamma", "dx", "dy", "dz", "cells_x", "cells_y",
                 "stop_rounding_m", "guard_alignment_cells", "mesh_alignment_cells", "duration"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            p.error(f"{name} must be finite and positive")
    if args.boost_gamma <= 1 or args.ballistic_time_safety_factor < 1:
        p.error("A moving relativistic boost and ballistic safety factor >= 1 are required")
    if args.self_test:
        self_test(args)
        return
    if not all((args.input_file, args.historical_plan, args.output)):
        p.error("--input-file, --historical-plan and --output are required")
    if not os.environ.get("SLURM_JOB_ID"):
        p.error("Full-input geometry scans must run inside a scheduled Slurm job; use --self-test for a tiny local fixture")
    output = Path(args.output)
    if output.exists():
        p.error(f"Refusing to overwrite {output}")
    result = analyze(args)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp-{os.getpid()}")
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    print(json.dumps({key: result[key] for key in ("pass", "stop_z_m", "recommended_lower_guard_cells",
        "recommended_upper_guard_cells", "estimated_cells_z", "expected_stop_box_s")}), flush=True)
    if not result["pass"]:
        raise SystemExit("Geometry preflight failed; inspect checks in the written JSON")


if __name__ == "__main__":
    main()
