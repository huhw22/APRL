#!/usr/bin/env python3
"""Prepare/check one bounded held-injection + retirement functionality case.

This is not a radiation normalization test or a paired baseline calculation.
The full-input geometry scan is a separate single-CPU step in preflight.sbatch.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import yaml

C = 299792458.0
N = 512
BINARY_SHA = "56cc1e0931103687dcc440c3ad390a923b2e74f5da0d68ef1806e3992ccbfeac"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(root, binary):
    if digest(binary) != BINARY_SHA:
        raise ValueError("Executable differs from accepted compute binary")
    case = root / "functional"
    case.mkdir(parents=True, exist_ok=False)
    dtype = np.dtype([("plane_position_m", "<f8", (2,)),
                      ("arrival_time_offset_s", "<f8"),
                      ("proper_velocity", "<f8", (3,)),
                      ("source_id", "<u8"), ("macro_weight", "<f8")])
    records = np.zeros(N, dtype=dtype)
    i = np.arange(N)
    records["plane_position_m"][:, 0] = ((i % 8) - 3.5) * .5e-6
    records["plane_position_m"][:, 1] = (((i // 8) % 8) - 3.5) * .5e-6
    records["arrival_time_offset_s"] = ((i // 64) - 3.5) * .5e-15
    records["proper_velocity"][:, 0] = ((i % 8) - 3.5) * 2e-4
    records["proper_velocity"][:, 1] = (((i // 8) % 8) - 3.5) * 2e-4
    records["proper_velocity"][:, 2] = np.sqrt(15.0)
    records["source_id"] = i + 1001
    records["macro_weight"] = 1.0
    particle = case / "particles.h5"
    with h5py.File(particle, "x") as f:
        g = f.create_group("particles")
        g.attrs["format_version"] = 4
        g.attrs["arrival_time_reference"] = "synthetic-zero-mean-event-time"
        g.create_dataset("records", data=records)
    card = {
        "units": {"length": "micrometer", "time": "femtosecond"},
        "runtime": {"mode": "throughput", "resource_monitor": {"enabled": True,
            "progress_interval_steps": 50, "calibration_steps": 2},
            "shutdown_watchdog": {"enabled": True, "timeout_seconds": 120,
                                  "finalize_timeout_seconds": 30}},
        "output": {"overwrite": False, "manifest": str(case / "run-manifest.yaml")},
        "mesh": {"field_solver": "cowan-z", "current_filter": "none",
            "ballistic_outside_transverse": False, "self_field_feedback": True,
            "cells": [16, 16, 4096], "cell_size": [2, 2, .5], "center": [0, 0, 0],
            "duration": 2500, "boost_gamma": 2,
            "particle_steps_per_undulator_period": 4},
        "boundary": {"type": "cpml", "cells": [3, 3, 8], "polynomial_order": 3,
            "target_reflection": 1e-8, "kappa_max": 8, "alpha_fraction": 0},
        "initial_self_field": {"enabled": True, "model": "relativistic-poisson",
            "relative_tolerance": 1e-6, "maximum_iterations": 2000,
            "maximum_mean_transverse_beta": 1e-7},
        "beam": {"reference": {"input_plane_z": 0, "initial_center_z": 1,
            "frozen_lab_plane_injection": True, "activation_plane_z": 0,
            "hold_transverse_until_activation": True},
            "input": {"type": "hdf5", "file": str(particle), "electrons": N,
                      "position_offset": [0, 0, 0]}},
        "sources": {"incident_waves": [], "magnetic_elements": [{
            "type": "uniform-dipole", "entrance_z": 5, "length": 1,
            "field_T": .01, "polarization_angle_rad": 0}]},
        "particle_retirement": {"enabled": True, "entrance_z": 20, "length": 10},
        "detectors": {"enabled": True, "directory": str(case / "detectors"),
            "field_planes": [{"name": "post-retirement-field", "z": 250,
                "rhythm": 1, "buffer_samples": 8, "compression": 0,
                "particle_background": {"enabled": False},
                "retirement_frequency_protection": {"enabled": True,
                    "minimum_photon_energy_eV": 50, "cycles": 2}}],
            "particle_planes": [{"name": name, "z": z, "buffer_records": 512,
                "compression": 0} for name, z in (("activation", 0), ("retirement-entry", 20))]},
        "stop": {"mode": "reference-center-z", "z": 1000},
        "trajectory": {"enabled": False},
        "energy_ledger": {"enabled": True, "directory": str(case / "ledger"),
            "filename": "energy-ledger.h5", "sample_interval_steps": 1,
            "buffer_records": 32, "compression": 0, "warning_relative_tolerance": .05}}
    (case / "run.yaml").write_text(yaml.safe_dump(card, sort_keys=False))
    plan = {"binary": str(binary), "binary_sha256": BINARY_SHA,
            "input_sha256": digest(particle), "particles": N,
            "purpose": "Single functionality test; no radiation validation or comparison",
            "reserved_nodes": 2, "ranks": 256, "threads_per_rank": 1}
    (root / "preflight-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print(case / "run.yaml")


def read_complete(path, name):
    with h5py.File(path, "r") as f:
        g = f[name]
        field = name == "field_plane"
        n = int(g["committed_samples" if field else "committed_records"][()])
        data = g["fields" if field else "records"][:]
        if int(g["complete"][()]) != 1 or n <= 0 or len(data) != n:
            raise ValueError(f"Incomplete output: {path}")
        times = g["time_s"][:] if field else None
        return data, times, dict(g.attrs)


def check(root, binary):
    case = root / "functional"
    plan = json.loads((root / "preflight-plan.json").read_text())
    if digest(binary) != BINARY_SHA or digest(case / "particles.h5") != plan["input_sha256"]:
        raise ValueError("Frozen binary or input changed")
    geometry = json.loads((root / "geometry.json").read_text())
    if geometry.get("pass") is not True:
        raise ValueError("Full-input geometry did not pass")
    log = (case / "aprl.log").read_text(errors="replace")
    if "Configured stop reached:" not in log or "reference" not in log.split("Configured stop reached:")[-1]:
        raise ValueError("Missing configured reference stop")
    retired = re.search(r"Experimental particle retirement: entries=(\d+).*?completed=(\d+), residual Q at taper exit=([^ ]+)", log)
    activated = re.search(r"Frozen lab-plane injection: activations=(\d+).*?still frozen=(\d+)", log)
    if not retired or tuple(map(int, retired.group(1, 2))) != (N, N) or float(retired.group(3)) != 0:
        raise ValueError("Not all retirement carriers completed with zero residual charge")
    if not activated or tuple(map(int, activated.groups())) != (N, 0):
        raise ValueError("Held injection did not activate the full cohort")
    manifest = yaml.safe_load((case / "run-manifest.yaml").read_text())
    if manifest["mpi_size"] != 256 or manifest["source_revision"] != "ac58e67e33d6b4ecfd12c3d00f4a9947d0e8500d":
        raise ValueError("Unexpected MPI layout or compute revision")
    with h5py.File(case / "particles.h5", "r") as f:
        source = f["particles/records"][:]
    plane_records = {}
    for name in ("activation", "retirement-entry"):
        a, _, _ = read_complete(case / "detectors" / (name + ".h5"), "particle_plane")
        a = a[np.argsort(a["source_id"])]
        if len(a) != N or not np.array_equal(a["source_id"], source["source_id"]):
            raise ValueError("Particle detector cohort differs from input")
        plane_records[name] = a
    activation_error = float(np.abs(plane_records["activation"]["position_m"][:, :2] - source["plane_position_m"]).max())
    if activation_error > 1e-12:
        raise ValueError("Held activation did not preserve input-plane transverse position")
    a = plane_records["retirement-entry"]
    beta_z = a["proper_velocity"][:, 2] / np.sqrt(1 + (a["proper_velocity"]**2).sum(axis=1))
    last_carrier_lab_time = float(np.max(a["time_s"] + 10e-6 / (C * beta_z)))
    fields, time, attrs = read_complete(case / "detectors/post-retirement-field.h5", "field_plane")
    if not np.all(np.diff(time) > 0) or len(time) != len(fields):
        raise ValueError("Invalid field time axis")
    for name in ("electric_V_per_m", "magnetic_T"):
        if not np.isfinite(fields[name]).all():
            raise ValueError("Nonfinite field sample")
    after = int(np.count_nonzero(time > last_carrier_lab_time + 20e-15))
    if after < 20 or b"enabled" not in str(attrs["retirement_frequency_protection"]).encode():
        raise ValueError("Insufficient guarded field recording after carrier exhaustion")
    ledger, _, _ = read_complete(case / "ledger/energy-ledger.h5", "energy_ledger")
    empty = np.flatnonzero(ledger["active_macroparticles"] == 0)
    if not len(empty) or int(ledger[-1]["removed_macroparticles"]) != N or len(ledger) - empty[0] < 20:
        raise ValueError("Field stepping did not continue after all real particles retired")
    report = {"pass": True, "functional_only": True, "particles": N,
        "retirement_entries": N, "completed_carriers": N, "residual_exit_charge_C": 0,
        "activation_xy_max_error_m": activation_error,
        "field_samples_after_last_carrier_plus_20fs": after,
        "last_carrier_lab_time_s_estimated_from_entry": last_carrier_lab_time,
        "field_last_time_lab_s": float(time[-1]),
        "ledger_samples_after_zero_active": int(len(ledger) - empty[0]),
        "geometry_pass": True, "binary_sha256": BINARY_SHA}
    (root / "preflight-result.json").write_text(json.dumps(report, indent=2) + "\n")
    (root / "PASS").write_text("Binding, production-input geometry, and single held-retirement functionality passed.\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("prepare", "check"))
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    args = p.parse_args()
    (prepare if args.mode == "prepare" else check)(args.root.resolve(), args.binary.resolve())
