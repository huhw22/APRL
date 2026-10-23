#!/usr/bin/env python3
"""Analyze the APRL 4 m SASE exact-unit-weight case."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable

try:
    import h5py
    import numpy as np
except ImportError as exc:  # pragma: no cover - cluster environment dependent
    raise SystemExit(
        "ERROR: analysis needs numpy and h5py; run prepare_python_env.sh first"
    ) from exc


C = 299_792_458.0
FACE_NAMES = ("x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus")
RESOURCE_RE = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
CELLS_RE = re.compile(
    r"resolved_mesh:\s*\n\s+cells:\s*\[([0-9]+),\s*([0-9]+),\s*([0-9]+)\]"
)
EXTENT_RE = re.compile(
    r"\s+extent_m:\s*\[([^,]+),\s*([^,]+),\s*([^\]]+)\]"
)
TRANSLATION_RE = re.compile(r"\s+translation_m:\s*([^\s]+)")
SELF_FIELD_PADDING_RE = re.compile(r"\s+self_field_padding_rms:\s*([^\s]+)")


def fail(message: str) -> None:
    raise RuntimeError(message)


def scalar(value: Any) -> int:
    array = np.asarray(value)
    if array.size != 1:
        fail("expected a scalar HDF5 completion/count value")
    return int(array.reshape(-1)[0])


def require_file(path: Path, description: str) -> Path:
    if not path.is_file():
        fail(f"missing {description}: {path}")
    return path


def group_complete(group: h5py.Group) -> bool:
    if "complete" in group:
        return scalar(group["complete"][()]) == 1
    if "complete" in group.attrs:
        return scalar(group.attrs["complete"]) == 1
    return False


def require_fields(dtype: np.dtype, names: Iterable[str], path: Path) -> None:
    available = set(dtype.names or ())
    missing = [name for name in names if name not in available]
    if missing:
        fail(f"{path} is missing record fields: {', '.join(missing)}")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        fail(f"refusing to write empty CSV: {path.name}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value: Any) -> float:
    result = float(value)
    return result if math.isfinite(result) else float("nan")


def load_plan(path: Path) -> dict[str, Any]:
    try:
        plan = json.loads(require_file(path, "campaign plan").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read campaign plan {path}: {exc}")
    if (not plan.get("cases") or "theory" not in plan or
            not plan.get("campaign_fingerprint")):
        fail(f"campaign plan is incomplete: {path}")
    return plan


def validate_completion(run: Path, expected_campaign: str, case_tag: str) -> None:
    marker = require_file(run / "COMPLETE", "completion marker")
    fields: dict[str, str] = {}
    for line in marker.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields[key] = value
    if fields.get("status") != "COMPLETE" or not fields.get("fingerprint"):
        fail(f"{case_tag}: malformed or legacy COMPLETE marker: {marker}")
    if fields.get("campaign_fingerprint") != expected_campaign:
        fail(
            f"{case_tag}: completed result belongs to campaign "
            f"{fields.get('campaign_fingerprint', '<missing>')}, but current plan is "
            f"{expected_campaign}; rerun the case or restore its submitted plan"
        )
    stored_run_fingerprint = require_file(
        run / "RUN_FINGERPRINT", "run fingerprint"
    ).read_text(encoding="utf-8").strip()
    if stored_run_fingerprint != fields["fingerprint"]:
        fail(f"{case_tag}: RUN_FINGERPRINT disagrees with COMPLETE")
    submitted = load_plan(run / "submitted-plan.json")
    if submitted["campaign_fingerprint"] != expected_campaign:
        fail(f"{case_tag}: submitted plan disagrees with current campaign plan")
    try:
        input_snapshot = json.loads(require_file(
            run / "input-noise.json", "SASE input metadata snapshot"
        ).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"{case_tag}: cannot read SASE input metadata snapshot: {exc}")
    if (input_snapshot.get("input_sha256") != fields.get("input_sha256") or
            input_snapshot.get("generator_sha256") !=
            fields.get("generator_sha256")):
        fail(f"{case_tag}: input-noise snapshot disagrees with COMPLETE provenance")


def read_manifest(path: Path) -> dict[str, Any]:
    text = require_file(path, "run manifest").read_text(
        encoding="utf-8", errors="replace"
    )
    cells = CELLS_RE.search(text)
    extent = EXTENT_RE.search(text)
    translation = TRANSLATION_RE.search(text)
    if not cells or not extent:
        fail(f"manifest has no resolved automatic mesh: {path}")
    if "initial_self_field_enabled: true" not in text:
        fail(f"manifest does not record the required initialized self field: {path}")
    padding = SELF_FIELD_PADDING_RE.search(text)
    if not padding or float(padding.group(1)) != 0.0:
        fail(f"manifest does not record self_field_padding_rms=0: {path}")
    return {
        "resolved_nx": int(cells.group(1)),
        "resolved_ny": int(cells.group(2)),
        "resolved_nz": int(cells.group(3)),
        "resolved_extent_z_m": float(extent.group(3)),
        "box_translation_m": float(translation.group(1)) if translation else float("nan"),
        "initial_self_field_enabled": True,
        "self_field_padding_rms": float(padding.group(1)),
    }


def read_ledger(path: Path, case_tag: str, boost_gamma: float) -> tuple[
    list[dict[str, Any]], dict[str, Any]
]:
    required = (
        "step", "time_box_s", "closure_valid", "particle_kinetic_active_J",
        "particle_kinetic_removed_J", "field_energy_interior_J",
        "outward_field_energy_J", "outward_field_energy_face_J",
        "prescribed_work_J", "balance_residual_J", "relative_balance_initial",
        "active_macroparticles", "removed_macroparticles", "mean_gamma_lab",
        "sigma_gamma_lab", "uncorrelated_sigma_gamma_lab",
    )
    try:
        with h5py.File(require_file(path, "energy ledger"), "r") as handle:
            group = handle["/energy_ledger"]
            if not group_complete(group):
                fail(f"energy ledger is incomplete: {path}")
            committed = scalar(group["committed_records"][()])
            dataset = group["records"]
            if committed < 2 or committed > dataset.shape[0]:
                fail(f"invalid committed record count {committed}: {path}")
            require_fields(dataset.dtype, required, path)
            records = dataset[:committed]
    except (OSError, KeyError) as exc:
        fail(f"cannot read energy ledger {path}: {exc}")

    beta = math.sqrt(1.0 - 1.0 / boost_gamma**2)
    advance = boost_gamma * beta * C * np.asarray(records["time_box_s"], dtype=np.float64)
    residual = np.abs(np.asarray(records["relative_balance_initial"], dtype=np.float64))
    log_residual = np.log10(np.maximum(residual, np.finfo(np.float64).tiny))
    field = np.asarray(records["field_energy_interior_J"], dtype=np.float64)
    if not (field[0] > 0.0) or not math.isfinite(float(field[0])):
        fail(
            f"{case_tag}: enabled initial self-field did not produce a "
            "positive finite initial field energy"
        )
    field_ratio = field / field[0]
    mean_gamma = np.asarray(records["mean_gamma_lab"], dtype=np.float64)
    sigma_gamma = np.asarray(records["sigma_gamma_lab"], dtype=np.float64)
    unc_sigma = np.asarray(records["uncorrelated_sigma_gamma_lab"], dtype=np.float64)
    active = np.asarray(records["active_macroparticles"], dtype=np.uint64)
    removed = np.asarray(records["removed_macroparticles"], dtype=np.uint64)
    initial_total = int(active[0]) + int(removed[0])
    if initial_total <= 0:
        fail(f"ledger has no initial particles: {path}")
    faces = np.asarray(records["outward_field_energy_face_J"], dtype=np.float64)
    if faces.ndim != 2 or faces.shape[1] != 6:
        fail(f"ledger face energy is not six-dimensional: {path}")
    delta_faces = faces - faces[0]
    outward = np.asarray(records["outward_field_energy_J"], dtype=np.float64)
    particle_active = np.asarray(records["particle_kinetic_active_J"], dtype=np.float64)
    particle_removed = np.asarray(records["particle_kinetic_removed_J"], dtype=np.float64)
    prescribed = np.asarray(records["prescribed_work_J"], dtype=np.float64)

    rows: list[dict[str, Any]] = []
    for index in range(records.size):
        row: dict[str, Any] = {
            "case": case_tag,
            "step": int(records["step"][index]),
            "time_box_s": float(records["time_box_s"][index]),
            "boost_reference_advance_m": float(advance[index]),
            "closure_valid": int(records["closure_valid"][index]),
            "abs_relative_residual": float(residual[index]),
            "log10_abs_relative_residual": float(log_residual[index]),
            "balance_residual_J": float(records["balance_residual_J"][index]),
            "field_energy_interior_J": float(field[index]),
            "field_energy_over_initial": float(field_ratio[index]),
            "delta_field_energy_J": float(field[index] - field[0]),
            "particle_kinetic_active_J": float(particle_active[index]),
            "particle_kinetic_removed_J": float(particle_removed[index]),
            "delta_particle_kinetic_active_J": float(particle_active[index] - particle_active[0]),
            "delta_particle_kinetic_removed_J": float(particle_removed[index] - particle_removed[0]),
            "outward_field_energy_J": float(outward[index] - outward[0]),
            "prescribed_work_J": float(prescribed[index] - prescribed[0]),
            "mean_gamma_lab": float(mean_gamma[index]),
            "delta_mean_gamma_lab": float(mean_gamma[index] - mean_gamma[0]),
            "sigma_gamma_lab": float(sigma_gamma[index]),
            "uncorrelated_sigma_gamma_lab": float(unc_sigma[index]),
            "active_macroparticles": int(active[index]),
            "removed_macroparticles": int(removed[index]),
            "removed_fraction": float(removed[index]) / float(initial_total),
        }
        for face_index, face_name in enumerate(FACE_NAMES):
            row[f"outward_{face_name}_J"] = float(delta_faces[index, face_index])
        rows.append(row)

    evolved_positive = field[(np.arange(field.size) > 0) & (field > 0.0)]
    first_positive = float(evolved_positive[0]) if evolved_positive.size else float("nan")
    final = rows[-1]
    summary = {
        "records": len(rows),
        "final_step": final["step"],
        "final_boost_reference_advance_m": final["boost_reference_advance_m"],
        "closure_valid_all": bool(np.all(records["closure_valid"] == 1)),
        "peak_abs_relative_residual": float(np.nanmax(residual)),
        "final_abs_relative_residual": final["abs_relative_residual"],
        "initial_field_energy_J": float(field[0]),
        "first_positive_field_energy_J": first_positive,
        "peak_field_energy_J": float(np.nanmax(field)),
        "final_field_energy_J": float(field[-1]),
        "peak_field_energy_over_initial": float(np.nanmax(field_ratio)),
        "final_field_energy_over_initial": float(field_ratio[-1]),
        "final_over_first_positive_field": (
            float(field[-1] / first_positive) if first_positive > 0.0 else float("nan")
        ),
        "initial_mean_gamma_lab": float(mean_gamma[0]),
        "final_mean_gamma_lab": float(mean_gamma[-1]),
        "delta_mean_gamma_lab": float(mean_gamma[-1] - mean_gamma[0]),
        "initial_sigma_gamma_lab": float(sigma_gamma[0]),
        "final_sigma_gamma_lab": float(sigma_gamma[-1]),
        "final_uncorrelated_sigma_gamma_lab": float(unc_sigma[-1]),
        "removed_fraction": final["removed_fraction"],
        "delta_particle_kinetic_active_J": final["delta_particle_kinetic_active_J"],
        "delta_particle_kinetic_removed_J": final["delta_particle_kinetic_removed_J"],
        "delta_particle_kinetic_total_J": (
            final["delta_particle_kinetic_active_J"] +
            final["delta_particle_kinetic_removed_J"]
        ),
        "delta_field_energy_J": final["delta_field_energy_J"],
        "outward_field_energy_J": final["outward_field_energy_J"],
        "prescribed_work_J": final["prescribed_work_J"],
        "final_balance_residual_J": final["balance_residual_J"],
        "outward_field_energy_faces_J": {
            name: float(delta_faces[-1, index])
            for index, name in enumerate(FACE_NAMES)
        },
    }
    return rows, summary


def read_particle_plane(
    path: Path, case_tag: str, label: str, wavelength_m: float,
    expected_records: int, harmonics: int = 3,
) -> dict[str, Any]:
    required = (
        "particle_id", "time_s", "proper_velocity", "charge_C",
        "crossing_direction",
    )
    omega = 2.0 * math.pi * C / wavelength_m
    count = 0
    raw_downstream_records = 0
    duplicate_records = 0
    invalid_id_records = 0
    seen = np.zeros(expected_records + 1, dtype=np.bool_)
    sum_weight = np.longdouble(0.0)
    sum_gamma = np.longdouble(0.0)
    sum_gamma2 = np.longdouble(0.0)
    sum_beta_t2 = np.longdouble(0.0)
    real = np.zeros(harmonics, dtype=np.longdouble)
    imag = np.zeros(harmonics, dtype=np.longdouble)
    cell_weight: dict[int, np.longdouble] = {}
    cell_real: list[dict[int, np.longdouble]] = [
        {} for _ in range(harmonics)
    ]
    cell_imag: list[dict[int, np.longdouble]] = [
        {} for _ in range(harmonics)
    ]
    phase_reference: float | None = None
    try:
        with h5py.File(require_file(path, f"{label} particle plane"), "r") as handle:
            group = handle["/particle_plane"]
            if not group_complete(group):
                fail(f"{label} particle plane is incomplete: {path}")
            committed = scalar(group["committed_records"][()])
            dataset = group["records"]
            if committed <= 0 or committed > dataset.shape[0]:
                fail(f"invalid committed record count for {label}: {path}")
            require_fields(dataset.dtype, required, path)
            selected = dataset.fields(list(required))
            for start in range(0, committed, 262_144):
                batch = selected[start:min(committed, start + 262_144)]
                batch = batch[batch["crossing_direction"] > 0]
                if batch.size == 0:
                    continue
                raw_downstream_records += int(batch.size)
                ids = np.asarray(batch["particle_id"], dtype=np.uint64)
                valid_ids = (ids >= 1) & (ids <= expected_records)
                invalid_id_records += int(np.count_nonzero(~valid_ids))
                batch = batch[valid_ids]
                ids = ids[valid_ids]
                if batch.size == 0:
                    continue
                candidate_indices = np.flatnonzero(~seen[ids])
                candidate_ids = ids[candidate_indices]
                if candidate_ids.size:
                    _, first = np.unique(candidate_ids, return_index=True)
                    keep = candidate_indices[np.sort(first)]
                else:
                    keep = np.empty(0, dtype=np.int64)
                duplicate_records += int(batch.size - keep.size)
                batch = batch[keep]
                ids = ids[keep]
                if batch.size == 0:
                    continue
                seen[ids] = True
                times = np.asarray(batch["time_s"], dtype=np.float64)
                if phase_reference is None:
                    phase_reference = float(times[0])
                phase = omega * (times - phase_reference)
                wavelength_cell = np.floor(phase / (2.0 * math.pi)).astype(np.int64)
                velocity = np.asarray(batch["proper_velocity"], dtype=np.float64)
                gamma = np.sqrt(1.0 + np.einsum("ij,ij->i", velocity, velocity))
                weight = np.abs(np.asarray(batch["charge_C"], dtype=np.float64))
                beta_t2 = (velocity[:, 0] ** 2 + velocity[:, 1] ** 2) / gamma**2
                count += int(batch.size)
                sum_weight += np.sum(weight, dtype=np.longdouble)
                sum_gamma += np.sum(weight * gamma, dtype=np.longdouble)
                sum_gamma2 += np.sum(weight * gamma**2, dtype=np.longdouble)
                sum_beta_t2 += np.sum(weight * beta_t2, dtype=np.longdouble)
                unique_cells, inverse = np.unique(wavelength_cell, return_inverse=True)
                weights_by_cell = np.bincount(inverse, weights=weight)
                for local_index, cell in enumerate(unique_cells):
                    key = int(cell)
                    cell_weight[key] = cell_weight.get(key, np.longdouble(0.0)) + \
                        np.longdouble(weights_by_cell[local_index])
                for harmonic in range(1, harmonics + 1):
                    cosine = weight * np.cos(harmonic * phase)
                    sine = weight * np.sin(harmonic * phase)
                    real[harmonic - 1] += np.sum(cosine, dtype=np.longdouble)
                    imag[harmonic - 1] += np.sum(sine, dtype=np.longdouble)
                    real_by_cell = np.bincount(inverse, weights=cosine)
                    imag_by_cell = np.bincount(inverse, weights=sine)
                    for local_index, cell in enumerate(unique_cells):
                        key = int(cell)
                        real_map = cell_real[harmonic - 1]
                        imag_map = cell_imag[harmonic - 1]
                        real_map[key] = real_map.get(key, np.longdouble(0.0)) + \
                            np.longdouble(real_by_cell[local_index])
                        imag_map[key] = imag_map.get(key, np.longdouble(0.0)) + \
                            np.longdouble(imag_by_cell[local_index])
    except (OSError, KeyError) as exc:
        fail(f"cannot read {label} particle plane {path}: {exc}")

    if count == 0 or sum_weight <= 0:
        fail(f"{label} particle plane has no downstream crossings: {path}")
    mean = float(sum_gamma / sum_weight)
    variance = max(0.0, float(sum_gamma2 / sum_weight) - mean**2)
    result: dict[str, Any] = {
        "case": case_tag,
        "plane": label,
        "records": count,
        "raw_downstream_records": raw_downstream_records,
        "duplicate_downstream_records": duplicate_records,
        "invalid_particle_id_records": invalid_id_records,
        "missing_particle_ids": expected_records - int(np.count_nonzero(seen)),
        "expected_macroparticles": expected_records,
        "record_fraction": count / expected_records,
        "absolute_charge_C": float(sum_weight),
        "mean_gamma": mean,
        "sigma_gamma": math.sqrt(variance),
        "rms_transverse_beta": math.sqrt(max(0.0, float(sum_beta_t2 / sum_weight))),
    }
    for harmonic in range(1, harmonics + 1):
        result[f"bunching_h{harmonic}"] = float(
            math.hypot(float(real[harmonic - 1]), float(imag[harmonic - 1])) /
            float(sum_weight)
        )
        local_power_numerator = np.longdouble(0.0)
        scaled_cell_power: list[float] = []
        for key, charge in cell_weight.items():
            if charge <= 0.0:
                continue
            cell_amplitude_squared = (
                cell_real[harmonic - 1].get(key, np.longdouble(0.0)) ** 2 +
                cell_imag[harmonic - 1].get(key, np.longdouble(0.0)) ** 2
            ) / (charge * charge)
            local_power_numerator += charge * cell_amplitude_squared
            electrons_in_cell = charge / np.longdouble(1.602_176_634e-19)
            scaled_cell_power.append(float(electrons_in_cell * cell_amplitude_squared))
        result[f"slice_rms_bunching_h{harmonic}"] = math.sqrt(max(
            0.0, float(local_power_numerator / sum_weight)
        ))
        result[f"mean_electrons_per_slice_times_bunching_squared_h{harmonic}"] = (
            float(np.mean(scaled_cell_power)) if scaled_cell_power else float("nan")
        )
    result["wavelength_slices"] = len(cell_weight)
    return result


def noise_result_from_attributes(
    attrs: Any, harmonics: int, description: str,
) -> dict[str, Any]:
    observed = np.asarray(attrs["observed_global_bunching_abs"], dtype=np.float64)
    expected = np.asarray(attrs["expected_global_bunching_rms"], dtype=np.float64)
    scaled = np.asarray(
        attrs["mean_electrons_per_cell_times_local_bunching_squared"],
        dtype=np.float64,
    )
    if min(observed.size, expected.size, scaled.size) < harmonics:
        fail(f"SASE input has fewer than {harmonics} noise harmonics: {description}")
    result: dict[str, Any] = {
        "input_metadata_available": True,
        "input_random_seed": int(attrs["random_seed"]),
        "input_wavelength_cells": int(attrs["wavelength_cells"]),
        "input_beamlets_per_wavelength": int(attrs["beamlets_per_wavelength"]),
        "input_electrons_per_beamlet": float(attrs["electrons_per_beamlet"]),
        "input_mean_current_A": float(attrs["mean_current_A"]),
    }
    for harmonic in range(1, harmonics + 1):
        result[f"input_global_bunching_h{harmonic}"] = float(
            observed[harmonic - 1]
        )
        result[f"input_expected_global_rms_h{harmonic}"] = float(
            expected[harmonic - 1]
        )
        result[f"input_scaled_local_power_h{harmonic}"] = float(
            scaled[harmonic - 1]
        )
    return result


def read_input_noise(
    path: Path, snapshot: Path | None = None, harmonics: int = 3,
) -> dict[str, Any]:
    """Read portable or HDF5 authoritative pre-run shot-noise diagnostics."""
    if snapshot is not None and snapshot.is_file():
        try:
            payload = json.loads(snapshot.read_text(encoding="utf-8"))
            if payload.get("schema") != "aprl-sase-input-metadata-v1":
                fail(f"unsupported SASE metadata snapshot: {snapshot}")
            result = noise_result_from_attributes(
                payload["attributes"], harmonics, str(snapshot)
            )
            result["input_sha256"] = str(payload["input_sha256"])
            result["input_generator_sha256"] = str(payload["generator_sha256"])
            return result
        except (OSError, KeyError, json.JSONDecodeError) as exc:
            fail(f"cannot read SASE metadata snapshot {snapshot}: {exc}")
    if not path.is_file():
        return {"input_metadata_available": False}
    try:
        with h5py.File(path, "r") as handle:
            group = handle["/particles"]
            if not group_complete(group):
                fail(f"SASE input is incomplete: {path}")
            return noise_result_from_attributes(group.attrs, harmonics, str(path))
    except (OSError, KeyError) as exc:
        fail(f"cannot read SASE input metadata {path}: {exc}")


def read_resource_log(path: Path, case_tag: str) -> tuple[
    list[dict[str, Any]], dict[str, Any]
]:
    text = require_file(path, "APRL log").read_text(encoding="utf-8", errors="replace")
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        if "[resource]" not in line:
            continue
        row: dict[str, Any] = {"case": case_tag}
        for key, raw in RESOURCE_RE.findall(line):
            try:
                row[key] = float(raw)
            except ValueError:
                row[key] = raw
        rows.append(row)
    if not rows:
        fail(f"no [resource] records found in {path}")

    def maximum(key: str) -> float:
        values = [finite_float(row[key]) for row in rows if key in row]
        values = [value for value in values if math.isfinite(value)]
        return max(values) if values else float("nan")

    def last_finite(key: str) -> float:
        for row in reversed(rows):
            if key in row and math.isfinite(finite_float(row[key])):
                return finite_float(row[key])
        return float("nan")

    summary = {
        "records": len(rows),
        "peak_rss_total_mib": maximum("peak_rss_total_mib"),
        "peak_rss_rank_max_mib": maximum("peak_rss_rank_max_mib"),
        "last_elapsed_s": last_finite("elapsed_s"),
        "last_measured_step_s": last_finite("measured_step_s"),
        "last_reported_step": last_finite("step"),
    }
    return rows, summary


def style_for(case: dict[str, Any]) -> tuple[str, str, str]:
    beamlets = int(case["beamlets_per_wavelength"])
    dz_nm = float(case["dz_lab_m"]) * 1.0e9
    label = f"{dz_nm:.2f} nm, B={beamlets}, weight=1"
    return "#2878b5", "-", label


def make_comparison_plot(
    case_results: list[dict[str, Any]], output: Path, warning: float,
) -> None:
    try:
        os.environ.setdefault("MPLCONFIGDIR", str(output.parent / ".matplotlib"))
        Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        fail(f"plotting needs matplotlib ({exc})")

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.8), constrained_layout=True)
    for result in case_results:
        case = result["case"]
        rows = result["ledger_rows"]
        color, linestyle, label = style_for(case)

        def column(name: str) -> np.ndarray:
            return np.asarray([row[name] for row in rows], dtype=np.float64)

        z = column("boost_reference_advance_m")
        residual = column("abs_relative_residual")
        plotted_log = np.full_like(residual, np.nan)
        # Step zero is the exact R=0 ledger reference.  Omit it only from the
        # plot; the raw zero and floored log value remain in the CSV/JSON.
        evolved = (np.asarray([row["step"] for row in rows]) > 0)
        evolved &= np.isfinite(residual) & (residual > 0.0)
        plotted_log[evolved] = np.log10(residual[evolved])
        axes[0, 0].plot(z, plotted_log, color=color, linestyle=linestyle,
                        linewidth=1.7, label=label)

        field_ratio = column("field_energy_over_initial")
        axes[0, 1].plot(z, field_ratio, color=color, linestyle=linestyle,
                        linewidth=1.7, label=label)
        axes[1, 0].plot(z, column("delta_mean_gamma_lab"), color=color,
                        linestyle=linestyle, linewidth=1.7, label=label)
        axes[1, 1].plot(z, column("uncorrelated_sigma_gamma_lab"), color=color,
                        linestyle=linestyle, linewidth=1.7, label=label)

    axes[0, 0].axhline(math.log10(warning), color="black", linestyle=":", linewidth=1)
    axes[0, 0].set_ylabel(r"$\log_{10}(|R|/E_0)$")
    axes[0, 0].set_title(
        f"Energy closure (dotted: {warning:g}; initial R=0 omitted)"
    )
    axes[0, 1].axhline(1.0, color="black", linestyle=":", linewidth=1)
    axes[0, 1].set_ylabel(r"$U_{field}/U_{field,0}$")
    axes[0, 1].set_title("Interior field energy relative to initialized field")
    axes[1, 0].axhline(0.0, color="black", linestyle=":", linewidth=1)
    axes[1, 0].set_ylabel(r"$\langle\gamma\rangle-\langle\gamma\rangle_0$")
    axes[1, 0].set_title("Particle mean energy")
    axes[1, 1].set_ylabel(r"uncorrelated $\sigma_\gamma$")
    axes[1, 1].set_title("Uncorrelated energy spread")
    for axis in axes.flat:
        axis.set_xlabel(r"boost-reference advance $\gamma_b\beta_bct'$ (m)")
        axis.grid(True, alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle("APRL 4 m SASE: exact unit macro-particle weight")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def make_bunching_plot(case_results: list[dict[str, Any]], output: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        fail(f"plotting needs matplotlib ({exc})")

    labels = [style_for(result["case"])[2] for result in case_results]
    x = np.arange(len(case_results), dtype=np.float64)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), constrained_layout=True)
    for harmonic, axis in enumerate(axes, start=1):
        entry_global = np.asarray([
            result["planes"]["entry"][f"bunching_h{harmonic}"]
            for result in case_results
        ], dtype=np.float64)
        exit_global = np.asarray([
            result["planes"]["exit"][f"bunching_h{harmonic}"]
            for result in case_results
        ], dtype=np.float64)
        entry_local = np.asarray([
            result["planes"]["entry"][f"slice_rms_bunching_h{harmonic}"]
            for result in case_results
        ], dtype=np.float64)
        exit_local = np.asarray([
            result["planes"]["exit"][f"slice_rms_bunching_h{harmonic}"]
            for result in case_results
        ], dtype=np.float64)
        axis.plot(x, entry_local, color="#6c757d", marker="o", linestyle="--",
                  linewidth=1.2, label="entry: 1λ slice RMS")
        axis.plot(x, exit_local, color="#7b2cbf", marker="s", linestyle="-",
                  linewidth=1.5, label="exit: 1λ slice RMS")
        axis.plot(x, entry_global, color="#adb5bd", marker=".", linestyle=":",
                  linewidth=1.0, alpha=0.8, label="entry: whole bunch")
        axis.plot(x, exit_global, color="#c77dff", marker=".", linestyle=":",
                  linewidth=1.0, alpha=0.8, label="exit: whole bunch")
        axis.set_yscale("log")
        axis.set_title(f"Harmonic {harmonic}")
        axis.set_ylabel(r"$|b_h|$")
        axis.set_xticks(x, labels, rotation=28, ha="right")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(loc="best", fontsize=8)
    fig.suptitle(
        "Fixed-lab-plane microbunching; local RMS is primary for SASE"
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def fmt(value: Any) -> str:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.8g}" if math.isfinite(float(value)) else "n/a"
    return str(value)


def sanitize(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: sanitize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    return value


def write_report(
    path: Path, plan: dict[str, Any], case_results: list[dict[str, Any]], warning: float,
) -> None:
    reference = plan["reference_case"]
    charge_ratio = float(plan["beam"]["electrons"]) / float(reference["electrons"])
    lines = [
        "# APRL 4 m SASE exact-unit-weight report", "",
        "This is one extreme macro-particle-resolution run at dz=1.35 nm, not a SASE "
        "ensemble result. It is compared with the already completed dz1p35-B0192 run "
        "only as a numerical reference; the two noise realizations are not pointwise paired.", "",
        "The input is a uniform 40 um top-hat. Exact unit weight and equal 256-rank "
        f"partitioning require `{int(plan['beam']['electrons']):,}` electrons/macroparticles, "
        f"which is `{100.0 * (charge_ratio - 1.0):.5g}%` above the old 18,000,000-electron "
        f"case; mean current is `{fmt(plan['beam']['mean_current_A'])} A`. The run initializes "
        "a relativistic-Poisson field on the resolved non-CPML interval with "
        "self_field_padding_rms=0. The state "
        "therefore satisfies the discrete Gauss constraint instead of starting from zero, "
        "but the omitted exterior Coulomb tail remains a finite-box approximation. "
        "No K=0 baseline or fixed field plane is used.", "",
        "The first energy-ledger record defines R=0 exactly. It remains in CSV/JSON but is "
        "omitted from the log-residual plot so it cannot compress the evolved curve.", "",
        "## Comparison", "",
        "| case | dz (nm) | macros | e/macro | resolved Nz | peak abs(R)/E0 | initial Ufield (J) | peak Ufield/U0 | final Ufield/U0 | delta mean gamma | final uncorr sigma | exit local-rms b1 | local-rms b1 gain | elapsed (h) | peak RSS (GiB) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in case_results:
        case = result["case"]
        ledger = result["ledger"]
        resources = result["resources"]
        lines.append(
            f"| {case['tag']} | {float(case['dz_lab_m']) * 1e9:.3f} | "
            f"{int(case['macroparticles']):,} | "
            f"{float(case['electrons_per_macroparticle']):.6g} | "
            f"{int(result['manifest']['resolved_nz']):,} | "
            f"{fmt(ledger['peak_abs_relative_residual'])} | "
            f"{fmt(ledger['initial_field_energy_J'])} | "
            f"{fmt(ledger['peak_field_energy_over_initial'])} | "
            f"{fmt(ledger['final_field_energy_over_initial'])} | "
            f"{fmt(ledger['delta_mean_gamma_lab'])} | "
            f"{fmt(ledger['final_uncorrelated_sigma_gamma_lab'])} | "
            f"{fmt(result['planes']['exit']['slice_rms_bunching_h1'])} | "
            f"{fmt(result['planes']['exit']['slice_rms_bunching_h1'] / result['planes']['entry']['slice_rms_bunching_h1'])} | "
            f"{fmt(resources['last_elapsed_s'] / 3600.0)} | "
            f"{fmt(resources['peak_rss_total_mib'] / 1024.0)} |"
        )
    result = case_results[0]
    ledger = result["ledger"]
    resources = result["resources"]
    residual_ratio = (
        float(ledger["peak_abs_relative_residual"]) /
        float(reference["peak_abs_relative_residual"])
    )
    field_ratio_raw = (
        float(ledger["final_field_energy_J"]) /
        float(reference["final_field_energy_J"])
    )
    field_ratio_charge_normalized = field_ratio_raw / (charge_ratio * charge_ratio)
    elapsed_ratio = (
        float(resources["last_elapsed_s"]) / float(reference["elapsed_s"])
    )
    rss_ratio = (
        float(resources["peak_rss_total_mib"]) /
        float(reference["peak_rss_total_mib"])
    )
    lines += [
        "", f"Dotted closure reference: `{warning:g}`.", "",
        "## Comparison with completed 1.35 nm B192 reference", "",
        f"Reference `{reference['tag']}` used {int(reference['macroparticles']):,} "
        f"macroparticles and {int(reference['electrons']):,} electrons. The field-energy "
        "ratio is also divided by the squared charge ratio so the 0.727% charge change "
        "does not masquerade as a numerical field change.", "",
        "| metric | weight-one result | B192 reference | ratio new/reference |",
        "|---|---:|---:|---:|",
        f"| peak abs(R)/E0 | {fmt(ledger['peak_abs_relative_residual'])} | "
        f"{fmt(reference['peak_abs_relative_residual'])} | {fmt(residual_ratio)} |",
        f"| final interior field energy (J) | {fmt(ledger['final_field_energy_J'])} | "
        f"{fmt(reference['final_field_energy_J'])} | {fmt(field_ratio_raw)} |",
        f"| final field energy normalized by Ne^2 | n/a | n/a | "
        f"{fmt(field_ratio_charge_normalized)} |",
        f"| elapsed time (h) | {fmt(resources['last_elapsed_s'] / 3600.0)} | "
        f"{fmt(float(reference['elapsed_s']) / 3600.0)} | {fmt(elapsed_ratio)} |",
        f"| peak total RSS (GiB) | {fmt(resources['peak_rss_total_mib'] / 1024.0)} | "
        f"{fmt(float(reference['peak_rss_total_mib']) / 1024.0)} | {fmt(rss_ratio)} |",
        "",
        f"Runtime growth is `{100.0 * (elapsed_ratio - 1.0):.4g}%`; peak-memory growth is "
        f"`{100.0 * (rss_ratio - 1.0):.4g}%`. A charge-normalized field ratio below 1 "
        "means the unit-weight run further suppresses the interior field relative to B192. "
        "The independent energy-closure criterion remains peak abs(R)/E0 <= 0.05.", "",
        "## Final energy exchange and closure", "",
        "The table exposes the signed ledger terms instead of inferring energy "
        "exchange from mean gamma alone. `PASS` requires every sampled closure record "
        "to be valid and the peak absolute relative residual not to exceed the "
        f"configured `{warning:g}` limit.", "",
        "| case | closure | all records valid | peak abs(R)/E0 | delta K particles (J) | delta Ufield (J) | outward field (J) | prescribed work (J) | final R (J) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in case_results:
        case = result["case"]
        ledger = result["ledger"]
        closure_pass = bool(ledger["closure_pass_at_warning"])
        lines.append(
            f"| {case['tag']} | {'PASS' if closure_pass else 'FAIL'} | "
            f"{str(bool(ledger['closure_valid_all'])).lower()} | "
            f"{fmt(ledger['peak_abs_relative_residual'])} | "
            f"{fmt(ledger['delta_particle_kinetic_total_J'])} | "
            f"{fmt(ledger['delta_field_energy_J'])} | "
            f"{fmt(ledger['outward_field_energy_J'])} | "
            f"{fmt(ledger['prescribed_work_J'])} | "
            f"{fmt(ledger['final_balance_residual_J'])} |"
        )
    lines += [
        "",
        "## Input shot noise", "",
        "The generator diagnostics below record the exact unit-weight B=768 shot-noise "
        "load. With only 8 physical electrons per beamlet, the generator's marginal "
        "small-phase warning is expected and must remain part of the interpretation.", "",
        "| case | electrons/beamlet | global b1 | expected global rms | mean Ne,slice times abs(b1)^2 |",
        "|---|---:|---:|---:|---:|",
    ]
    for result in case_results:
        noise = result["input_noise"]
        lines.append(
            f"| {result['case']['tag']} | "
            f"{fmt(noise.get('input_electrons_per_beamlet', float('nan')))} | "
            f"{fmt(noise.get('input_global_bunching_h1', float('nan')))} | "
            f"{fmt(noise.get('input_expected_global_rms_h1', float('nan')))} | "
            f"{fmt(noise.get('input_scaled_local_power_h1', float('nan')))} |"
        )
    lines += [
        "", "## Microbunching", "",
        "The one-wavelength slice RMS is the primary SASE metric; whole-bunch b can "
        "cancel when slices have unrelated phase. Each particle ID contributes only its "
        "first downstream crossing.", "",
        "| case | plane | unique/expected | duplicate | missing | whole-bunch b1 | 1lambda slice-rms b1 | mean Ne,slice times abs(b1)^2 | mean gamma | sigma gamma |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in case_results:
        for plane_name in ("entry", "exit"):
            plane = result["planes"][plane_name]
            lines.append(
                f"| {result['case']['tag']} | {plane_name} | "
                f"{plane['records']}/{plane['expected_macroparticles']} | "
                f"{plane['duplicate_downstream_records']} | {plane['missing_particle_ids']} | "
                f"{fmt(plane['bunching_h1'])} | "
                f"{fmt(plane['slice_rms_bunching_h1'])} | "
                f"{fmt(plane['mean_electrons_per_slice_times_bunching_squared_h1'])} | "
                f"{fmt(plane['mean_gamma'])} | "
                f"{fmt(plane['sigma_gamma'])} |"
            )
    lines += [
        "", "Treat the comparison as one-seed numerical evidence only; a publishable "
        "macro-particle conclusion needs several independent seeds. Harmonic-2/3 bunching "
        "columns are particle diagnostics only; "
        "the configured radiation-resolution contract covers 1.25 times the fundamental, "
        "not emitted harmonic-2/3 radiation. Automatic z is observation-bounded, so detector/field claims "
        "still require a guard or explicit-box convergence check.", "",
        "Artifacts: `comparison-timeseries.csv`, `case-summary.csv`, "
        "`reference-comparison.csv`, `particle-plane-summary.csv`, `resource-timeseries.csv`, "
        "`comparison-diagnostics.png`, `microbunching.png`, and `summary.json`.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    work = Path(os.environ.get("CAMPAIGN_WORK_DIR", str(root / "work"))).resolve()
    plan = load_plan(work / "metadata" / "campaign-plan.json")
    results_dir = work / "results"
    warning = float(os.environ.get("ENERGY_RESIDUAL_WARNING", "0.05"))
    wavelength = float(plan["theory"]["resonant_wavelength_m"])
    boost_gamma = float(plan["theory"]["boost_gamma"])

    case_results: list[dict[str, Any]] = []
    all_ledger_rows: list[dict[str, Any]] = []
    all_plane_rows: list[dict[str, Any]] = []
    all_resource_rows: list[dict[str, Any]] = []
    case_summary_rows: list[dict[str, Any]] = []

    for case in plan["cases"]:
        # Do not trust absolute paths embedded when the cards were rendered on
        # the cluster.  Results remain analyzable after the work tree is copied.
        run = work / "runs" / case["tag"]
        if not (run / "COMPLETE").is_file():
            continue
        validate_completion(
            run, str(plan["campaign_fingerprint"]), str(case["tag"])
        )
        beamlets = int(case["beamlets_per_wavelength"])
        input_noise = read_input_noise(
            work / "inputs" / f"sase-B{beamlets:04d}.h5",
            run / "input-noise.json",
        )
        manifest = read_manifest(run / "run-manifest.yaml")
        ledger_rows, ledger = read_ledger(
            run / "ledger" / "energy-ledger.h5", case["tag"], boost_gamma
        )
        ledger["closure_warning_limit"] = warning
        ledger["closure_pass_at_warning"] = (
            bool(ledger["closure_valid_all"]) and
            float(ledger["peak_abs_relative_residual"]) <= warning
        )
        entry = read_particle_plane(
            run / "detectors" / "entry-particles.h5", case["tag"], "entry",
            wavelength, int(case["macroparticles"]),
        )
        exit_plane = read_particle_plane(
            run / "detectors" / "exit-particles.h5", case["tag"], "exit",
            wavelength, int(case["macroparticles"]),
        )
        resource_rows, resources = read_resource_log(run / "aprl.log", case["tag"])
        result = {
            "case": case,
            "manifest": manifest,
            "ledger": ledger,
            "ledger_rows": ledger_rows,
            "planes": {"entry": entry, "exit": exit_plane},
            "input_noise": input_noise,
            "resources": resources,
        }
        case_results.append(result)
        all_ledger_rows.extend(ledger_rows)
        all_plane_rows.extend((entry, exit_plane))
        all_resource_rows.extend(resource_rows)
        ledger_flat = {
            key: value for key, value in ledger.items()
            if key != "outward_field_energy_faces_J"
        }
        face_flat = {
            f"outward_{name}_J": value
            for name, value in ledger["outward_field_energy_faces_J"].items()
        }
        case_summary_rows.append({
            "case": case["tag"],
            "dz_lab_nm": float(case["dz_lab_m"]) * 1e9,
            "beamlets_per_wavelength": case["beamlets_per_wavelength"],
            "macroparticles": case["macroparticles"],
            "electrons_per_macroparticle": case["electrons_per_macroparticle"],
            **manifest,
            **ledger_flat,
            **face_flat,
            **input_noise,
            "entry_bunching_h1": entry["bunching_h1"],
            "entry_bunching_h2": entry["bunching_h2"],
            "entry_bunching_h3": entry["bunching_h3"],
            "exit_bunching_h1": exit_plane["bunching_h1"],
            "exit_bunching_h2": exit_plane["bunching_h2"],
            "exit_bunching_h3": exit_plane["bunching_h3"],
            "entry_slice_rms_bunching_h1": entry["slice_rms_bunching_h1"],
            "entry_slice_rms_bunching_h2": entry["slice_rms_bunching_h2"],
            "entry_slice_rms_bunching_h3": entry["slice_rms_bunching_h3"],
            "exit_slice_rms_bunching_h1": exit_plane["slice_rms_bunching_h1"],
            "exit_slice_rms_bunching_h2": exit_plane["slice_rms_bunching_h2"],
            "exit_slice_rms_bunching_h3": exit_plane["slice_rms_bunching_h3"],
            "exit_duplicate_downstream_records": exit_plane["duplicate_downstream_records"],
            "exit_missing_particle_ids": exit_plane["missing_particle_ids"],
            "exit_record_fraction": exit_plane["record_fraction"],
            **resources,
        })

    if not case_results:
        fail(f"no completed cases found under {work / 'runs'}")
    results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(results_dir / "comparison-timeseries.csv", all_ledger_rows)
    write_csv(results_dir / "particle-plane-summary.csv", all_plane_rows)
    write_csv(results_dir / "resource-timeseries.csv", all_resource_rows)
    write_csv(results_dir / "case-summary.csv", case_summary_rows)
    current = case_results[0]
    reference = plan["reference_case"]
    charge_ratio = float(plan["beam"]["electrons"]) / float(reference["electrons"])
    raw_field_ratio = (
        float(current["ledger"]["final_field_energy_J"]) /
        float(reference["final_field_energy_J"])
    )
    reference_comparison = {
        "case": current["case"]["tag"],
        "reference_case": reference["tag"],
        "charge_ratio": charge_ratio,
        "peak_residual_ratio": (
            float(current["ledger"]["peak_abs_relative_residual"]) /
            float(reference["peak_abs_relative_residual"])
        ),
        "final_field_energy_ratio_raw": raw_field_ratio,
        "final_field_energy_ratio_charge_squared_normalized": (
            raw_field_ratio / (charge_ratio * charge_ratio)
        ),
        "elapsed_ratio": (
            float(current["resources"]["last_elapsed_s"]) /
            float(reference["elapsed_s"])
        ),
        "peak_rss_ratio": (
            float(current["resources"]["peak_rss_total_mib"]) /
            float(reference["peak_rss_total_mib"])
        ),
    }
    write_csv(results_dir / "reference-comparison.csv", [reference_comparison])
    make_comparison_plot(
        case_results, results_dir / "comparison-diagnostics.png", warning
    )
    make_bunching_plot(case_results, results_dir / "microbunching.png")
    summary = {
        "scope": "single-seed 4 m SASE exact-unit-weight test at dz=1.35 nm; no K=0 baseline",
        "completed_cases": [result["case"]["tag"] for result in case_results],
        "planned_cases": [case["tag"] for case in plan["cases"]],
        "plan": plan,
        "reference_comparison": reference_comparison,
        "cases": [
            {key: value for key, value in result.items() if key != "ledger_rows"}
            for result in case_results
        ],
    }
    (results_dir / "summary.json").write_text(
        json.dumps(sanitize(summary), ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_report(results_dir / "report.md", plan, case_results, warning)
    print(
        f"Analyzed {len(case_results)}/{len(plan['cases'])} completed cases: "
        f"{results_dir}"
    )


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, KeyError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
