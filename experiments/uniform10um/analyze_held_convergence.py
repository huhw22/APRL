#!/usr/bin/env python3
"""Compare held-start free-drift runs on one fixed laboratory field aperture.

The field plane stores the *total* Maxwell field, including particle near field.
Its E and Yee-half-time B are read one sample at a time.  The boosted-frame
energy ledger is reported separately and is not treated as a laboratory flux.
"""

import argparse
import csv
import json
import math
from contextlib import ExitStack
from pathlib import Path

import h5py
import numpy as np


C = 299792458.0
MU0 = 4.0e-7 * math.pi


def _scalar(group, name):
    if name not in group:
        raise ValueError(f"{group.name}: missing {name}")
    return int(group[name][()])


def _times(dataset, count, label):
    if count < 2 or count > len(dataset):
        raise ValueError(f"{label}: invalid committed count {count}")
    values = np.asarray(dataset[:count], dtype=np.float64)
    if not np.all(np.isfinite(values)) or not np.all(np.diff(values) > 0):
        raise ValueError(f"{label}: times must be finite and strictly increasing")
    return values


def _ratio(numerator, denominator):
    return float(numerator / denominator) if denominator > 0 else None


def _format_optional(value):
    return f"{value:.6g}" if value is not None else "undefined"


def _integral(values, times):
    return float(np.sum((values[:-1] + values[1:]) * (np.diff(times) / 2),
                        dtype=np.float64))


def _weights(times):
    weights = np.empty(len(times), dtype=np.float64)
    weights[0] = (times[1] - times[0]) / 2
    weights[-1] = (times[-1] - times[-2]) / 2
    if len(times) > 2:
        weights[1:-1] = (times[2:] - times[:-2]) / 2
    return weights


def _coordinates(group, axis):
    first = float(group.attrs[f"{axis}_first_m"])
    step = float(group.attrs[f"d{axis}_m"])
    count = int(group.attrs[f"n{axis}"])
    if not math.isfinite(step) or step <= 0 or count < 1:
        raise ValueError(f"{group.name}: invalid {axis} grid")
    return first + np.arange(count, dtype=np.float64) * step, step


def _grid_slice(source, target, step, label):
    """Find the same cell centres, even when a variant widens its box."""
    if len(target) < 1:
        raise ValueError(f"{label}: empty ROI")
    indices = np.rint((target - source[0]) / step).astype(np.int64)
    tolerance = max(abs(step) * 1e-7, 1e-12)
    if (np.any(indices < 0) or np.any(indices >= len(source)) or
            np.any(np.diff(indices) != 1)):
        raise ValueError(f"{label}: fixed ROI is not contained in variant grid")
    if not np.allclose(source[indices], target, rtol=0, atol=tolerance):
        raise ValueError(f"{label}: cell centres do not align; spatial resampling is not implicit")
    return slice(int(indices[0]), int(indices[-1]) + 1)


class Plane:
    def __init__(self, stack, name, run_dir, field_name):
        self.name = name
        self.run_dir = Path(run_dir).resolve()
        path = self.run_dir / "detectors" / (field_name + ".h5")
        self.file = stack.enter_context(h5py.File(path, "r"))
        self.group = self.file["field_plane"]
        self.complete = bool(_scalar(self.group, "complete"))
        if not self.complete:
            raise ValueError(f"{path}: field plane is incomplete")
        self.count = _scalar(self.group, "committed_samples")
        self.times = _times(self.group["time_s"], self.count, str(path))
        self.fields = self.group["fields"]
        if self.fields.ndim != 3 or self.fields.shape[0] < self.count:
            raise ValueError(f"{path}: invalid field dimensions")
        if not {"electric_V_per_m", "magnetic_T"}.issubset(self.fields.dtype.names or ()):
            raise ValueError(f"{path}: missing E/B components")
        self.x, self.dx = _coordinates(self.group, "x")
        self.y, self.dy = _coordinates(self.group, "y")
        if self.fields.shape[1:] != (len(self.y), len(self.x)):
            raise ValueError(f"{path}: field dimensions disagree with grid attributes")
        self.z = float(self.group.attrs["plane_z_m"])
        self.x_slice = slice(0, len(self.x))
        self.y_slice = slice(0, len(self.y))

    def select_roi(self, x_target, y_target, dx, dy):
        if not math.isclose(self.dx, dx, rel_tol=1e-7, abs_tol=1e-12):
            raise ValueError(f"{self.name}: dx differs; choose a common grid")
        if not math.isclose(self.dy, dy, rel_tol=1e-7, abs_tol=1e-12):
            raise ValueError(f"{self.name}: dy differs; choose a common grid")
        self.x_slice = _grid_slice(self.x, x_target, self.dx, self.name + "/x")
        self.y_slice = _grid_slice(self.y, y_target, self.dy, self.name + "/y")

    def frame(self, index):
        raw = self.fields[index, self.y_slice, self.x_slice]
        return (np.asarray(raw["electric_V_per_m"], dtype=np.float64),
                np.asarray(raw["magnetic_T"], dtype=np.float64))

    def at(self, time):
        right = int(np.searchsorted(self.times, time, side="left"))
        if right < self.count and self.times[right] == time:
            return self.frame(right)
        if right == 0 or right == self.count:
            raise ValueError(f"{self.name}: time outside committed plane window")
        left = right - 1
        alpha = (time - self.times[left]) / (self.times[right] - self.times[left])
        e0, b0 = self.frame(left)
        e1, b1 = self.frame(right)
        return ((1 - alpha) * e0 + alpha * e1,
                (1 - alpha) * b0 + alpha * b1)


def _power(e, b, area):
    sz = (e[..., 0] * b[..., 1] - e[..., 1] * b[..., 0]) / MU0
    return (float(np.sum(sz, dtype=np.float64) * area),
            float(np.sum(np.maximum(sz, 0), dtype=np.float64) * area))


def _edge_summary(signed, forward, threshold):
    peak_signed = float(np.max(np.abs(signed)))
    peak_forward = float(np.max(forward))
    ratios = {
        "signed_start_over_peak_abs": _ratio(abs(signed[0]), peak_signed),
        "signed_end_over_peak_abs": _ratio(abs(signed[-1]), peak_signed),
        "forward_start_over_peak": _ratio(forward[0], peak_forward),
        "forward_end_over_peak": _ratio(forward[-1], peak_forward),
    }
    ratios["tail_unresolved"] = any(
        ratios[key] is not None and ratios[key] > threshold
        for key in ("signed_end_over_peak_abs", "forward_end_over_peak")
    )
    ratios["head_unresolved"] = any(
        ratios[key] is not None and ratios[key] > threshold
        for key in ("signed_start_over_peak_abs", "forward_start_over_peak")
    )
    return ratios


def _native_power(plane, area, threshold):
    signed = np.empty(plane.count, dtype=np.float64)
    forward = np.empty(plane.count, dtype=np.float64)
    for i in range(plane.count):
        signed[i], forward[i] = _power(*plane.frame(i), area)
    return {
        "samples": plane.count,
        "first_lab_time_s": float(plane.times[0]),
        "last_lab_time_s": float(plane.times[-1]),
        "signed_energy_J": _integral(signed, plane.times),
        "forward_energy_J": _integral(forward, plane.times),
        "edge": _edge_summary(signed, forward, threshold),
    }


def _ledger(run_dir):
    path = Path(run_dir) / "ledger" / "energy-ledger.h5"
    with h5py.File(path, "r") as file:
        group = file["energy_ledger"]
        if not _scalar(group, "complete"):
            raise ValueError(f"{path}: ledger is incomplete")
        count = _scalar(group, "committed_records")
        dataset = group["records"]
        required = ("time_box_s", "field_energy_interior_J", "outward_field_energy_J")
        if not set(required).issubset(dataset.dtype.names or ()):
            raise ValueError(f"{path}: missing required ledger columns")
        if count < 2 or count > len(dataset):
            raise ValueError(f"{path}: invalid committed record count")
        # Select only scalar columns, not particle-moment or face-vector payloads.
        records = dataset.fields(required)[:count]
        values = {key: np.asarray(records[key], dtype=np.float64) for key in required}
        times = values["time_box_s"]
        if not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0):
            raise ValueError(f"{path}: invalid box-time sequence")
        if not all(np.all(np.isfinite(values[key])) for key in required[1:]):
            raise ValueError(f"{path}: nonfinite ledger energy")
        summary = {
            "records": count,
            "first_box_time_s": float(times[0]),
            "last_box_time_s": float(times[-1]),
            "field_energy_initial_J": float(values["field_energy_interior_J"][0]),
            "field_energy_peak_J": float(np.max(values["field_energy_interior_J"])),
            "field_energy_final_J": float(values["field_energy_interior_J"][-1]),
            "outward_field_energy_final_J": float(values["outward_field_energy_J"][-1]),
        }
        return values, summary


def _named_path(raw):
    if "=" not in raw:
        raise argparse.ArgumentTypeError("use NAME=RUN_DIR")
    name, path = raw.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("use NAME=RUN_DIR")
    return name, Path(path)


def analyze(reference, variants, roi_half_width_um=160.0, field_name="near-field",
            edge_warning_fraction=0.05, time_start_s=None, time_end_s=None):
    if not variants:
        raise ValueError("at least one --variant is required")
    names = [reference[0]] + [name for name, _ in variants]
    if len(set(names)) != len(names):
        raise ValueError("run names must be unique")
    if (not math.isfinite(roi_half_width_um) or roi_half_width_um <= 0 or
            not math.isfinite(edge_warning_fraction) or
            not 0 <= edge_warning_fraction <= 1):
        raise ValueError("ROI width must be positive and edge fraction in [0,1]")
    if (time_start_s is not None and not math.isfinite(time_start_s) or
            time_end_s is not None and not math.isfinite(time_end_s)):
        raise ValueError("requested lab-time bounds must be finite")
    with ExitStack() as stack:
        runs = [Plane(stack, name, path, field_name) for name, path in [reference] + variants]
        ref = runs[0]
        half_width_m = roi_half_width_um * 1e-6
        x_target = ref.x[np.abs(ref.x) <= half_width_m + 1e-12]
        y_target = ref.y[np.abs(ref.y) <= half_width_m + 1e-12]
        if len(x_target) < 2 or len(y_target) < 2:
            raise ValueError("ROI has fewer than two cells on an axis")
        for run in runs:
            if not math.isclose(run.z, ref.z, abs_tol=1e-9):
                raise ValueError(f"{run.name}: field-plane z differs from reference")
            run.select_roi(x_target, y_target, ref.dx, ref.dy)
        area = ref.dx * ref.dy
        native = {run.name: _native_power(run, area, edge_warning_fraction) for run in runs}
        start = max(run.times[0] for run in runs)
        end = min(run.times[-1] for run in runs)
        if time_start_s is not None:
            start = max(start, time_start_s)
        if time_end_s is not None:
            end = min(end, time_end_s)
        common_times = ref.times[(ref.times >= start) & (ref.times <= end)]
        if len(common_times) < 2:
            raise ValueError("fewer than two reference samples in common lab-time interval")
        weights = _weights(common_times)
        powers = {run.name: {"signed": [], "forward": []} for run in runs}
        norms = {name: {"e_diff": 0.0, "b_diff": 0.0, "em_diff": 0.0}
                 for name in names[1:]}
        ref_e_norm = ref_b_norm = ref_em_norm = 0.0
        trace = []
        for time, weight in zip(common_times, weights):
            e_ref, b_ref = ref.at(float(time))
            e_norm = float(np.sum(e_ref * e_ref, dtype=np.float64))
            b_norm = float(np.sum(b_ref * b_ref, dtype=np.float64))
            ref_e_norm += weight * e_norm
            ref_b_norm += weight * b_norm
            ref_em_norm += weight * (e_norm + C * C * b_norm)
            for run in runs:
                e, b = (e_ref, b_ref) if run is ref else run.at(float(time))
                signed, forward = _power(e, b, area)
                powers[run.name]["signed"].append(signed)
                powers[run.name]["forward"].append(forward)
                trace.append({"kind": "lab_plane", "run": run.name,
                              "time_lab_s": float(time), "time_box_s": "",
                              "signed_power_W": signed, "forward_power_W": forward,
                              "field_energy_interior_J": "", "outward_field_energy_J": ""})
                if run is not ref:
                    de = e - e_ref
                    db = b - b_ref
                    e_diff = float(np.sum(de * de, dtype=np.float64))
                    b_diff = float(np.sum(db * db, dtype=np.float64))
                    norms[run.name]["e_diff"] += weight * e_diff
                    norms[run.name]["b_diff"] += weight * b_diff
                    norms[run.name]["em_diff"] += weight * (e_diff + C * C * b_diff)
        lab = {}
        for run in runs:
            name = run.name
            signed = np.asarray(powers[name]["signed"])
            forward = np.asarray(powers[name]["forward"])
            lab[name] = {
                "native": native[name],
                "matched_signed_energy_J": _integral(signed, common_times),
                "matched_forward_energy_J": _integral(forward, common_times),
                "matched_edge": _edge_summary(signed, forward, edge_warning_fraction),
            }
        ledgers = {}
        for run in runs:
            values, summary = _ledger(run.run_dir)
            ledgers[run.name] = {"values": values, "summary": summary}
        box_start = max(item["values"]["time_box_s"][0] for item in ledgers.values())
        box_end = min(item["values"]["time_box_s"][-1] for item in ledgers.values())
        if box_end <= box_start:
            raise ValueError("ledger box-time windows do not overlap")
        for run in runs:
            name = run.name
            values = ledgers[name]["values"]
            times = values["time_box_s"]
            summary = ledgers[name]["summary"]
            summary["common_box_time_s"] = float(box_end)
            summary["common_field_energy_J"] = float(np.interp(
                box_end, times, values["field_energy_interior_J"]))
            summary["common_outward_field_energy_J"] = float(np.interp(
                box_end, times, values["outward_field_energy_J"]))
            for i in range(len(times)):
                trace.append({"kind": "boosted_ledger", "run": name,
                              "time_lab_s": "", "time_box_s": float(times[i]),
                              "signed_power_W": "", "forward_power_W": "",
                              "field_energy_interior_J": float(values["field_energy_interior_J"][i]),
                              "outward_field_energy_J": float(values["outward_field_energy_J"][i])})
        comparisons = {}
        ref_forward = lab[ref.name]["matched_forward_energy_J"]
        ref_signed = lab[ref.name]["matched_signed_energy_J"]
        ref_ledger = ledgers[ref.name]["summary"]
        for run in runs[1:]:
            name = run.name
            differences = norms[name]
            field_tail_unresolved = bool(
                lab[name]["native"]["edge"]["tail_unresolved"]
                or lab[ref.name]["native"]["edge"]["tail_unresolved"]
                or lab[name]["matched_edge"]["tail_unresolved"]
                or lab[ref.name]["matched_edge"]["tail_unresolved"])
            comparisons[name] = {
                "relative_E_spacetime_L2": _ratio(math.sqrt(differences["e_diff"]),
                                                    math.sqrt(ref_e_norm)),
                "relative_B_spacetime_L2": _ratio(math.sqrt(differences["b_diff"]),
                                                    math.sqrt(ref_b_norm)),
                "relative_combined_E_cB_spacetime_L2": _ratio(
                    math.sqrt(differences["em_diff"]), math.sqrt(ref_em_norm)),
                "matched_signed_energy_difference_J": lab[name]["matched_signed_energy_J"] - ref_signed,
                "matched_forward_energy_difference_J": lab[name]["matched_forward_energy_J"] - ref_forward,
                "matched_forward_energy_relative_difference": (
                    (lab[name]["matched_forward_energy_J"] - ref_forward) / ref_forward
                    if ref_forward != 0 else None),
                "common_box_field_energy_difference_J": (
                    ledgers[name]["summary"]["common_field_energy_J"]
                    - ref_ledger["common_field_energy_J"]),
                "common_box_outward_field_energy_difference_J": (
                    ledgers[name]["summary"]["common_outward_field_energy_J"]
                    - ref_ledger["common_outward_field_energy_J"]),
                "tail_unresolved": field_tail_unresolved,
                "transmitted_energy_status": (
                    "inconclusive_tail_unresolved" if field_tail_unresolved
                    else "recorded_window_only_not_radiation_isolated"),
            }
        report = {
            "reference": ref.name,
            "plane_z_m": ref.z,
            "fixed_roi": {"requested_half_width_um": roi_half_width_um,
                          "x_first_m": float(x_target[0]), "x_last_m": float(x_target[-1]),
                          "y_first_m": float(y_target[0]), "y_last_m": float(y_target[-1]),
                          "nx": len(x_target), "ny": len(y_target),
                          "dx_m": ref.dx, "dy_m": ref.dy},
            "common_lab_window": {"start_s": float(common_times[0]),
                                  "end_s": float(common_times[-1]),
                                  "reference_samples": len(common_times)},
            "common_box_window": {"start_s": float(box_start), "end_s": float(box_end)},
            "runs": {run.name: {"path": str(run.run_dir), "lab_plane": lab[run.name],
                                "boosted_ledger": ledgers[run.name]["summary"]}
                     for run in runs},
            "comparisons": comparisons,
            "caveats": [
                "The raw field plane is total Maxwell field including particle near field, not isolated radiation.",
                "B is stored at Yee leapfrog half time; E x B has a time-stagger discretization error, especially when dz/dt differ.",
                "Ledger U_field is boosted-frame interior-volume energy. Different transverse box volumes make absolute U_field an unsuitable convergence criterion.",
                "A non-decayed field-plane tail makes transmitted-energy convergence inconclusive; extend the observation window.",
            ],
        }
        return report, trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=_named_path, metavar="NAME=RUN_DIR")
    parser.add_argument("--variant", required=True, action="append", type=_named_path,
                        metavar="NAME=RUN_DIR")
    parser.add_argument("--roi-half-width-um", type=float, default=160.0)
    parser.add_argument("--field-name", default="near-field")
    parser.add_argument("--edge-warning-fraction", type=float, default=0.05)
    parser.add_argument("--time-start-s", type=float)
    parser.add_argument("--time-end-s", type=float)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--csv-out", type=Path, required=True)
    args = parser.parse_args()
    try:
        report, trace = analyze(args.reference, args.variant, args.roi_half_width_um,
                                args.field_name, args.edge_warning_fraction,
                                args.time_start_s, args.time_end_s)
    except (OSError, KeyError, ValueError) as exc:
        parser.error(str(exc))
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    with args.csv_out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(trace[0]))
        writer.writeheader()
        writer.writerows(trace)
    print(f"reference={report['reference']} ROI={report['fixed_roi']['nx']}x{report['fixed_roi']['ny']} "
          f"lab samples={report['common_lab_window']['reference_samples']}")
    for name, comparison in report["comparisons"].items():
        print(f"{name}: E_L2={_format_optional(comparison['relative_E_spacetime_L2'])} "
              f"B_L2={_format_optional(comparison['relative_B_spacetime_L2'])} "
              f"forward_energy_delta={comparison['matched_forward_energy_difference_J']:.6g} J "
              f"tail_unresolved={comparison['tail_unresolved']}")
    print(f"JSON: {args.json_out}\nCSV: {args.csv_out}")


if __name__ == "__main__":
    main()
