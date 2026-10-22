#!/usr/bin/env python3
"""Summarize APRL field-energy histories without h5py.

The shipped energy_ledger_report executable reads one committed HDF5 record
at a time. This script invokes it for every record by default so that peaks
and positive one-step jumps are exact at the ledger sampling cadence.
"""

import argparse
import bisect
import csv
import json
import re
import subprocess
from pathlib import Path


PAIR = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=([^\s]+)")
OUTSIDE = re.compile(
    r"boosted t=0 transverse ballistic positions outside the physical "
    r"aperture: (\d+)/(\d+) \(([^%]+)%\); outside the full mesh: "
    r"(\d+)/(\d+) \(([^%]+)%\)"
)
LAB_OUTSIDE = re.compile(
    r"reconstructed initial lab snapshot has (\d+) particles outside"
)
ENTRIES = re.compile(
    r"physical-aperture entries=(\d+), currently outside=(\d+)"
)
ACTIVATIONS = re.compile(
    r"Frozen lab-plane injection: activations=(\d+),.*?still frozen=(\d+)"
)
INPUT_RECORDS = re.compile(r"Macro-particle normalization: records=(\d+)")
ALL_PASSED = re.compile(
    r"Configured stop reached: all (\d+) valid particles passed the last "
    r"interaction boundary"
)
ESCAPED = re.compile(r"excluded escaped particles=(\d+)")
FIELDS = (
    "step", "time_box_s", "closure_valid", "active_macroparticles",
    "removed_macroparticles", "boundary_injection_entries",
    "field_energy_interior_J", "outward_field_energy_J",
    "boundary_source_work_J", "boundary_source_abs_work_J",
    "balance_residual_raw_J", "particle_kinetic_active_J",
)


def value(text):
    try:
        return float(text)
    except ValueError:
        return text


def parse_report(text):
    result = {key: value(raw) for key, raw in PAIR.findall(text.splitlines()[0])}
    section = None
    for line in text.splitlines()[1:]:
        match = re.match(r"^(initial|selected|final):", line)
        if match:
            section = match.group(1)
            result[section] = {}
        elif line and not line.startswith(" "):
            section = None
        if section:
            result[section].update(
                (key, value(raw)) for key, raw in PAIR.findall(line)
            )
    return result


def reporter_call(executable, ledger, index=None):
    command = [str(executable), str(ledger)]
    if index is not None:
        command.append(str(index))
    completed = subprocess.run(command, check=False, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    if completed.returncode:
        raise RuntimeError(
            f"{' '.join(command)} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}")
    return parse_report(completed.stdout)


def log_diagnostics(run_dir):
    path = run_dir / "run.log"
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    result = {}
    match = OUTSIDE.search(text)
    if match:
        result["initial_outside_physical_count"] = int(match.group(1))
        result["initial_particles"] = int(match.group(2))
        result["initial_outside_physical_pct"] = float(match.group(3))
        result["initial_outside_mesh_count"] = int(match.group(4))
        result["initial_outside_mesh_pct"] = float(match.group(6))
    match = LAB_OUTSIDE.search(text)
    if match:
        result["initial_lab_outside_physical_count"] = int(match.group(1))
    match = ENTRIES.search(text)
    if match:
        result["physical_aperture_entries"] = int(match.group(1))
        result["final_outside_physical_count"] = int(match.group(2))
    match = ACTIVATIONS.search(text)
    if match:
        result["activated_macroparticles"] = int(match.group(1))
        result["still_frozen_macroparticles"] = int(match.group(2))
    match = INPUT_RECORDS.search(text)
    if match:
        result["input_macroparticles"] = int(match.group(1))
    match = ALL_PASSED.search(text)
    if match:
        result["valid_particles_past_last_element"] = int(match.group(1))
        escaped = ESCAPED.search(text[match.start():].splitlines()[0])
        result["escaped_macroparticles"] = int(escaped.group(1)) if escaped else 0
        result["all_input_particles_past_last_element"] = (
            result["valid_particles_past_last_element"] ==
            result.get("input_macroparticles") and
            result["escaped_macroparticles"] == 0 and
            result.get("still_frozen_macroparticles") == 0
        )
    if "no valid particles remain before the last interaction" in text:
        result["all_input_particles_past_last_element"] = False
    result["duration_exhausted"] = "mesh.duration was exhausted" in text
    result["normal_stop"] = "Configured stop reached:" in text
    return result


def inspect_run(run_dir, executable, stride):
    ledger = run_dir / "ledger" / "energy-ledger.h5"
    first = reporter_call(executable, ledger)
    count = int(first["committed_records"])
    indices = sorted(set([0, count - 1] + list(range(1, count - 1, stride))))
    by_index = {0: first["initial"], count - 1: first["final"]}
    for index in indices:
        if index not in by_index:
            by_index[index] = reporter_call(executable, ledger, index)["selected"]
    trace = [{"record_index": index, **by_index[index]} for index in indices]
    energies = [row["field_energy_interior_J"] for row in trace]
    jumps = [right - left for left, right in zip(energies, energies[1:])]
    positive = max([0.0] + jumps)
    initial = trace[0]
    final = trace[-1]
    summary = {
        "run_dir": str(run_dir),
        "ledger": str(ledger),
        "format_version": int(first["format_version"]),
        "complete": bool(first["complete"]),
        "committed_records": count,
        "scanned_records": len(trace),
        "scan_stride": stride,
        "initial_field_J": energies[0],
        "peak_field_J": max(energies),
        "peak_field_step": int(trace[energies.index(max(energies))]["step"]),
        "final_field_J": energies[-1],
        "max_positive_sample_jump_J": positive,
        "max_positive_sample_jump_over_initial_field":
            positive / energies[0] if energies[0] > 0 else None,
        "final_field_over_initial_field":
            energies[-1] / energies[0] if energies[0] > 0 else None,
        "final_outward_field_energy_J": final["outward_field_energy_J"],
        "final_boundary_source_abs_work_J":
            final.get("boundary_source_abs_work_J"),
        "final_boundary_source_work_J": final["boundary_source_work_J"],
        "final_raw_residual_J": final["balance_residual_raw_J"],
        "initial_active_macroparticles": int(initial["active_macroparticles"]),
        "final_active_macroparticles": int(final["active_macroparticles"]),
        "closure_valid_all_scanned":
            all(row["closure_valid"] == 1 for row in trace),
        "elapsed_box_time_s": final["time_box_s"],
        **log_diagnostics(run_dir),
    }
    return summary, trace


def interpolate(trace, time):
    times = [row["time_box_s"] for row in trace]
    at = bisect.bisect_left(times, time)
    if at == 0:
        return trace[0]["field_energy_interior_J"]
    if at >= len(trace):
        return trace[-1]["field_energy_interior_J"]
    left, right = trace[at - 1], trace[at]
    a, b = left["time_box_s"], right["time_box_s"]
    if b <= a:
        raise ValueError("ledger record times must increase strictly")
    fraction = (time - a) / (b - a)
    return (left["field_energy_interior_J"] + fraction *
            (right["field_energy_interior_J"] -
             left["field_energy_interior_J"]))


def compare(left_summary, left_trace, right_summary, right_trace):
    common_end = min(left_trace[-1]["time_box_s"],
                     right_trace[-1]["time_box_s"])
    times = sorted(set([0.0, common_end] + [
        row["time_box_s"] for trace in (left_trace, right_trace)
        for row in trace if row["time_box_s"] <= common_end
    ]))
    differences = [interpolate(left_trace, time) -
                   interpolate(right_trace, time) for time in times]
    scale = max(left_summary["peak_field_J"],
                right_summary["peak_field_J"])
    return {
        "left": Path(left_summary["run_dir"]).name,
        "right": Path(right_summary["run_dir"]).name,
        "last_common_time_s": common_end,
        "field_difference_at_last_common_time_J": differences[-1],
        "max_abs_field_difference_J": max(map(abs, differences)),
        "max_abs_field_difference_over_larger_peak":
            max(map(abs, differences)) / scale if scale > 0 else None,
        "matching_method": "piecewise-linear interpolation in box time",
    }


def print_summary(summary):
    name = Path(summary["run_dir"]).name
    print(f"{name}: records={summary['committed_records']} "
          f"complete={int(summary['complete'])} "
          f"normal_stop={summary.get('normal_stop', 'unknown')}")
    print("  U_field [J] initial / peak / final: "
          f"{summary['initial_field_J']:.9e} / "
          f"{summary['peak_field_J']:.9e} / "
          f"{summary['final_field_J']:.9e}")
    print("  max positive sampled ΔU [J]: "
          f"{summary['max_positive_sample_jump_J']:.9e}")
    if "initial_outside_physical_pct" in summary:
        print("  initial outside physical / mesh: "
              f"{summary['initial_outside_physical_pct']:.6g}% / "
              f"{summary['initial_outside_mesh_pct']:.6g}%")
    if "still_frozen_macroparticles" in summary:
        print("  activated / still frozen: "
              f"{summary['activated_macroparticles']} / "
              f"{summary['still_frozen_macroparticles']}")
    if "all_input_particles_past_last_element" in summary:
        print("  all input particles past last element: "
              f"{summary['all_input_particles_past_last_element']}")
        if not summary["all_input_particles_past_last_element"]:
            print("  WARNING: final state is not a complete whole-bunch "
                  "crossing of the last element")
    print("  final boundary-source |work| [J]: "
          f"{summary['final_boundary_source_abs_work_J']}")
    if not summary["closure_valid_all_scanned"]:
        print("  closure_valid=0: raw residual is a diagnostic, not a "
              "whole-bunch conservation proof")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path,
                        help="run directories containing ledger/energy-ledger.h5")
    parser.add_argument("--reporter", type=Path,
                        default=Path(__file__).resolve().parents[2] /
                        "build/energy_ledger_report")
    parser.add_argument("--stride", type=int, default=1,
                        help="scan every Nth interior record; N>1 makes "
                        "peak and jump estimates approximate")
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--csv-out", type=Path,
                        help="optional sampled energy-history table")
    args = parser.parse_args()
    if args.stride < 1:
        parser.error("--stride must be positive")
    if not args.reporter.is_file():
        parser.error("energy_ledger_report is missing; build it or pass --reporter")
    results = [inspect_run(run.resolve(), args.reporter.resolve(), args.stride)
               for run in args.runs]
    for summary, _ in results:
        print_summary(summary)
    comparisons = []
    for index in range(1, len(results)):
        comp = compare(*results[0], *results[index])
        comparisons.append(comp)
        print(f"{comp['left']} − {comp['right']} at matched box time "
              f"{comp['last_common_time_s']:.6e} s: "
              f"ΔU={comp['field_difference_at_last_common_time_J']:.9e} J; "
              "max |ΔU| / larger peak="
              f"{comp['max_abs_field_difference_over_larger_peak']:.6g}")
    if args.csv_out:
        args.csv_out.parent.mkdir(parents=True, exist_ok=True)
        with args.csv_out.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["run", "record_index", *FIELDS])
            for summary, trace in results:
                for row in trace:
                    writer.writerow([Path(summary["run_dir"]).name,
                                     row["record_index"],
                                     *(row.get(field, "") for field in FIELDS)])
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "runs": [summary for summary, _ in results],
            "comparisons": comparisons,
            "caveat": "Frozen/injection closure_valid=0; compare energy "
                      "histories as diagnostics, not exact conservation.",
        }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
