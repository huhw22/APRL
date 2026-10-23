#!/usr/bin/env python3
"""Summarize frozen generic-particle reports without pooling physical cases or modes."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


CPU_BUDGET = 128
MIN_THROUGHPUT_REPEATS = 3
STAGES = ("field_total", "particle_prepare", "particle_deposit", "particle_migration_total",
          "diagnostics", "control", "field_communication", "particle_communication")
NESTED = {"field_communication": "field_total", "particle_communication": "particle_migration_total"}
METRICS = ("reported_wall_s", "loop_wall_s", "process_wall_s", "loop_seconds_per_step",
           "logical_particle_updates_per_s", "logical_grid_cell_updates_per_s")


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def number(value, positive=False):
    result = float(value)
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError("expected finite " + ("positive" if positive else "nonnegative") + " measurement")
    return result


def stats(values):
    return {"min": min(values), "median": statistics.median(values), "max": max(values)}


def profile_summary(run):
    """APRL reportPerformance emits interval totals, then resets each stage timer."""
    stages, workload, initialization = {}, [], []
    for record in run.get("performance", []):
        phase = record["phase"]
        if phase == "initialization":
            initialization.append({key: number(record[key]) for key in ("mean_s", "max_s")})
        elif phase in ("interval", "final"):
            stage = record["stage"]
            if stage not in STAGES:
                raise ValueError("unknown profiling stage: " + stage)
            row = {"phase": phase, "step": int(record["step"]), "steps": int(record["steps"]),
                   "elapsed_s": number(record["elapsed_s"], True),
                   "mean_s": number(record["mean_s"]), "max_s": number(record["max_s"]),
                   "max_rank": int(record["max_rank"])}
            if row["steps"] < 1 or not 0 <= row["max_rank"] < run["ranks"]:
                raise ValueError("invalid profiling interval length or rank")
            stages.setdefault(stage, []).append(row)
        elif phase == "workload":
            keys = ("step", "ranks", "threads", "particles_total", "particles_min", "particles_max",
                    "active_total", "active_min", "active_max", "frozen_total", "frozen_min", "frozen_max")
            row = {key: int(record[key]) for key in keys}
            if row["ranks"] != run["ranks"] or row["threads"] != run["threads"] or min(row.values()) < 0:
                raise ValueError("workload rank/thread count differs from launch")
            for prefix in ("particles", "active", "frozen"):
                low, total, high = (row[prefix + suffix] for suffix in ("_min", "_total", "_max"))
                if not low * row["ranks"] <= total <= high * row["ranks"]:
                    raise ValueError("invalid workload rank range")
            workload.append(row)
        else:
            raise ValueError("unknown profiling phase: " + phase)
    result = {"run": run["name"], "profile_interval": run["profile_interval"],
              "initialization": initialization, "stages": [], "workload_samples": workload}
    for stage, rows in stages.items():
        result["stages"].append({
            "stage": stage, "nested_in": NESTED.get(stage), "interval_count": len(rows),
            "sampled_steps": sum(row["steps"] for row in rows),
            "interval_rank_mean_s": stats([row["mean_s"] for row in rows]),
            "interval_rank_max_s": stats([row["max_s"] for row in rows]),
            "rank_mean_s_per_sampled_step": sum(row["mean_s"] for row in rows) / sum(row["steps"] for row in rows),
            "interval_rank_max_s_per_step": stats([row["max_s"] / row["steps"] for row in rows]),
            "intervals": rows})
    if workload:
        result["workload_summary"] = {
            "samples": len(workload),
            "active_particles_per_rank_observed_min": min(row["active_min"] for row in workload),
            "active_particles_per_rank_observed_max": max(row["active_max"] for row in workload),
            "active_total": stats([row["active_total"] for row in workload]),
            "particles_per_rank_observed_min": min(row["particles_min"] for row in workload),
            "particles_per_rank_observed_max": max(row["particles_max"] for row in workload),
            "frozen_total": stats([row["frozen_total"] for row in workload])}
    return result


def best_tested_node_comparison(table, allocation, throughput_eligible):
    baseline = [row for row in table if row["label"] == "baseline"]
    candidate = [row for row in table if row["label"] == "candidate"]
    result = {
        "eligible": False,
        "selection": "best tested layouts on one exclusive node; combines implementation and layout selection",
        "selection_metric": "minimum median APRL total time (reported_wall_s) within each executable",
        "allocated_node_cpus": allocation.get("allocated_cpus"),
        "same_selected_layouts_for_all_timing_metrics": True,
        "baseline": None, "candidate": None,
        "total_speedup": None, "loop_speedup": None, "process_wall_speedup": None,
        "exclusion_reasons": []}
    if not throughput_eligible:
        result["exclusion_reasons"].append("complete, passed, repeated throughput evidence is required")
    if allocation.get("nodes") != 1 or allocation.get("allocated_cpus") != CPU_BUDGET or \
            allocation.get("exclusive_requested_by_batch_script") is not True:
        result["exclusion_reasons"].append("comparison requires the recorded exclusive allocation of one 128-CPU node")
    if not baseline or not candidate:
        result["exclusion_reasons"].append("both baseline and candidate measured layouts are required")
    if any(row["completed_repeats"] < MIN_THROUGHPUT_REPEATS for row in baseline + candidate):
        result["exclusion_reasons"].append("every tested baseline and candidate layout requires at least three completed repeats")
    if result["exclusion_reasons"]:
        return result
    for label, choices in (("baseline", baseline), ("candidate", candidate)):
        selected = min(choices, key=lambda row: (row["statistics"]["reported_wall_s"]["median"],
                                                row["ranks"], row["threads"]))
        result[label] = {key: selected[key] for key in ("ranks", "threads", "step_cpus", "completed_repeats")}
        result[label].update(
            total_median_s=selected["statistics"]["reported_wall_s"]["median"],
            loop_median_s=selected["statistics"]["loop_wall_s"]["median"],
            process_wall_median_s=selected["statistics"]["process_wall_s"]["median"])
    result["eligible"] = True
    for metric in ("total", "loop", "process_wall"):
        result[metric + "_speedup"] = result["baseline"][metric + "_median_s"] / result["candidate"][metric + "_median_s"]
    return result


def analyze(report, plan):
    runs, comparisons = report.get("runs", []), report.get("comparisons", [])
    errors = []
    mode = report.get("mode")
    if mode not in ("calibration", "throughput", "profile"):
        errors.append("unrecognized report mode")
    expected = plan.get("runs", []) if plan else []
    names = [run.get("name") for run in runs]
    expected_names = [run.get("name") for run in expected]
    if not plan:
        errors.append("adjacent plan.json is missing; planned completion cannot be verified")
    elif plan.get("card_sha256") != report.get("card_sha256") or plan.get("mode") != mode:
        errors.append("plan and report physical card hash/mode differ")
    if len(set(names)) != len(names):
        errors.append("duplicate run names")
    if any(run.get("card_sha256") != report.get("card_sha256") for run in runs):
        errors.append("run physical card hash differs from report")
    planned = {run["name"]: run for run in expected}
    for run in runs:
        reference = planned.get(run.get("name"))
        if reference is None or any(run.get(key) != reference.get(key) for key in
                                    ("label", "ranks", "threads", "repeat", "profile_interval", "step_cpus")):
            errors.append("run layout/repeat/profile differs from plan: " + str(run.get("name")))
    validated = [run for run in runs if run.get("valid") is True and run.get("returncode") == 0
                 and run.get("status") == "passed"]
    complete = bool(expected) and set(names) == set(expected_names) and len(validated) == len(expected)
    covered = {row.get("candidate") for row in comparisons if row.get("mode") == "same-rank-bitwise"}
    same_rank = bool(comparisons) and all(row.get("pass") is True for row in comparisons) and covered == set(expected_names)
    probes = report.get("probes", [])
    probe_layouts = {(row.get("ranks"), row.get("threads")) for row in probes}
    expected_probes = {(row.get("ranks"), row.get("threads")) for row in plan.get("probes", [])} if plan else set()
    bindings = bool(probes) and probe_layouts == expected_probes and all(
        row.get("status") == "passed" and row.get("returncode") == 0 for row in probes)
    cross = report.get("cross_rank_comparisons", [])
    cross_result = None if not cross else all(row.get("pass") is True for row in cross)
    gate_pass = complete and same_rank and bindings and report.get("pass") is True and \
        report.get("status") == "passed" and report.get("card_unchanged") is True and \
        report.get("implementation_equivalence_pass") is True and not errors
    failed = report.get("status") == "failed" or bool(report.get("error")) or any(
        row.get("returncode") not in (None, 0) for row in runs + probes)
    status = "passed" if gate_pass else ("failed" if failed else "pending_or_incomplete")
    if errors:
        status = "invalid"
    groups = {}
    for run in validated:
        groups.setdefault((run["label"], run["ranks"], run["threads"]), []).append(run)
    table, profiles = [], []
    for (label, ranks, threads), rows in groups.items():
        try:
            measurements = {key: stats([number(row[key], True) for row in rows]) for key in METRICS}
            table.append({"label": label, "ranks": ranks, "threads": threads, "step_cpus": ranks * threads,
                          "completed_repeats": len(rows), "steps": sorted({row["steps"] for row in rows}),
                          "statistics": measurements})
        except (ValueError, KeyError, TypeError) as exc:
            errors.append("invalid timing measurement: " + str(exc))
    fixed = [row for row in table if row["step_cpus"] == CPU_BUDGET]
    reference = next((row for row in fixed if row["label"] == "baseline" and row["ranks"] == CPU_BUDGET), None)
    profile_off = bool(runs) and all(row.get("profile_interval") == 0 for row in runs)
    budget_valid = all(run.get("step_cpus") == run.get("ranks", 0) * run.get("threads", 0)
                       and 0 < run["step_cpus"] <= CPU_BUDGET for run in runs) and all(
        run.get("step_cpus") == CPU_BUDGET for run in runs if run.get("label") == "candidate")
    repeats_sufficient = bool(fixed) and all(row["completed_repeats"] >= MIN_THROUGHPUT_REPEATS for row in fixed)
    eligible = bool(gate_pass and mode == "throughput" and profile_off and budget_valid and reference
                    and repeats_sufficient and not errors)
    for row in fixed:
        row["fixed_budget_total_speedup"] = None
        row["fixed_budget_loop_speedup"] = None
        if eligible:
            for key, metric in (("fixed_budget_total_speedup", "reported_wall_s"),
                                ("fixed_budget_loop_speedup", "loop_wall_s")):
                row[key] = reference["statistics"][metric]["median"] / row["statistics"][metric]["median"]
    if mode == "profile":
        for run in validated:
            if run.get("profile_interval", 0) > 0:
                try:
                    profile = profile_summary(run)
                    if not profile["stages"] or not profile["workload_samples"]:
                        raise ValueError("enabled profiling has no stage/workload records")
                    profiles.append(profile)
                except (ValueError, TypeError, KeyError, ZeroDivisionError) as exc:
                    errors.append(run["name"] + ": invalid profiling records: " + str(exc))
    reasons = []
    if mode != "throughput":
        reasons.append(str(mode) + " runs are separate from repeated profiling-disabled throughput evidence")
    if not gate_pass:
        reasons.append("planned runs and scientific/binding/card gates are not all complete and passed")
    if not profile_off:
        reasons.append("profiling is enabled or no completed runs are available")
    if not budget_valid or reference is None:
        reasons.append("fixed 128-CPU candidate layouts or baseline 128x1 reference are missing")
    if not repeats_sufficient:
        reasons.append("each compared 128-CPU layout requires at least three completed repeats for a throughput conclusion")
    reasons.extend(errors)
    return {"mode": mode, "status": "invalid" if errors else status, "host": report.get("host"),
            "card_sha256": report.get("card_sha256"), "case": report.get("case"),
            "binary_sha256": {key: value.get("sha256") for key, value in report.get("binaries", {}).items()},
            "allocation": report.get("allocation"), "runs_planned": len(expected), "runs_seen": len(runs),
            "runs_validated": len(validated), "planned_repeats": plan.get("repeats") if plan else None,
            "gates": {"planned_runs_complete": complete, "bindings_pass": bindings,
                      "same_rank_bitwise_pass": same_rank, "card_unchanged": report.get("card_unchanged"),
                      "cross_rank_tolerance_pass": cross_result, "cross_rank_comparisons": len(cross),
                      "cross_rank_is_gating": False},
            "throughput_eligible": eligible, "throughput_exclusion_reasons": reasons,
            "minimum_throughput_repeats_per_layout": MIN_THROUGHPUT_REPEATS,
            "fixed_128_cpu_table": fixed,
            "best_tested_exclusive_node_comparison": best_tested_node_comparison(
                table, report.get("allocation") or {}, eligible),
            "same_rank_controls_using_fewer_cpus": [row for row in table if row["step_cpus"] < CPU_BUDGET],
            "profiles": profiles, "report_error": report.get("error"), "analysis_errors": errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path, help="scientific/report.json paths")
    parser.add_argument("--output-dir", type=Path, required=True, help="new directory; existing directories are refused")
    args = parser.parse_args()
    paths = [path.resolve() for path in args.reports]
    if len(set(paths)) != len(paths):
        parser.error("duplicate input report paths are not allowed")
    snapshots = []
    try:
        for path in paths:
            payload = path.read_bytes()
            plan_path = path.with_name("plan.json")
            plan_payload = plan_path.read_bytes() if plan_path.is_file() else None
            snapshots.append((path, payload, plan_payload, json.loads(payload),
                              json.loads(plan_payload) if plan_payload is not None else None))
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=False)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    sources = output / "sources"
    sources.mkdir()
    summary = {"schema_version": 1, "reports_by_mode": {"calibration": [], "throughput": [], "profile": []},
               "unrecognized_reports": [], "notes": [
                   "Each report is summarized independently; no samples are pooled across reports, card hashes, binary hashes, or modes.",
                   "Source report/plan bytes are snapshotted with SHA256 so running-job snapshots remain reproducible.",
                   "Partial validated timings are descriptive only. Speedups require a completed, passed throughput report with profiling disabled.",
                   "Every compared 128-CPU layout requires at least three completed repeats before throughput speedups are reported.",
                   "Total time is APRL reported_wall_s including initialization; loop time includes final diagnostics/output closure.",
                   "Logical particle/grid rates normalize field steps and are not FLOPS or measured kernel operation counts.",
                   "Cross-rank differences are non-gating; baseline 64x1/32x1 controls do not use the fixed 128-CPU budget.",
                   "The separate best-tested node comparison selects baseline and candidate layouts by median APRL total time across all tested layouts, including baseline controls. It combines implementation and layout selection; loop/process speedups use those same selected layouts.",
                   "Profile mean_s/max_s are rank statistics for each reporting interval. Interval lengths may differ.",
                   "field_communication nests in field_total; particle_communication nests in particle_migration_total. No nested stages or stage maxima are added."]}
    for index, (path, payload, plan_payload, report, plan) in enumerate(snapshots, 1):
        report_snapshot = sources / f"report-{index:03d}.json"
        report_snapshot.write_bytes(payload)
        entry = analyze(report, plan)
        entry["source"] = {"path": str(path), "sha256": sha256(payload), "snapshot": str(report_snapshot.relative_to(output))}
        if plan_payload is not None:
            plan_snapshot = sources / f"plan-{index:03d}.json"
            plan_snapshot.write_bytes(plan_payload)
            entry["plan_source"] = {"path": str(path.with_name("plan.json")), "sha256": sha256(plan_payload),
                                    "snapshot": str(plan_snapshot.relative_to(output))}
        target = summary["reports_by_mode"].get(entry["mode"], summary["unrecognized_reports"])
        target.append(entry)
    write_json(output / "summary.json", summary)
    print("Summary:", output / "summary.json")
    for mode, entries in summary["reports_by_mode"].items():
        for entry in entries:
            print(mode, entry["status"], f"{entry['runs_validated']}/{entry['runs_planned']} validated runs;",
                  "throughput eligible:" + str(entry["throughput_eligible"]), entry["source"]["path"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
