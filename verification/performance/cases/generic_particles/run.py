#!/usr/bin/env python3
"""Bounded, fixed-128-CPU calibration of immutable generated-particle cases."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time

import yaml


CPU_BUDGET = 128
RTOL = 1e-9
SCALED_ATOL = 1e-12
ENV_KEYS = (
    "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PLACES", "OMP_PROC_BIND", "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "APRL_PROFILE_INTERVAL",
    "SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_JOB_NUM_NODES", "SLURM_NTASKS",
    "SLURM_CPUS_PER_TASK", "SLURM_JOB_CPUS_PER_NODE", "SLURM_CPUS_ON_NODE", "SRUN_MPI",
    "OMPI_MCA_btl", "OMPI_MCA_pml", "UCX_TLS", "UCX_NET_DEVICES", "FI_PROVIDER",
)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def parse_layouts(value):
    layouts = []
    for item in value.split(","):
        match = re.fullmatch(r"([1-9][0-9]*)x([1-9][0-9]*)", item.strip())
        if not match:
            raise argparse.ArgumentTypeError("layouts must be comma-separated RxT pairs, such as 128x1,64x2,32x4")
        layout = tuple(map(int, match.groups()))
        if math.prod(layout) != CPU_BUDGET:
            raise argparse.ArgumentTypeError("each candidate layout must use exactly 128 CPUs")
        if layout not in layouts:
            layouts.append(layout)
    return layouts


def positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(name + " must be a positive integer")
    return value


def read_case(path):
    payload = path.read_bytes()
    card = yaml.safe_load(payload)
    beam = card["beam"]["input"]
    if beam["type"] != "generated-gaussian" or "file" in beam:
        raise ValueError("case must use built-in generated-gaussian particles without an input file")
    particles = positive_int(beam["macroparticles"], "beam.input.macroparticles")
    seed = beam["random_seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 0xffffffff:
        raise ValueError("case must specify a fixed unsigned 32-bit random_seed")
    mesh = card["mesh"]
    if mesh.get("self_field_feedback", True) is not True:
        raise ValueError("self_field_feedback must remain enabled")
    cells = mesh["cells"]
    if not isinstance(cells, list) or len(cells) != 3:
        raise ValueError("case must specify three fixed mesh cell counts")
    cells = [positive_int(value, "mesh.cells") for value in cells]
    if card.get("output", {}).get("overwrite", False) is not False:
        raise ValueError("case must disable output overwrite")
    if card["output"]["manifest"] != "run-manifest.yaml":
        raise ValueError("case manifest must be relative run-manifest.yaml")
    detectors, ledger = card["detectors"], card["energy_ledger"]
    if detectors.get("enabled") is not True or detectors["directory"] != "detectors":
        raise ValueError("case must enable relative detectors output")
    if [plane["name"] for plane in detectors["field_planes"]] != ["field"] or \
            [plane["name"] for plane in detectors["particle_planes"]] != ["particles"]:
        raise ValueError("case requires one field and one particles diagnostic plane")
    if ledger.get("enabled") is not True or ledger["directory"] != "ledger" or \
            ledger["filename"] != "energy-ledger.h5":
        raise ValueError("case must enable relative ledger/energy-ledger.h5")
    if card.get("trajectory", {}).get("enabled", False):
        raise ValueError("trajectory output is outside this benchmark's frozen diagnostic contract")
    monitor = card["runtime"]["resource_monitor"]
    if monitor.get("enabled") is not True:
        raise ValueError("resource monitoring must be enabled to measure the simulation loop")
    if card["stop"]["mode"] != "reference-center-z":
        raise ValueError("case must have a fixed reference-center-z physical endpoint")
    return payload, {"expected_macroparticles": particles, "represented_electrons": beam["electrons"],
                     "random_seed": seed, "mesh_cells": cells, "grid_cells": math.prod(cells),
                     "self_field_feedback": True, "physical_card": card}


def allocation():
    cpus = os.environ.get("SLURM_JOB_CPUS_PER_NODE", os.environ.get("SLURM_CPUS_ON_NODE", "0"))
    match = re.fullmatch(r"([0-9]+)(?:\(x1\))?", cpus)
    nodes = os.environ.get("SLURM_JOB_NUM_NODES", os.environ.get("SLURM_NNODES", "0"))
    if not os.environ.get("SLURM_JOB_ID") or nodes != "1" or not match or int(match.group(1)) < CPU_BUDGET:
        raise ValueError("execution requires an existing single-node Slurm allocation of at least 128 CPUs")
    return {"nodes": 1, "allocated_cpus": int(match.group(1)), "candidate_cpu_budget": CPU_BUDGET,
            "allocation_ntasks": os.environ.get("SLURM_NTASKS"),
            "allocation_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
            "exclusive_requested_by_batch_script": True}


def srun(plugin, ranks, threads):
    if ranks * threads > CPU_BUDGET:
        raise ValueError("step exceeds fixed CPU budget")
    return ["timeout", "--kill-after=30", "600", "srun", "--exact", "--mpi=" + plugin,
            "--nodes=1", "--ntasks=" + str(ranks), "--ntasks-per-node=" + str(ranks),
            "--cpus-per-task=" + str(threads), "--cpu-bind=verbose,cores", "--kill-on-bad-exit=1"]


def execute(command, directory, threads, interval, logfile):
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", OMP_PLACES="cores",
                       OMP_PROC_BIND="close", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       NUMEXPR_NUM_THREADS="1", APRL_PROFILE_INTERVAL=str(interval))
    write_json(logfile.with_suffix(".launch.json"), {"command": command, "cwd": str(directory),
               "environment": {key: environment.get(key) for key in ENV_KEYS}})
    start = time.monotonic()
    with logfile.open("w") as stream:
        result = subprocess.run(command, cwd=directory, env=environment, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    return result.returncode, time.monotonic() - start


def load_inspector(workspace):
    path = workspace / "source/aprl/verification/performance/run.py"
    spec = importlib.util.spec_from_file_location("aprl_fixed_physics_inspector", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def inspect(inspector, run, case):
    import h5py
    import numpy as np
    directory = Path(run["directory"])
    count = case["expected_macroparticles"]
    summary, arrays = inspector.inspect_run(directory, count)
    ids = arrays["particle_plane.source_id"]
    if not np.array_equal(ids, np.arange(1, count + 1, dtype=ids.dtype)):
        raise ValueError("particle source IDs must be exactly 1 through the generated macroparticle count")
    if np.any(arrays["energy_ledger.active_macroparticles"] != count):
        raise ValueError("ledger active macroparticle count changed")
    if np.any(arrays["energy_ledger.closure_valid"] != 1):
        raise ValueError("energy ledger closure is not valid with self-field feedback")
    if np.any(arrays["energy_ledger.boundary_injection_entries"] != 0):
        raise ValueError("generic generated beam unexpectedly used boundary injection")
    if run["label"] == "candidate" and summary["reported_threads"] != run["threads"]:
        raise ValueError("candidate did not confirm requested effective thread count")
    with h5py.File(directory / "detectors/particles.h5", "r") as file:
        if int(file["particle_plane"].attrs["mpi_size"]) != run["ranks"]:
            raise ValueError("particle output MPI size differs from requested rank count")
    loop_wall = summary.get("loop_wall_s", 0)
    if not math.isfinite(loop_wall) or loop_wall <= 0 or summary["reported_wall_s"] <= 0:
        raise ValueError("missing or invalid completed-loop timing")
    steps = summary["steps"]
    summary.update(loop_seconds_per_step=loop_wall / steps,
                   logical_particle_updates_per_s=count * steps / loop_wall,
                   logical_grid_cell_updates_per_s=case["grid_cells"] * steps / loop_wall,
                   particle_source_ids_complete=True, macroparticle_conservation_pass=True,
                   energy_ledger_closure_valid=True)
    return summary, arrays


def compare(inspector, reference, actual, strict):
    try:
        return inspector.compare(reference, actual, RTOL, SCALED_ATOL, strict_bitwise=strict)
    except (ValueError, KeyError) as exc:
        return {"pass": False, "error": str(exc)}


def summarize(report):
    groups = {}
    metrics = ("process_wall_s", "reported_wall_s", "loop_wall_s", "loop_seconds_per_step",
               "logical_particle_updates_per_s", "logical_grid_cell_updates_per_s")
    for row in report["runs"]:
        if row.get("valid"):
            groups.setdefault((row["label"], row["ranks"], row["threads"]), []).append(row)
    report["timing_summary"] = [
        {"label": label, "ranks": ranks, "threads": threads, "step_cpus": ranks * threads,
         "repeats_completed": len(rows), "profile_interval": rows[0]["profile_interval"],
         "statistics": {metric: {"median": statistics.median(row[metric] for row in rows),
                                  "min": min(row[metric] for row in rows),
                                  "max": max(row[metric] for row in rows)} for metric in metrics}}
        for (label, ranks, threads), rows in groups.items()]
    if report["mode"] == "throughput" and report.get("pass"):
        medians = {(row["label"], row["ranks"], row["threads"]): row["statistics"]
                   for row in report["timing_summary"]}
        reference = medians[("baseline", CPU_BUDGET, 1)]
        report["throughput_comparisons"] = [
            {"candidate_ranks": ranks, "candidate_threads": threads, "fixed_cpu_budget": CPU_BUDGET,
             "fixed_budget_loop_speedup": reference["loop_wall_s"]["median"] / stats["loop_wall_s"]["median"],
             "fixed_budget_total_speedup": reference["reported_wall_s"]["median"] / stats["reported_wall_s"]["median"],
             "same_rank_loop_speedup": medians[("baseline", ranks, 1)]["loop_wall_s"]["median"] / stats["loop_wall_s"]["median"],
             "same_rank_baseline_cpus": ranks}
            for (label, ranks, threads), stats in medians.items() if label == "candidate"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mode", choices=("calibration", "throughput", "profile"), default="calibration")
    parser.add_argument("--layouts", type=parse_layouts, default=parse_layouts("128x1,64x2,32x4"))
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--srun-mpi", default=os.environ.get("SRUN_MPI", "pmix_v3"))
    parser.add_argument("--plan", action="store_true", help="freeze cards and commands only; no simulation or binding probe")
    args = parser.parse_args()
    if not args.workspace.is_absolute() or not args.case.is_absolute():
        parser.error("--workspace and --case must be absolute paths")
    repeats = args.repeats if args.repeats is not None else (3 if args.mode == "throughput" else 1)
    if repeats < 1:
        parser.error("--repeats must be positive")
    workspace, root, card_path = args.workspace.resolve(), args.output_dir.resolve(), args.case.resolve()
    try:
        card_bytes, case = read_case(card_path)
        if max(CPU_BUDGET, *(ranks for ranks, _ in args.layouts)) > case["mesh_cells"][2] // 2:
            raise ValueError("longitudinal mesh requires at least two cells per MPI rank")
        resource = None if args.plan else allocation()
    except (ValueError, KeyError, TypeError, OSError, yaml.YAMLError) as exc:
        parser.error(str(exc))
    root.mkdir(parents=True, exist_ok=False)
    card_hash = hashlib.sha256(card_bytes).hexdigest()
    (root / "case.yaml").write_bytes(card_bytes)
    binaries = {label: workspace / "build" / label / "aprl" for label in ("baseline", "candidate")}
    probe_executable = workspace / "build/cluster-tools/binding_probe"
    controls = [(CPU_BUDGET, 1)] + [(ranks, 1) for ranks, _ in args.layouts if ranks != CPU_BUDGET]
    controls = list(dict.fromkeys(controls))
    configurations = [("baseline", *layout) for layout in controls]
    configurations += [("candidate", *layout) for layout in args.layouts]
    # Pair baseline/candidate in the first repetition; reverse the complete order thereafter.
    configurations.sort(key=lambda row: (-row[1], row[0] == "candidate"))
    probes = [{"ranks": ranks, "threads": threads, "step_cpus": ranks * threads,
               "command": srun(args.srun_mpi, ranks, threads) + [str(probe_executable), str(ranks), str(threads)]}
              for ranks, threads in dict.fromkeys(args.layouts + controls)]
    runs = []
    for repeat in range(1, repeats + 1):
        ordered = configurations if repeat % 2 else list(reversed(configurations))
        for label, ranks, threads in ordered:
            name = f"rep{repeat:02d}-{label}-r{ranks}-t{threads}"
            directory = root / name
            directory.mkdir()
            (directory / "run.yaml").write_bytes(card_bytes)
            interval = 64 if args.mode == "profile" and label == "candidate" else 0
            runs.append({"name": name, "label": label, "ranks": ranks, "threads": threads,
                         "step_cpus": ranks * threads, "repeat": repeat, "profile_interval": interval,
                         "directory": str(directory), "card_sha256": card_hash,
                         "command": srun(args.srun_mpi, ranks, threads) +
                         [str(binaries[label]), str(directory / "run.yaml")]})
    plan = {"mode": args.mode, "repeats": repeats, "card": str(card_path), "card_sha256": card_hash,
            "case": case, "cpu_budget": CPU_BUDGET, "step_timeout_s": 600, "kill_after_s": 30,
            "allocation": resource, "probes": probes, "runs": runs,
            "external_particle_input": None, "generator": None,
            "comparison_policy": {"same_rank": "strict-bitwise-gating", "cross_rank": "reported-non-gating",
                                  "rtol": RTOL, "scaled_atol": SCALED_ATOL}}
    write_json(root / "plan.json", plan)
    if args.plan:
        print("Frozen plan:", root / "plan.json")
        return 0
    report = {"pass": False, "mode": args.mode, "host": platform.node(), "allocation": resource,
              "case": case, "card_sha256": card_hash, "probes": [], "runs": [], "comparisons": [],
              "cross_rank_comparisons": [], "environment": {key: os.environ.get(key) for key in ENV_KEYS},
              "timing_notes": [
                  "process_wall_s includes launcher startup; reported_wall_s includes APRL initialization.",
                  "loop_wall_s includes final diagnostics and output closure; seconds per step uses completed ledger steps.",
                  "Particle and grid cell rates normalize logical completed field steps; they are not FLOPS or measured kernel operations.",
                  "Single-thread baseline controls below 128 ranks use fewer CPUs; only baseline 128x1 is the fixed-budget performance reference.",
                  "Profile mode uses candidate interval 64 and baseline interval 0; its timings are not throughput evidence.",
                  "field_communication is nested in field_total, and particle_communication in migration_total; stages include MPI waits and maxima must not be added.",
                  "Cross-rank differences are preserved using rtol=1e-9 and scaled_atol=1e-12 and never gate same-rank equivalence."]}
    report_path = root / "report.json"
    write_json(report_path, report)
    try:
        inspector, inspector_path = load_inspector(workspace)
        report["inspector"] = {"path": str(inspector_path), "sha256": digest(inspector_path)}
        report["binaries"] = {label: {"path": str(path), "sha256": digest(path)} for label, path in binaries.items()}
        report["binding_probe"] = {"path": str(probe_executable), "sha256": digest(probe_executable)}
        prepared_record = workspace / "tests/cluster_acceptance/build-record.json"
        report["build_record"] = {"path": str(prepared_record), "sha256": digest(prepared_record),
                                  "record": json.loads(prepared_record.read_text())}
        for probe in probes:
            row = {**probe, "status": "running"}
            report["probes"].append(row)
            write_json(report_path, report)
            code, elapsed = execute(probe["command"], root, probe["threads"], 0,
                                    root / f"binding-r{probe['ranks']}-t{probe['threads']}.log")
            row.update(returncode=code, process_wall_s=elapsed, status="passed" if code == 0 else "failed")
            write_json(report_path, report)
            if code:
                raise ValueError("binding probe failed for " + str((probe["ranks"], probe["threads"])))
        references, reference_names = {}, {}
        for run in runs:
            directory = Path(run["directory"])
            row = {**run, "valid": False, "status": "running"}
            report["runs"].append(row)
            write_json(report_path, report)
            if digest(directory / "run.yaml") != card_hash:
                raise ValueError(run["name"] + ": physical card changed before launch")
            print("Running", run["name"], flush=True)
            code, elapsed = execute(run["command"], directory, run["threads"], run["profile_interval"], directory / "run.log")
            row.update(returncode=code, process_wall_s=elapsed, status="inspecting" if code == 0 else "failed")
            write_json(report_path, report)
            if code:
                raise ValueError(run["name"] + ": nonzero process exit (124 means step timeout)")
            if digest(directory / "run.yaml") != card_hash:
                raise ValueError(run["name"] + ": physical card changed during execution")
            summary, arrays = inspect(inspector, run, case)
            row.update(summary, valid=True, status="passed")
            ranks = run["ranks"]
            if ranks not in references:
                if run["label"] != "baseline":
                    raise ValueError("first same-rank reference must be baseline")
                references[ranks], reference_names[ranks] = arrays, run["name"]
            comparison = {"reference": reference_names[ranks], "candidate": run["name"],
                          "mode": "same-rank-bitwise", **compare(inspector, references[ranks], arrays, True)}
            report["comparisons"].append(comparison)
            if ranks != CPU_BUDGET:
                report["cross_rank_comparisons"].append({
                    "reference": reference_names[CPU_BUDGET], "candidate": run["name"],
                    "mode": "cross-rank-tolerance-non-gating", "rtol": RTOL, "scaled_atol": SCALED_ATOL,
                    **compare(inspector, references[CPU_BUDGET], arrays, False)})
            del arrays
            write_json(report_path, report)
            if not comparison["pass"]:
                raise ValueError(run["name"] + ": same-rank bitwise equivalence failed")
        report["card_unchanged"] = digest(card_path) == card_hash and digest(root / "case.yaml") == card_hash
        report["implementation_equivalence_pass"] = all(row["pass"] for row in report["comparisons"])
        report["cross_rank_tolerance_pass"] = all(row["pass"] for row in report["cross_rank_comparisons"])
        report["pass"] = report["card_unchanged"] and report["implementation_equivalence_pass"]
        report["status"] = "passed" if report["pass"] else "failed"
    except (ValueError, OSError, KeyError, TypeError, ImportError, RuntimeError) as exc:
        report.update(error=str(exc), status="failed")
    summarize(report)
    write_json(report_path, report)
    print("Generic particle case:", "PASS" if report["pass"] else "FAIL", report_path)
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
