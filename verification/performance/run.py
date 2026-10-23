#!/usr/bin/env python3
"""Controlled SASE-loader smoke/performance comparison; never changes physics between runs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time

ENV_KEYS = ("APRL_PROFILE_INTERVAL", "OMP_NUM_THREADS", "OMP_PLACES", "OMP_PROC_BIND", "OMP_DYNAMIC",
            "SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK",
            "OMPI_MCA_btl", "OMPI_MCA_pml", "OMPI_MCA_btl_tcp_if_include",
            "OMPI_MCA_btl_tcp_if_exclude", "UCX_TLS", "UCX_NET_DEVICES", "FI_PROVIDER")
PRODUCTS = (("ledger/energy-ledger.h5", "energy_ledger", "records", "committed_records"),
            ("detectors/particles.h5", "particle_plane", "records", "committed_records"),
            ("detectors/field.h5", "field_plane", "fields", "committed_samples"))


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def csv_positive(value):
    try:
        values = list(dict.fromkeys(int(part) for part in value.split(",")))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated positive integers") from exc
    if not values or min(values) < 1:
        raise argparse.ArgumentTypeError("expected positive integers")
    return values


def card(input_file, electrons):
    # Paths are deliberately relative to each run's cwd; all cards are byte-identical.
    return f'''# Fixed small SASE-loading computational check, not a production gain model.
units: {{length: micrometer, time: femtosecond}}
runtime: {{mode: throughput, resource_monitor: {{enabled: true, calibration_steps: 1, progress_interval_steps: 10}}}}
output: {{overwrite: false, manifest: run-manifest.yaml}}
mesh:
  field_solver: cowan-z
  current_filter: none
  cells: [24, 24, 128]
  cell_size: [0.25, 0.25, 0.1]
  center: [0, 0, 0]
  duration: 30
  boost_gamma: 2
  particle_steps_per_undulator_period: 16
boundary: {{type: pec}}
initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-9
  maximum_iterations: 2000
  maximum_mean_transverse_beta: 0.01
beam:
  reference: {{initial_center_z: -1.0}}
  input:
    type: hdf5
    file: {json.dumps(str(input_file))}
    electrons: {electrons}
    position_offset: [0, 0, 0]
sources:
  incident_waves: []
  magnetic_elements:
    - {{type: planar-undulator, strength_parameter: 0.1, period: 0.5, periods: 4, entrance_z: 0, gaussian_fringe: false}}
detectors:
  enabled: true
  directory: detectors
  field_planes:
    - {{name: field, z: 2.5, rhythm: 0.1, buffer_samples: 8, compression: 0}}
  particle_planes:
    - {{name: particles, z: 3.0, buffer_records: 4096, compression: 0}}
stop: {{mode: after-last-element}}
trajectory: {{enabled: false}}
energy_ledger:
  enabled: true
  directory: ledger
  filename: energy-ledger.h5
  sample_interval_steps: 1
  buffer_records: 256
  compression: 0
'''


def launch_command(args, executable, ranks, threads, config):
    if args.launcher == "none":
        if ranks != 1:
            raise ValueError("--launcher none requires --ranks 1")
        prefix = []
    elif args.launcher == "mpiexec":
        prefix = [args.mpiexec, "-n", str(ranks), "--bind-to", args.binding]
        if args.binding == "core":
            prefix += ["--map-by", "slot:PE=" + str(threads)]
    else:
        prefix = ["srun", "--mpi=" + args.srun_mpi, "--ntasks=" + str(ranks),
                  "--cpus-per-task=" + str(threads), "--cpu-bind=cores",
                  "--kill-on-bad-exit=1"]
        allocated_nodes = os.environ.get("SLURM_JOB_NUM_NODES") or os.environ.get("SLURM_NNODES")
        if allocated_nodes:
            nodes = min(ranks, int(allocated_nodes))
            prefix += ["--nodes=" + str(nodes), "--ntasks-per-node=" + str((ranks + nodes - 1) // nodes)]
    return prefix + args.launcher_arg + [str(executable), str(config)]


def run_command(command, cwd, env, log):
    start = time.monotonic()
    with log.open("w") as stream:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    return result.returncode, time.monotonic() - start


def numeric_fields(array, prefix):
    import numpy as np
    if array.dtype.names:
        for name in array.dtype.names:
            yield from numeric_fields(array[name], prefix + "." + name)
    elif np.issubdtype(array.dtype, np.number):
        yield prefix, array


def inspect_run(directory, expected_particles):
    import h5py
    import numpy as np
    text = (directory / "run.log").read_text(errors="replace")
    if "Configured stop reached:" not in text:
        raise ValueError("missing normal physical-stop marker")
    if "duration was exhausted" in text or "stopped after a complete field step" in text:
        raise ValueError("run stopped before its physical endpoint")
    arrays = {}
    for filename, group_name, data_name, commit_name in PRODUCTS:
        with h5py.File(directory / filename, "r") as file:
            group = file[group_name]
            if int(group["complete"][()]) != 1:
                raise ValueError(filename + ": completion marker is false")
            data = group[data_name][:]
            if not len(data) or int(group[commit_name][()]) != len(data):
                raise ValueError(filename + ": empty/uncommitted data")
            if group_name == "particle_plane":
                if len(data) != expected_particles or len(np.unique(data["source_id"])) != expected_particles:
                    raise ValueError("particle endpoint count/unique source IDs do not match input")
                data = data[np.argsort(data["source_id"], kind="stable")]
            arrays.update(numeric_fields(data, group_name))
            if group_name == "field_plane":
                times = group["time_s"][:]
                if len(times) != len(data) or np.any(np.diff(times) <= 0):
                    raise ValueError(filename + ": invalid field time axis")
                arrays[group_name + ".time_s"] = times
    for name, array in arrays.items():
        if not np.all(np.isfinite(array)):
            raise ValueError(name + ": non-finite output")
    steps = arrays["energy_ledger.step"]
    if steps[0] != 0 or steps[-1] <= 0 or np.any(np.diff(steps.astype(np.int64)) < 0):
        raise ValueError("invalid ledger step sequence")
    if np.any(arrays["energy_ledger.removed_macroparticles"] != 0):
        raise ValueError("unexpected particle losses")
    wall = re.search(r"Elapsed wall time \[s\]:\s*([0-9.eE+-]+)", text)
    if not wall:
        raise ValueError("missing completed wall-time measurement")
    thread_match = re.search(r"Compute threads per MPI rank:\s*(\d+)", text)
    performance = [dict(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)=([^\s,]+)", line))
                   for line in text.splitlines() if line.startswith("[performance]")]
    progress = [dict(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)=([^\s,]+)", line))
                for line in text.splitlines() if line.startswith("[resource] phase=final")]
    summary = {"steps": int(steps[-1]), "reported_wall_s": float(wall.group(1)),
               "particle_records": expected_particles, "reported_threads": int(thread_match.group(1)) if thread_match else None, "performance": performance,
               "final_resource": progress,
               "thread_startup": [line for line in text.splitlines() if "[parallel]" in line or "threads=" in line or "Compute threads per MPI rank:" in line]}
    if progress:
        summary["loop_wall_s"] = float(progress[-1]["elapsed_s"])
        summary["loop_mean_step_s"] = float(progress[-1]["measured_step_s"])
    summary["wall_per_step_including_initialization_s"] = summary["reported_wall_s"] / summary["steps"]
    return summary, arrays


def compare(reference, candidate, rtol, scaled_atol, strict_bitwise=False):
    import numpy as np
    if reference.keys() != candidate.keys():
        raise ValueError("output schemas differ")
    report = {}
    for name, expected in reference.items():
        actual = candidate[name]
        if expected.shape != actual.shape or expected.dtype != actual.dtype:
            raise ValueError(name + ": shape/type differs")
        if expected.dtype.kind in "iub":
            passed = np.array_equal(expected, actual)
            report[name] = {"pass": bool(passed), "exact": True}
        else:
            scale = float(np.max(np.abs(expected)))
            difference = float(np.max(np.abs(actual - expected)))
            # Relative to this physical field's scale, not a unit-blind SI absolute tolerance.
            absolute = scaled_atol * scale
            passed = bool(np.all(np.abs(actual - expected) <= absolute + rtol * np.abs(expected)))
            bitwise = np.ascontiguousarray(expected).tobytes() == np.ascontiguousarray(actual).tobytes()
            report[name] = {"pass": bool(bitwise if strict_bitwise else passed), "bitwise_equal": bitwise, "max_abs_difference": difference,
                            "reference_scale": scale, "scaled_max_difference": difference / scale if scale else difference}
    return {"pass": all(row["pass"] for row in report.values()), "fields": report}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ranks", type=csv_positive, default=[1, 2])
    parser.add_argument("--threads", type=csv_positive, default=[1, 2])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--beamlets", type=int, default=16, help="fixed across every run; changes the input record count")
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--launcher", choices=("none", "mpiexec", "srun"), default="mpiexec")
    parser.add_argument("--binding", choices=("core", "none"), default="core", help="Open MPI binding; none is for an unbound local check")
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--srun-mpi", default=os.environ.get("SRUN_MPI", "pmix_v3"))
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--scaled-atol", type=float, default=1e-12)
    parser.add_argument("--require-cross-rank", action="store_true", help="also fail on baseline-inherited cross-rank differences")
    parser.add_argument("--profile-interval", type=int, default=16, help="candidate stage timing interval; 0 disables profiling")
    parser.add_argument("--plan", action="store_true", help="write plan/cards only; do not generate particles or run APRL")
    args = parser.parse_args()
    if min(args.repeats, args.beamlets) < 1 or args.rtol < 0 or args.scaled_atol < 0 or args.profile_interval < 0:
        parser.error("repeats/beamlets must be positive; tolerances must be nonnegative")
    if max(args.ranks) > 64:
        parser.error("small 128-cell z mesh supports at most 64 ranks; use a separate production card for larger scaling")
    if args.launcher == "none" and args.ranks != [1]:
        parser.error("--launcher none requires --ranks 1")
    if args.launcher == "srun" and not args.plan:
        allocated_tasks = int(os.environ.get("SLURM_NTASKS", "0"))
        allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
        if not allocated_tasks or max(args.ranks) > allocated_tasks or max(args.threads) > allocated_cpus:
            parser.error("srun requires an existing allocation fitting max(ranks) and max(threads) CPUs per task")
    args.baseline, args.candidate, args.generator = (p.resolve() for p in (args.baseline, args.candidate, args.generator))
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    input_file = root / "particles.h5"
    generator = [str(args.generator), "--output", str(input_file), "--electrons", "100000",
                 "--length-m", "4e-7", "--resonant-wavelength-m", "2e-8", "--gamma", "50",
                 "--sigma-gamma", "0.002", "--sigma-x-m", "1e-7", "--sigma-y-m", "1e-7",
                 "--emittance-n-x-m", "1e-10", "--emittance-n-y-m", "1e-10",
                 "--harmonics", "3", "--beamlets-per-wavelength", str(args.beamlets), "--seed", str(args.seed)]
    runs = []
    for repeat in range(args.repeats):
        for ranks in args.ranks:
            variants = [("baseline", args.baseline, 1)] + [("candidate", args.candidate, n) for n in args.threads]
            # Reverse the order every repetition to reduce systematic warm-cache/order bias.
            if repeat % 2:
                variants.reverse()
            for label, executable, threads in variants:
                name = f"{label}-r{ranks}-t{threads}-rep{repeat + 1}"
                directory = root / name
                directory.mkdir()
                config = directory / "run.yaml"
                config.write_text(card(input_file, 100000), encoding="utf-8")
                runs.append({"name": name, "label": label, "ranks": ranks, "threads": threads,
                             "repeat": repeat + 1, "directory": str(directory),
                             "command": launch_command(args, executable, ranks, threads, config)})
    plan = {"format": "aprl-performance-v1", "generator": generator, "runs": runs,
            "rtol": args.rtol, "scaled_atol": args.scaled_atol,
            "card_sha256": digest(root / runs[0]["name"] / "run.yaml"),
            "profile_interval": args.profile_interval, "require_cross_rank": args.require_cross_rank,
            "description": "Fixed HDF5-v2 SASE-style loading, short undulator, coupled self-fields; not a production gain/convergence claim"}
    (root / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print("Plan:", root / "plan.json", flush=True)
    if args.plan:
        return 0
    import h5py
    import numpy as np  # noqa: F401
    for executable in (args.baseline, args.candidate, args.generator):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError("not executable: " + str(executable))
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE", OMP_PLACES=env.get("OMP_PLACES", "cores"),
               OMP_PROC_BIND=env.get("OMP_PROC_BIND", "false" if args.binding == "none" and args.launcher == "mpiexec" else "close"))
    code, _ = run_command(generator, root, env, root / "generator.log")
    if code:
        raise ValueError("generator failed; see generator.log")
    with h5py.File(input_file, "r") as file:
        if int(file["particles"].attrs["complete"]) != 1:
            raise ValueError("generated input is incomplete")
        expected_particles = len(file["particles/records"])
    report = {"host": platform.node(), "platform": platform.platform(),
              "timing_notes": ["reported_wall_s includes initialization", "loop_wall_s/loop_mean_step_s come from the final resource report and include final diagnostics/output closure", "communication stages are nested in total stages and include MPI wait time"],
              "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
              "input_sha256": digest(input_file), "input_records": expected_particles,
              "binary_sha256": {"baseline": digest(args.baseline), "candidate": digest(args.candidate), "generator": digest(args.generator)},
              "environment": {key: env.get(key) for key in ENV_KEYS}, "runs": [], "comparisons": []}
    loaded = {}
    for run in runs:
        directory = Path(run["directory"])
        env["OMP_NUM_THREADS"] = str(run["threads"])
        env["APRL_PROFILE_INTERVAL"] = str(args.profile_interval if run["label"] == "candidate" else 0)
        print("Running", run["name"], flush=True)
        (directory / "launch.json").write_text(json.dumps({"command": run["command"],
            "environment": {key: env.get(key) for key in ENV_KEYS}}, indent=2) + "\n")
        code, elapsed = run_command(run["command"], directory, env, directory / "run.log")
        row = {**run, "returncode": code, "process_wall_s": elapsed}
        try:
            if code:
                raise ValueError("nonzero process exit")
            summary, arrays = inspect_run(directory, expected_particles)
            if run["label"] == "candidate" and run["threads"] > 1 and summary["reported_threads"] != run["threads"]:
                raise ValueError("candidate did not confirm the requested effective thread count")
            row.update(summary, valid=True)
            loaded[run["name"]] = arrays
        except (ValueError, KeyError, OSError) as exc:
            row.update(valid=False, error=str(exc))
        report["runs"].append(row)
        (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    references = {run["ranks"]: run["name"] for run in reversed(runs) if run["label"] == "baseline"}
    canonical = references[args.ranks[0]]
    for run in runs:
        name = run["name"]
        reference_name = references[run["ranks"]]
        if name not in loaded or reference_name not in loaded:
            continue
        try:
            result = compare(loaded[reference_name], loaded[name], args.rtol, args.scaled_atol, strict_bitwise=True)
        except ValueError as exc:
            result = {"pass": False, "error": str(exc)}
        report["comparisons"].append({"reference": reference_name, "candidate": name,
                                       "mode": "same-rank-bitwise", **result})
        if run["ranks"] != args.ranks[0] and canonical in loaded:
            try:
                cross = compare(loaded[canonical], loaded[name], args.rtol, args.scaled_atol)
            except ValueError as exc:
                cross = {"pass": False, "error": str(exc)}
            report["comparisons"].append({"reference": canonical, "candidate": name,
                                           "mode": "cross-rank-tolerance", **cross})
    valid = all(row["valid"] for row in report["runs"])
    expected_comparisons = len(runs) + sum(run["ranks"] != args.ranks[0] for run in runs)
    report["implementation_equivalence_pass"] = valid and len(report["comparisons"]) == expected_comparisons and all(c["pass"] for c in report["comparisons"] if c["mode"] == "same-rank-bitwise")
    report["cross_rank_consistency_pass"] = valid and all(c["pass"] for c in report["comparisons"] if c["mode"] == "cross-rank-tolerance")
    report["pass"] = report["implementation_equivalence_pass"] and (not args.require_cross_rank or report["cross_rank_consistency_pass"])
    groups = {}
    for row in report["runs"]:
        if row["valid"]:
            key = f'{row["label"]}-r{row["ranks"]}-t{row["threads"]}'
            groups.setdefault(key, []).append(row["reported_wall_s"])
    report["timing_summary"] = {key: {"samples_s": values, "median_s": statistics.median(values),
                                              "min_s": min(values), "max_s": max(values)} for key, values in groups.items()}
    report["input_unchanged"] = digest(input_file) == report["input_sha256"]
    report["pass"] = report["pass"] and report["input_unchanged"]
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"pass": report["pass"], "implementation_equivalence_pass": report["implementation_equivalence_pass"], "cross_rank_consistency_pass": report["cross_rank_consistency_pass"], "timing": report["timing_summary"]}, indent=2))
    print("Report:", root / "report.json")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, ImportError) as error:
        raise SystemExit("ERROR: " + str(error))
