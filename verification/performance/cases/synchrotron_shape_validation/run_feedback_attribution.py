"""Run one shape's same-binary full/no-feedback synchrotron controls."""
import json
import os
from pathlib import Path
import re
import sys

import h5py
import numpy as np

from launcher import launch, binding
from prepare_paired import require, save, sha


def verify_case(directory, plan):
    no_feedback = plan["mode"] == "no-feedback"
    log = (directory / "aprl.log").read_text()
    markers = ["Configured stop reached:", "phase=solve-complete",
               "phase=mpi-finalize-complete"]
    if no_feedback:
        markers += ["self_field_feedback is disabled", "[one-way-source]"]
    for marker in markers:
        require(marker in log, "Missing completion/control marker " + marker)
    require(not re.search(
        r"Simulation failure|Segmentation fault|OpenFabrics device|DUE TO TIME LIMIT",
        log), "Simulation error marker")
    require(("self_field_feedback is disabled" in log) == no_feedback,
            "Feedback-mode log marker differs")
    require(("[one-way-source]" in log) == no_feedback,
            "One-way diagnostic marker differs")
    with h5py.File(directory / "detectors/screen.h5") as handle:
        group = handle["field_plane"]
        require(int(group["complete"][()]) == 1, "Screen field is incomplete")
        field_samples = int(group["committed_samples"][()])
        require(field_samples == len(group["time_s"]) == len(group["fields"]) and
                field_samples > 32, "Screen field commit count differs")
    for name in ["entry", "bend-exit", "exit"]:
        with h5py.File(directory / "detectors" / (name + ".h5")) as handle:
            group = handle["particle_plane"]
            require(int(group["complete"][()]) == 1,
                    "Particle plane incomplete: " + name)
            count = int(group["committed_records"][()])
            require(count == len(group["records"]) == plan["expected_particles"],
                    "Particle plane cohort differs: " + name)
    with h5py.File(directory / "ledger/energy-ledger.h5") as handle:
        group = handle["energy_ledger"]
        records = group["records"][:]
        require(int(group["complete"][()]) == 1 and
                int(group["committed_records"][()]) == len(records),
                "Energy ledger incomplete")
        require(int(records["step"][-1]) == plan["expected_steps"],
                "Energy ledger final step differs")
        if no_feedback:
            require(np.all(records["closure_valid"] == 0),
                    "No-feedback ledger unexpectedly claims closure")
            require(np.all(records["particle_grid_work_J"] == 0),
                    "Grid work is nonzero although feedback is disabled")
    paths = sorted((directory / "ledger").glob("discrete-field-rank-*.csv"))
    require([path.name for path in paths] ==
            [f"discrete-field-rank-{rank:05d}.csv" for rank in range(16)],
            "Missing strict-ledger rank")
    for path in paths:
        require(path.with_suffix(".complete").read_text().split() ==
                ["1", str(plan["expected_steps"])], "Unclosed strict ledger")
    one_way = sorted((directory / "reference").glob("one-way-source-rank-*.csv"))
    if no_feedback:
        require(len(one_way) == 16, "Missing one-way passive-work rank")
        for path in one_way:
            require(path.with_suffix(".complete").read_text().split() ==
                    ["1", str(plan["expected_steps"])],
                    "Unclosed one-way passive-work rank")
    else:
        require(not one_way, "Full-feedback case wrote one-way diagnostics")
    return dict(screen_samples=field_samples,
                field_steps=plan["expected_steps"],
                particle_records_per_plane=plan["expected_particles"],
                strict_ledger_ranks=16,
                particle_grid_work_exactly_zero=(no_feedback),
                one_way_work_ranks=(16 if no_feedback else 0))


def main(root):
    root = Path(root).resolve()
    campaign = json.loads((root / "campaign.json").read_text())
    frozen = json.loads((root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed")
    require(os.environ.get("SLURM_NTASKS") == "16", "Expected 16 MPI ranks")
    task = int(os.environ["SLURM_ARRAY_TASK_ID"])
    shapes = campaign["shapes"]
    require(0 <= task < len(shapes), "Invalid array task")
    shape = shapes[task]
    require(sha(campaign["binary"]) == campaign["common"]["binary_sha256"],
            "Binary changed")
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("APRL_")}
    env.update(OMP_NUM_THREADS="8", OMP_DYNAMIC="FALSE", OMP_PLACES="cores",
               OMP_PROC_BIND="close", OMPI_MCA_btl="self,vader,tcp",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", APRL_PROFILE_INTERVAL="0",
               HDF5_USE_FILE_LOCKING="FALSE",
               MPLCONFIGDIR=str(root / ("mpl-cache-" + str(task))))
    srun = ["srun", "--mpi=" + env.get("SRUN_MPI", "pmix_v3"), "--exact",
            "--nodes=1", "--ntasks=16", "--ntasks-per-node=16",
            "--cpus-per-task=8", "--cpu-bind=cores", "--kill-on-bad-exit=1",
            "--wait=30"]
    binding_log = root / f"binding-{task}.log"
    launch(srun + [campaign["binding_probe"], "16", "8"], root, env, 60,
           binding_log.name)
    save(root / f"BINDING_{task}.json", binding(binding_log, 1))
    results = []
    for mode in campaign["modes"]:
        for kind in ["straight", "bend"]:
            case = next(item for item in campaign["cases"]
                        if item["shape"] == shape and item["mode"] == mode and
                        item["kind"] == kind)
            directory = Path(case["directory"])
            plan = json.loads((directory / "PLAN.json").read_text())
            require(plan["mode"] == mode, "Prepared mode differs")
            record = launch(
                srun + [campaign["binary"], str(directory / "run.yaml")],
                directory, env, case["timeout_seconds"])
            checks = verify_case(directory, plan)
            complete = dict(complete=True, shape=shape, mode=mode, kind=kind,
                            job_id=os.environ["SLURM_JOB_ID"], launch=record,
                            checks=checks)
            save(directory / "PIC_COMPLETE.json", complete)
            results.append(complete)
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed during run")
    complete = dict(complete=True, task=task, shape=shape,
                    job_id=os.environ["SLURM_JOB_ID"], cases=results)
    save(root / f"ARRAY_COMPLETE_{task}.json", complete)
    print(json.dumps(complete, indent=2), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
