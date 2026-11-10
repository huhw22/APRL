"""Run one member of the parallel single-electron control."""
import json
import os
from pathlib import Path
import re
import sys

import h5py

from launcher import binding, launch
from prepare import require, save, sha


def main(root, kind):
    root = Path(root).resolve()
    require(kind in ("straight", "bend"), "Unknown case")
    campaign = json.loads((root / "campaign.json").read_text())
    frozen = json.loads((root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed")
    require(os.environ.get("SLURM_NTASKS") == "16", "Expected 16 MPI ranks")
    case = next(value for value in campaign["cases"] if value["kind"] == kind)
    directory = Path(case["directory"])
    plan = json.loads((directory / "PLAN.json").read_text())
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("APRL_")}
    env.update(OMP_NUM_THREADS="8", OMP_DYNAMIC="FALSE", OMP_PLACES="cores",
               OMP_PROC_BIND="close", OMPI_MCA_btl="self,vader,tcp",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", APRL_PROFILE_INTERVAL="0",
               HDF5_USE_FILE_LOCKING="FALSE")
    srun = ["srun", "--mpi=" + env.get("SRUN_MPI", "pmix_v3"), "--exact",
            "--nodes=1", "--ntasks=16", "--ntasks-per-node=16",
            "--cpus-per-task=8", "--cpu-bind=cores", "--kill-on-bad-exit=1",
            "--wait=30"]
    bind_log = directory / "binding.log"
    launch(srun + [campaign["binding_probe"], "16", "8"], directory, env, 60,
           bind_log.name)
    save(directory / "BINDING.json", binding(bind_log, 1))
    record = launch(srun + [campaign["binary"], str(directory / "run.yaml")],
                    directory, env, case["timeout_seconds"])
    log = (directory / "aprl.log").read_text()
    for marker in ["Configured stop reached:", "phase=solve-complete",
                   "phase=mpi-finalize-complete",
                   "self_field_feedback is disabled", "[one-way-source]"]:
        require(marker in log, "Missing marker " + marker)
    require(not re.search(r"Simulation failure|Segmentation fault|DUE TO TIME LIMIT",
                          log), "Simulation error marker")
    samples = {}
    for plane in plan["field_planes"]:
        path = directory / "detectors" / (plane["name"] + ".h5")
        with h5py.File(path, "r") as handle:
            group = handle["field_plane"]
            count = int(group["committed_samples"][()])
            require(int(group["complete"][()]) == 1 and
                    count == len(group["time_s"]) == len(group["fields"]) and
                    count > 32, "Incomplete field plane " + plane["name"])
            samples[plane["name"]] = count
    steps = []
    for rank in range(16):
        complete = directory / "reference" / \
            f"one-way-source-rank-{rank:05d}.complete"
        words = complete.read_text().split()
        require(len(words) == 2 and words[0] == "1", "Unclosed source rank")
        steps.append(int(words[1]))
    require(len(set(steps)) == 1 and 0 < steps[0] <= plan["expected_steps"],
            "Inconsistent actual step count")
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed during run")
    save(directory / "PIC_COMPLETE.json", {
        "complete": True,
        "kind": kind,
        "job_id": os.environ["SLURM_JOB_ID"],
        "launch": record,
        "field_samples": samples,
        "actual_completed_steps": steps[0],
    })


if __name__ == "__main__":
    require(len(sys.argv) == 3,
            "usage: run_single_electron_case.py ROOT straight|bend")
    main(sys.argv[1], sys.argv[2])

