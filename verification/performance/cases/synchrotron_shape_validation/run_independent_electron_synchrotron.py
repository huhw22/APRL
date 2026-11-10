"""Run or safely resume the independent-electron bend/straight pair."""
import json
import os
from pathlib import Path
import re
import sys

import h5py

from launcher import launch, binding
from prepare_paired import require, save, sha


def verify_case(directory, plan):
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
    completed_steps = []
    for rank in range(16):
        path = directory / "reference" / f"one-way-source-rank-{rank:05d}.csv"
        words = path.with_suffix(".complete").read_text().split()
        require(len(words) == 2 and words[0] == "1", "Unclosed one-way rank")
        completed_steps.append(int(words[1]))
    require(len(set(completed_steps)) == 1 and
            0 < completed_steps[0] <= plan["expected_steps"],
            "Inconsistent actual one-way step count")
    return {"field_samples": samples,
            "actual_completed_steps": completed_steps[0]}


def main(root):
    root = Path(root).resolve()
    campaign = json.loads((root / "campaign.json").read_text())
    frozen = json.loads((root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed")
    require(os.environ.get("SLURM_NTASKS") == "16", "Expected 16 MPI ranks")
    require(sha(campaign["binary"]) == campaign["binary_sha256"], "Binary changed")
    env = {key: value for key, value in os.environ.items() if not key.startswith("APRL_")}
    env.update(OMP_NUM_THREADS="8", OMP_DYNAMIC="FALSE", OMP_PLACES="cores",
               OMP_PROC_BIND="close", OMPI_MCA_btl="self,vader,tcp",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", APRL_PROFILE_INTERVAL="0",
               HDF5_USE_FILE_LOCKING="FALSE")
    srun = ["srun", "--mpi=" + env.get("SRUN_MPI", "pmix_v3"), "--exact",
            "--nodes=1", "--ntasks=16", "--ntasks-per-node=16",
            "--cpus-per-task=8", "--cpu-bind=cores", "--kill-on-bad-exit=1", "--wait=30"]
    launch(srun + [campaign["binding_probe"], "16", "8"], root, env, 60,
           "binding-resume.log")
    save(root / "BINDING_RESUME.json", binding(root / "binding-resume.log", 1))
    results = []
    for case in campaign["cases"]:
        directory = Path(case["directory"])
        plan = json.loads((directory / "PLAN.json").read_text())
        if (directory / "aprl.log").is_file():
            verification = verify_case(directory, plan)
            results.append({"kind": case["kind"], "launch": None,
                            "reused_completed_case": True, **verification})
            continue
        record = launch(srun + [campaign["binary"], str(directory / "run.yaml")],
                        directory, env, case["timeout_seconds"])
        verification = verify_case(directory, plan)
        results.append({"kind": case["kind"], "launch": record,
                        "reused_completed_case": False, **verification})
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed during run")
    save(root / "PIC_COMPLETE.json", {"complete": True, "cases": results,
                                       "job_id": os.environ["SLURM_JOB_ID"],
                                       "resumed_after_validator_only_failure": True})
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
