"""Run, close, and analyze the two bend-only boundary-exit controls."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import h5py

from launcher import binding, launch
from prepare import require, save, sha


def verify_case(directory, plan):
    log = (directory / "aprl.log").read_text()
    for marker in ("Configured stop reached:", "phase=solve-complete",
                   "phase=mpi-finalize-complete",
                   "self_field_feedback is disabled", "[one-way-source]"):
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
    particle_counts = {}
    for name in ("entry", "bend-exit", "exit"):
        with h5py.File(directory / "detectors" / (name + ".h5"), "r") as handle:
            group = handle["particle_plane"]
            require(int(group["complete"][()]) == 1,
                    "Incomplete particle plane " + name)
            particle_counts[name] = int(group["committed_records"][()])
    expected = int(plan["expected_particles"])
    require(particle_counts["bend-exit"] >= .99 * expected,
            "Particles were lost before leaving the magnet")
    require(particle_counts["exit"] < .05 * expected,
            "Boundary-exit geometry did not remove the bent bunch")
    return {"field_samples": samples,
            "particle_plane_records": particle_counts,
            "exit_survival_fraction": particle_counts["exit"] / expected}


def main(root):
    root = Path(root).resolve()
    campaign = json.loads((root / "campaign.json").read_text())
    frozen = json.loads((root / "FROZEN.json").read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed")
    require(os.environ.get("SLURM_NTASKS") == "16", "Expected 16 MPI ranks")
    require(sha(campaign["binary"]) == campaign["binary_sha256"],
            "Binary changed")
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
    launch(srun + [campaign["binding_probe"], "16", "8"], root, env, 60,
           "binding.log")
    save(root / "BINDING.json", binding(root / "binding.log", 1))
    results = []
    for case in campaign["cases"]:
        directory = Path(case["directory"])
        plan = json.loads((directory / "PLAN.json").read_text())
        record = launch(srun + [campaign["binary"], str(directory / "run.yaml")],
                        directory, env, case["timeout_seconds"])
        results.append({"label": case["label"], "launch": record,
                        **verify_case(directory, plan)})
    require(all(sha(path) == digest for path, digest in frozen.items()),
            "Frozen input changed during run")
    save(root / "PIC_COMPLETE.json", {
        "complete": True, "cases": results,
        "job_id": os.environ["SLURM_JOB_ID"],
    })
    analysis = subprocess.run(
        [sys.executable, str(root / "scripts/analyze_100k_bend_boundary_exit.py"),
         str(root)], cwd=root, env=env, check=False)
    require(analysis.returncode == 0, "Boundary-exit analysis failed")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
