#!/usr/bin/env bash
# Shared initialization for manually submitted, single-node acceptance jobs.
set -Eeuo pipefail

fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

init_run() {
  local stage="$1"
  [[ -n "${SLURM_JOB_ID:-}" ]] || fail 'Submit this script yourself with sbatch; no login-node execution.'
  [[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}}" == 1 ]] || fail 'This acceptance stage requires one node.'
  workspace="${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}"
  [[ "$workspace" == /* && -f "$workspace/cluster.env" ]] || fail 'Set APRL_WORKSPACE to the prepared workspace.'
  source "$workspace/cluster.env"
  export APRL_VERIFICATION_BUILD_DIR="$workspace/build/verification"
  local prepared="$workspace/tests/cluster_acceptance"
  [[ -f "$prepared/manifest.sha256" && -f "$prepared/build-record.json" ]] || fail 'Build manifest is missing.'
  mkdir -p "$workspace/results/cluster_acceptance"
  run_dir="$(mktemp -d "$workspace/results/cluster_acceptance/${stage}-${SLURM_JOB_ID}-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")"
  log_dir="$run_dir"
  export workspace run_dir log_dir
  exec > >(tee "$run_dir/driver.log") 2>&1
  printf 'Result directory: %s\n' "$run_dir"
  printf '%s\n' "$stage" > "$run_dir/STAGE"
  cp "$prepared/manifest.sha256" "$prepared/build-record.json" "$run_dir/"
  { date -u --iso-8601=seconds; hostname; } > "$run_dir/preflight-host.txt"
  sha256sum --check --quiet "$run_dir/manifest.sha256" || fail 'Prepared artifacts changed; rebuild/review before running.'
  "$PYTHON_BIN" "$prepared/check-runtime.py" \
    "$run_dir/build-record.json" "$run_dir/runtime-audit.json"
  export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
  export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  export APRL_PROFILE_INTERVAL=0
  {
    date -u --iso-8601=seconds
    hostname
    lscpu
    scontrol show job "$SLURM_JOB_ID"
    scontrol show hostnames "$SLURM_JOB_NODELIST"
    mpirun --version
    srun --mpi=list
    for executable in "$APRL_BASELINE" "$APRL_CANDIDATE" "$SASE_GENERATOR"; do
      printf '\nExecutable: %s\n' "$executable"
      ldd "$executable"
    done
    "$PYTHON_BIN" - <<'PY'
import json, os, sys, h5py, numpy
keys = ["SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK",
        "SLURM_JOB_CPUS_PER_NODE", "SLURM_CPUS_ON_NODE",
        "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PLACES", "OMP_PROC_BIND", "OMP_THREAD_LIMIT",
        "APRL_PROFILE_INTERVAL", "SRUN_MPI", "LD_LIBRARY_PATH",
        "OMPI_MCA_btl", "OMPI_MCA_pml", "UCX_TLS", "UCX_NET_DEVICES"]
print(json.dumps({"python": sys.executable, "numpy": numpy.__version__, "h5py": h5py.__version__,
                  "h5py_hdf5": h5py.version.hdf5_version,
                  "environment": {k: os.environ.get(k) for k in keys},
                  "batch_affinity": sorted(os.sched_getaffinity(0))}, indent=2))
PY
  } > "$run_dir/environment.log" 2>&1
}

require_pass() {
  local result="${1:-}" stage="$2"
  [[ "$result" == /* && -f "$result/PASS" && -f "$result/STAGE" ]] || fail "Missing completed $stage result directory."
  [[ "$(cat "$result/STAGE")" == "$stage" ]] || fail "Expected $stage evidence: $result"
  cmp -s "$result/manifest.sha256" "$run_dir/manifest.sha256" || fail 'Evidence uses different binaries or scripts.'
  printf '%s\n' "$result" > "$run_dir/${stage}-evidence.txt"
}
