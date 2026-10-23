#!/usr/bin/env bash
# Explicitly authorized gnode1-only, bounded CPU calibration.
set -Eeuo pipefail
workspace=/public/home/wfollowm/APRL_Workspace
[[ "$(hostname)" == gnode1 ]] || { echo 'gnode1 only' >&2; exit 1; }
[[ "${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-}" == 1 ]] || { echo 'Explicit gnode1 CPU authorization flag required' >&2; exit 1; }
[[ "${SLURM_JOB_ID:-}" == '' ]] || { echo 'Use the Slurm workflow inside Slurm' >&2; exit 1; }
source "$workspace/cluster.env"
[[ -x "$APRL_BASELINE" && -x "$APRL_CANDIDATE" ]] || { echo 'Missing prepared APRL binaries' >&2; exit 1; }
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another gnode1 calibration holds the lock' >&2; exit 1; }
load_one=$(awk '{print $1}' /proc/loadavg)
awk -v value="$load_one" 'BEGIN { exit !(value <= 2) }' || { echo "gnode1 is busy: load $load_one" >&2; exit 1; }
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
exec "$PYTHON_BIN" "$workspace/tests/cpu_calibration/run.py" "$@"
