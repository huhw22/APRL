#!/usr/bin/env bash
set -Eeuo pipefail
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || {
  echo 'Requires gnode1 and APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1' >&2; exit 2;
}
[[ -z ${SLURM_JOB_ID:-} ]] || { echo 'Do not fabricate Slurm allocation variables' >&2; exit 2; }
source "$workspace/cluster.env"
# Both variants use the same two accepted, read-only FFTW runtime libraries.
# prepare records their exact hashes; verify_frozen checks them before launch.
followup_runtime="$workspace/build/cowan-accepted-r2-20261025/runtime-libs"
for followup_library in libfftw3.so.3 libfftw3_threads.so.3; do
  [[ -r "$followup_runtime/$followup_library" ]] || {
    echo "Missing accepted runtime library: $followup_runtime/$followup_library" >&2; exit 2;
  }
done
export LD_LIBRARY_PATH="$followup_runtime${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU calibration holds the workspace lock' >&2; exit 2; }
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
scripts=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec timeout --signal=TERM --kill-after=20 7200 "$PYTHON_BIN" "$scripts/run.py" "$@"
