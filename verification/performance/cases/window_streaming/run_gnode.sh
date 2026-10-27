#!/usr/bin/env bash
# User-authorized, bounded CPU-only feasibility probes. No production jobs.
set -Eeuo pipefail
[[ $# == 2 || ( $# == 3 && $3 == --analytic-only ) ]] || { echo 'usage: run_gnode.sh BUILD FRESH_RESULT [--analytic-only]'; exit 2; }
mode=${3:-paired-probes}
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || exit 2
[[ -z ${SLURM_JOB_ID:-} ]] || exit 2
build=$(realpath "$1")
result=$(realpath "$2")
[[ -d "$result" && ! -e "$result/STARTED" ]] || { echo 'Use a fresh prepared result directory'; exit 2; }
source "$workspace/cluster.env"
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU test holds the lock'; exit 2; }
[[ -z $(pgrep -x aprl || true) ]] || { echo 'gnode1 already has APRL work'; exit 2; }
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$result/EXIT_CODE"; date --iso-8601=seconds > "$result/ENDED"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$result/STARTED"
hostname > "$result/host.txt"
sha256sum "$build/verification/aprl_pic_dependency_probe" \
  "$build/verification/aprl_window_streaming_probe" \
  "$build/verification/aprl_window_analytic_probe" "$build/libdirect_eb_core.a" \
  > "$result/executables.sha256"
launcher=(mpiexec -n 1 --map-by slot:PE=1 --bind-to core --nooversubscribe)
timeout --signal=TERM --kill-after=10s 60s "${launcher[@]}" \
  "$workspace/build/node-survey/binding_probe" 1 1 </dev/null \
  > "$result/binding.log" 2>&1
grep -Fq '[binding] result=PASS failures=0' "$result/binding.log"
if [[ $mode == --analytic-only ]]; then
  timeout --signal=TERM --kill-after=10s 180s "${launcher[@]}" \
    "$build/verification/aprl_window_analytic_probe" "$result/analytic" </dev/null \
    > "$result/analytic.log" 2> "$result/analytic.err"
  echo 'Independent analytic boundary probe completed; scientific review remains separate.'
  exit 0
fi
timeout --signal=TERM --kill-after=10s 180s "${launcher[@]}" \
  "$build/verification/aprl_pic_dependency_probe" </dev/null \
  > "$result/pic-dependency.json" 2> "$result/pic-dependency.err"
timeout --signal=TERM --kill-after=10s 180s "${launcher[@]}" \
  "$build/verification/aprl_window_streaming_probe" "$result/window" </dev/null \
  > "$result/window.log" 2> "$result/window.err"
echo 'Both bounded probes completed; scientific review remains separate.'
