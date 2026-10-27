#!/usr/bin/env bash
# Bounded release-at-t0 experiment. It does not enable a production window.
set -Eeuo pipefail
[[ $# == 5 ]] || { echo 'usage: run_startup.sh BUILD FRESH_RESULT MODE BETA GRID'; exit 2; }
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || exit 2
[[ -z ${SLURM_JOB_ID:-} ]] || exit 2
build=$(realpath "$1")
result=$(realpath "$2")
mode=$3
beta=$4
grid=$5
[[ -d "$result" && ! -e "$result/STARTED" ]] || { echo 'Use a fresh result directory'; exit 2; }
source "$workspace/cluster.env"
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU test holds the lock'; exit 2; }
[[ -z $(pgrep -x aprl || true) ]] || { echo 'gnode1 has APRL work'; exit 2; }
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$result/EXIT_CODE"; date --iso-8601=seconds > "$result/ENDED"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$result/STARTED"
hostname > "$result/host.txt"
sha256sum "$build/verification/aprl_startup_self_field_probe" "$build/libdirect_eb_core.a" > "$result/executables.sha256"
launcher=(mpiexec -n 1 --map-by slot:PE=1 --bind-to core --nooversubscribe)
timeout --signal=TERM --kill-after=10s 60s "${launcher[@]}" \
  "$workspace/build/node-survey/binding_probe" 1 1 </dev/null \
  > "$result/binding.log" 2>&1
grep -Fq '[binding] result=PASS failures=0' "$result/binding.log"
timeout --signal=TERM --kill-after=10s 180s "${launcher[@]}" \
  "$build/verification/aprl_startup_self_field_probe" "$result" "$mode" "$beta" "$grid" </dev/null \
  > "$result/probe.log" 2> "$result/probe.err"
