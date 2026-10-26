#!/usr/bin/env bash
# Bounded regression suite; no production input or layout/performance scan.
set -Eeuo pipefail
[[ $# == 2 ]] || { echo 'usage: run_gnode.sh BUILD RESULT'; exit 2; }
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || exit 2
[[ -z ${SLURM_JOB_ID:-} ]] || exit 2
build=$(realpath "$1")
result=$(realpath "$2")
[[ ! -e "$result/STARTED" ]] || { echo 'Use a fresh result directory'; exit 2; }
source "$workspace/cluster.env"
# Shared copies of the build host FFTW libraries, as in accepted regressions.
[[ -f "$build/runtime-libs/libfftw3.so.3" && -f "$build/runtime-libs/libfftw3_threads.so.3" ]] || exit 2
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU test holds the lock'; exit 2; }
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$result/EXIT_CODE"; date --iso-8601=seconds > "$result/ENDED"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$result/STARTED"
hostname > "$result/host.txt"
timeout --signal=TERM --kill-after=10s 60s mpiexec -n 4 --map-by slot:PE=2 \
  --bind-to core --nooversubscribe --report-bindings \
  env OMP_NUM_THREADS=2 "$workspace/build/node-survey/binding_probe" 4 2 \
  > "$result/binding.log" 2>&1
grep -q '\[binding\] result=PASS failures=0' "$result/binding.log"
timeout --signal=TERM --kill-after=20s 900s taskset -c 0-15 \
  ctest --test-dir "$build" -L required --output-on-failure -j 1 \
  > "$result/ctest.log" 2>&1
cp -r "$build/Testing" "$result/Testing"
