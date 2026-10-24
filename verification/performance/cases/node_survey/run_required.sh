#!/usr/bin/env bash
# Required regression on the explicitly authorized, unmanaged CPU test host.
set -Eeuo pipefail
[[ $# == 2 ]] || { echo 'Usage: run_required.sh BUILD_DIRECTORY NEW_RESULT_DIRECTORY' >&2; exit 2; }
build_dir=$(realpath "$1")
result_dir=$2
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || {
  echo 'Requires gnode1 and explicit unscheduled CPU-test authorization' >&2; exit 2;
}
[[ -z ${SLURM_JOB_ID:-} && "$result_dir" == /* && ! -e "$result_dir" ]] || {
  echo 'Use a fresh absolute result directory; do not fabricate Slurm variables' >&2; exit 2;
}
[[ "$build_dir" == "$workspace/build/"* && -x "$build_dir/aprl" ]] || {
  echo 'Requires an independent verification build inside the workspace' >&2; exit 2;
}
grep -Fxq 'MPIEXEC_EXECUTABLE:FILEPATH=/usr/bin/mpiexec' "$build_dir/CMakeCache.txt"
grep -Fxq 'MPIEXEC_PREFLAGS:STRING=--host;gnode1:16;--map-by;slot:PE=2;--bind-to;core;--nooversubscribe' "$build_dir/CMakeCache.txt"
source "$workspace/cluster.env"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU test holds the calibration lock' >&2; exit 2; }
awk '{exit !($1 <= 2)}' /proc/loadavg || {
  echo 'gnode1 has nontrivial load; inspect resources before starting regression' >&2; exit 2;
}
mkdir "$result_dir"
mkdir "$build_dir/.gnode1-required-started"
finish() {
  code=$?
  trap - EXIT
  printf '%s\n' "$code" > "$result_dir/EXIT_CODE"
  if (( code == 0 )); then
    printf 'All registered required tests passed on gnode1\n' > "$result_dir/PASS"
  else
    printf 'Required regression failed with exit %s\n' "$code" > "$result_dir/FAIL"
  fi
  exit "$code"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
exec >"$result_dir/driver.log" 2>&1
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
date --iso-8601=seconds
hostname
lscpu --parse=CPU,CORE,SOCKET,NODE,ONLINE > "$result_dir/cpu-topology.csv"
sha256sum "$build_dir/aprl" "$build_dir/CMakeCache.txt" > "$result_dir/build.sha256"
cp "$build_dir/CMakeCache.txt" "$result_dir/CMakeCache.txt"
OMP_NUM_THREADS=4 timeout --signal=TERM --kill-after=10 60 \
  mpirun --host gnode1:16 -np 4 --map-by slot:PE=4 --bind-to core \
  --nooversubscribe --report-bindings \
  "$workspace/build/node-survey/binding_probe" 4 4 > "$result_dir/binding.log" 2>&1
grep -q '^\[binding\] result=PASS failures=0' "$result_dir/binding.log"
# A new build owns its CTest output paths. Single-process tests may create up
# to four threads; MPI tests use the explicit CPU mapping recorded in its cache.
timeout --signal=TERM --kill-after=20 1200 taskset -c 0-15 \
  ctest --test-dir "$build_dir" -L required -j 1 --stop-on-failure --no-tests=error --output-on-failure
cp -a "$build_dir/Testing" "$result_dir/Testing"
date --iso-8601=seconds
