#!/usr/bin/env bash
# Bounded CPU-only checks on the explicitly authorized small-test host.
set -euo pipefail
test "$(hostname -s)" = gnode1
test "${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-0}" = 1
test -z "${SLURM_JOB_ID:-}"
test "$#" = 2
check_build=$(realpath "$1")
check_result=$(realpath -m "$2")
test ! -e "$check_result"
cd /public/home/wfollowm/APRL_Workspace
source cluster.env
exec 9>results/node_survey/.cpu-calibration.lock
flock -n 9
if pgrep -x 'aprl|mpirun|mpiexec|orted' > /dev/null; then
  echo 'Existing MPI/compute process; refusing to overlap tests.' >&2
  exit 1
fi
mkdir -p "$check_result"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$check_build/verification/runtime-libs:${LD_LIBRARY_PATH:-}"
OMP_NUM_THREADS=2 timeout --signal=TERM --kill-after=5 45 mpirun -n 4 --map-by slot:PE=2 \
  --bind-to core --nooversubscribe build/node-survey/binding_probe 4 2 \
  > "$check_result/binding.log" 2>&1
grep -F '[binding] result=PASS failures=0' "$check_result/binding.log"
cd "$check_build"
oracle_status=0
timeout --signal=TERM --kill-after=15 600 ctest -j1 --output-on-failure \
  -L cowan-initializer > "$check_result/oracle-ctest.log" 2>&1 || oracle_status=$?
cp Testing/Temporary/LastTest.log "$check_result/oracle-detail.log"
test "$oracle_status" = 0
regression_status=0
timeout --signal=TERM --kill-after=15 1200 ctest -j1 --output-on-failure \
  -L required > "$check_result/required-ctest.log" 2>&1 || regression_status=$?
cp Testing/Temporary/LastTest.log "$check_result/required-detail.log"
test "$regression_status" = 0
printf 'Independent initializer and required regression checks passed.\n' > "$check_result/PASS"
