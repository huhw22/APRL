#!/usr/bin/env bash
#SBATCH --job-name=aprl-cowan-acceptance
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Pass new prepared acceptance directory}")
build=$(realpath "${2:?Pass independent verification build directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
[[ "$build" == "$APRL_WORKSPACE/build/"* && -x "$build/aprl" ]] || exit 2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0 HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
ulimit -s unlimited
cd "$root"
finish() {
  rc=$?
  trap - EXIT
  cp -a "$build/Testing" "$root/required/Testing" || rc=1
  echo "$rc" > "$root/JOB_EXIT_CODE"
  date --iso-8601=seconds > "$root/JOB_ENDED_AT"
  if (( rc == 0 )); then
    echo 'Required regression, reproduced old race and fixed 256-rank normal/bounded outputs accepted.' > "$root/ACCEPTANCE_PASS"
  fi
  exit "$rc"
}
mkdir "$build/.duration-acceptance-started"
mkdir -p required old-race
trap finish EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
scontrol show hostnames "$SLURM_JOB_NODELIST" > allocation-nodes.txt
sha256sum --check frozen.sha256 > frozen-verification.log
sha256sum --check build-frozen.sha256 > build-verification.log
OMP_NUM_THREADS=2 timeout --signal=TERM --kill-after=10 90 \
  srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=4 --ntasks-per-node=4 \
  --cpus-per-task=2 --cpu-bind=cores --kill-on-bad-exit=1 --wait=10 \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 4 2 > binding-required.log 2>&1
ctest_cpus=$("$PYTHON_BIN" - <<'PY'
import re
from pathlib import Path
s=Path('binding-required.log').read_text()
assert '[binding] result=PASS failures=0' in s
rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+).*?status=PASS',s)
assert len(rows)==8 and len({r[2] for r in rows})==1 and len({r[3] for r in rows})==8
print(','.join(sorted({r[3] for r in rows},key=int)))
PY
)
echo 'Reproduce pre-fix duration-close race with delayed root writer'
set +e
timeout --signal=TERM --kill-after=10 90 cmake \
  -DAPRL_EXECUTABLE="$APRL_WORKSPACE/build/cowan-interior-20261025/aprl" \
  -DDURATION_DELAY_LIBRARY="$build/verification/libaprl_duration_close_delay.so" \
  -DOUTPUT_CHECKER="$build/verification/aprl_output_checks" \
  -DDURATION_CARD="$APRL_SOURCE_DIR/verification/cards/duration-close.yaml" \
  -DDURATION_WORK="$root/old-race" \
  -P "$build/verification/check_duration_close.cmake" > old-race/driver.log 2>&1
old_rc=$?
set -e
echo "$old_rc" > old-race/DRIVER_EXIT_CODE
"$PYTHON_BIN" - "$old_rc" <<'PY'
import json,sys
from pathlib import Path
assert sys.argv[1]=='1', 'The deterministic old regression must fail, not time out'
s=Path('old-race/run.log').read_text()
assert 'DURATION_CLOSE_DELAY_BEGIN' in s and 'mesh.duration was exhausted' in s
assert 'DURATION_CLOSE_DELAY_COMPLETE' not in s
Path('old-race/REPRODUCED.json').write_text(json.dumps({'expected_failure_reproduced':True,'root_delayed_close_started':True,'root_delayed_close_completed':False})+'\n')
PY
echo "Run full required suite $(date --iso-8601=seconds)"
timeout --signal=TERM --kill-after=20 900 taskset -c "$ctest_cpus" \
  ctest --test-dir "$build" -L required -j 1 --stop-on-failure \
  --no-tests=error --output-on-failure --output-junit "$root/required/junit.xml" \
  > required/ctest.log 2>&1
echo 0 > required/EXIT_CODE
echo 'All registered required checks passed' > required/PASS
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30)
timeout --signal=TERM --kill-after=10 90 "${step[@]}" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 256 1 > binding-256.log 2>&1
for variant in normal-control normal-fixed bounded-1 bounded-2; do
  binary="$build/aprl"
  if [[ $variant == normal-control ]]; then binary="$APRL_WORKSPACE/build/candidate/aprl"; fi
  extra=()
  if [[ $variant == bounded-* ]]; then
    extra=(env "LD_PRELOAD=$build/verification/libaprl_duration_close_delay.so" APRL_TEST_DELAY_H5_CLOSE=1)
  fi
  echo "256-rank $variant $(date --iso-8601=seconds)"
  set +e
  timeout --signal=TERM --kill-after=20 240 "${step[@]}" "${extra[@]}" "$binary" "$root/$variant/run.yaml" > "$root/$variant/aprl.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$root/$variant/EXIT_CODE"
  [[ $rc != 124 && $rc != 137 ]] || exit "$rc"
  "$PYTHON_BIN" "$root/scripts/cases.py" check --root "$root" --variant "$variant" > "$root/$variant/check.log" 2>&1
done
"$PYTHON_BIN" "$root/scripts/cases.py" check --root "$root" > scientific-comparison.log 2>&1
echo "Accepted $(date --iso-8601=seconds)"
