#!/usr/bin/env bash
#SBATCH --job-name=aprl-transport-unit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
build=$(realpath "${1:?candidate build}")
root=$(realpath "${2:?new result directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 1 && ${SLURM_NTASKS:-0} == 16 && ${SLURM_CPUS_PER_TASK:-0} == 8 ]]
[[ ! -e "$root/STARTED" ]]
date --iso-8601=seconds > "$root/STARTED"
trap 'code=$?; printf "%s\n" "$code" > "$root/EXIT_CODE"' EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
cd "$build"
ctest -L required --show-only=json-v1 > "$root/required-tests.json"
# This target is currently built but not registered as a required CTest.
# Stop for review if future registration could cause duplicate MPI launches.
"$PYTHON_BIN" - "$root/required-tests.json" <<'PY'
import json,sys
tests=json.load(open(sys.argv[1]))['tests']
assert tests, 'No required tests found'
for executable in ('aprl_boris_work_checks',):
 assert any(executable in ' '.join(t.get('command', [])) for t in tests), 'Missing new required module: '+executable
matches=[t['name'] for t in tests if 'moving_window_mpi_checks' in ' '.join(t.get('command',[]))]
assert not matches, 'Moving-window MPI checks now registered; remove duplicate manual runs: '+repr(matches)
assert not any('aprl_energy_plane_checks' in ' '.join(t.get('command', [])) for t in tests), 'Energy-plane MPI module now registered; review duplicate launch policy'
print('Required tests:',len(tests),'moving_window_mpi_checks: separate explicit MPI coverage')
PY
timeout --signal=TERM --kill-after=10 120 srun --mpi="${SRUN_MPI:-pmix_v3}" --exact \
  --nodes=1 --ntasks=2 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
  "$build/verification/aprl_moving_window_mpi_checks" "$root/moving-2rank.json" > "$root/moving-2rank.log" 2>&1
timeout --signal=TERM --kill-after=5 60 "$build/verification/aprl_pending_particle_statistics_checks" > "$root/pending-statistics.log" 2>&1
timeout --signal=TERM --kill-after=15 600 ctest -L required --output-on-failure -j 1 > "$root/required.log" 2>&1
timeout --signal=TERM --kill-after=10 120 srun --mpi="${SRUN_MPI:-pmix_v3}" --exact \
  --nodes=1 --ntasks=4 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
  "$build/verification/aprl_moving_window_mpi_checks" "$root/moving-4rank.json" > "$root/moving-4rank.log" 2>&1
for ranks in 2 4; do
  [[ ! -e "$root/energy-plane-${ranks}rank" ]]
  timeout --signal=TERM --kill-after=10 90 srun --mpi="${SRUN_MPI:-pmix_v3}" --exact \
    --nodes=1 --ntasks="$ranks" --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
    "$build/verification/aprl_energy_plane_checks" "$root/energy-plane-${ranks}rank" > "$root/energy-plane-${ranks}rank.log" 2>&1
done
"$PYTHON_BIN" - "$root" <<'PY'
import json,sys
from pathlib import Path
r=Path(sys.argv[1])
for ranks in (2,4):
 d=json.loads((r/f'moving-{ranks}rank.json').read_text())
 assert d['status']=='PASS' and d['ranks']==ranks
 d=json.loads((r/f'energy-plane-{ranks}rank/summary.json').read_text())
 assert d['status']=='PASS' and d['ranks']==ranks
 for key in ('energy_planes_independent', 'pure_B_endpoint_energy_invariant', 'terminal_clip_and_half_open_crossings', 'macro_weight_once', 'signed_transport_and_physical_aperture', 'field_step_integration_independent_of_write_rhythm', 'no_unavailable_gap_bridge'):
  assert d[key] is True, key
PY
"$PYTHON_BIN" "$APRL_SOURCE_DIR/verification/performance/cases/energy_transport_audit/test_review.py" > "$root/reader-tests.log" 2>&1
printf 'PASS\n' > "$root/PASS"
