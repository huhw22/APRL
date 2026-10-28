#!/usr/bin/env bash
#SBATCH --job-name=aprl-moving-preflight
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Prepared campaign}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
mkdir "$root/preflight"
cd "$root"
trap 'code=$?; printf "%s\n" "$code" > preflight/EXIT_CODE' EXIT
scontrol show job "$SLURM_JOB_ID" > preflight/allocation.txt
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
for case in undulator_4p5m free_drift_0p5m; do
  timeout --signal=TERM --kill-after=15 240 srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$PYTHON_BIN" "$root/scripts/geometry.py" --verify-plan "$root/$case/geometry-plan.json" \
    --output "$root/preflight/$case-geometry.json" > "$root/preflight/$case-geometry.log" 2>&1
done
timeout --signal=TERM --kill-after=30 900 "$PYTHON_BIN" "$root/scripts/preflight.py" --root "$root" > preflight/smoke-driver.log 2>&1
"$PYTHON_BIN" - "$root" <<'PY'
import hashlib,json,sys
from pathlib import Path
r=Path(sys.argv[1]);plan=json.loads((r/'plan.json').read_text())
smoke=json.loads((r/'preflight/SMOKE_PASS.json').read_text());assert smoke['pass']
geos={name:json.loads((r/'preflight'/f'{name}-geometry.json').read_text()) for name in plan['cases']}
for item in geos.values(): assert item['pass']
cards={str(Path(path)/'run.yaml'):hashlib.sha256((Path(path)/'run.yaml').read_bytes()).hexdigest() for path in plan['cases'].values()}
inputs={json.loads((Path(path)/'run-plan.json').read_text())['input_sha256'] for path in plan['cases'].values()}
assert len(inputs)==1
(r/'preflight/PASS.json').write_text(json.dumps({'pass':True,'binary_sha256':plan['binary_sha256'],
 'input_sha256':next(iter(inputs)),'card_sha256':cards,'geometry_reports':list(geos),
 'scope':'Full-input kinematics/reference bounds and 16x8 short functional smoke; not full-scale speed or radiation convergence.'},indent=2)+'\n')
PY
