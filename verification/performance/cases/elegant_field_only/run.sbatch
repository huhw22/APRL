#!/usr/bin/env bash
#SBATCH --job-name=elegant-held-field
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=144:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
[[ $# == 1 ]] || exit 2
root=$(realpath "$1")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
[[ ! -e "$root/JOB_STARTED_AT" && ! -e "$root/run-manifest.yaml" ]] || { echo 'Use a fresh frozen production directory'; exit 2; }
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0 HDF5_USE_FILE_LOCKING=FALSE
# Preserve the previously accepted two-node environment, without system tuning.
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
ulimit -s unlimited
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/EXIT_CODE"; date --iso-8601=seconds > "$root/JOB_ENDED_AT"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
sha256sum --check frozen.sha256 > frozen-check.log
"$PYTHON_BIN" - "$root" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
p=Path(sys.argv[1]); plan=json.loads((p/'run-plan.json').read_text())
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
    return h.hexdigest()
assert (Path(plan['preflight_result'])/'PASS').is_file(), 'Held field-only smoke has not passed'
smoke=json.loads((Path(plan['preflight_result'])/'plan.json').read_text())
assert smoke['binary_sha256']==plan['binary_sha256'], 'Smoke used a different program'
assert plan.get('case_variant','undulator') in ('undulator','free-drift'), 'Unknown field-only case variant'
assert smoke.get('case_variant','undulator')==plan.get('case_variant','undulator'), 'Smoke tested a different field-only variant'
assert (plan['nodes'],plan['ranks'],plan['threads_per_rank'])==(2,256,1)
binary=Path(plan['binary']); assert binary.is_absolute() and os.access(binary,os.X_OK)
assert digest(binary)==plan['binary_sha256'], 'Frozen executable changed'
f=Path(plan['input_file']);st=f.stat()
assert (st.st_size,st.st_mtime_ns)==(plan['input_size_bytes'],plan['input_mtime_ns']), 'Particle input changed'
assert digest(f)==plan['input_sha256'], 'Particle input hash changed'
keys=['SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_NTASKS','SLURM_CPUS_PER_TASK',
      'OMP_NUM_THREADS','OMP_PLACES','OMP_PROC_BIND','APRL_PROFILE_INTERVAL',
      'OMPI_MCA_btl','SRUN_MPI','HDF5_USE_FILE_LOCKING']
(p/'execution-environment.json').write_text(json.dumps({'environment':{k:os.environ.get(k) for k in keys},
    'binary':str(binary),'binary_sha256_verified':plan['binary_sha256'],
    'input_sha256_verified':plan['input_sha256']},indent=2)+'\n')
PY
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1)
timeout --signal=TERM --kill-after=30s 120s "${step[@]}" --wait=30 \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 256 1 > binding.log 2>&1
grep -q '\[binding\] result=PASS failures=0' binding.log
echo "Geometry preflight $(date --iso-8601=seconds)"
timeout --signal=TERM --kill-after=30s 900s srun --mpi="$SRUN_MPI" --exact \
  --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --cpu-bind=cores \
  --kill-on-bad-exit=1 --wait=30 "$PYTHON_BIN" "$root/scripts/geometry.py" "$root" > geometry.log 2>&1
"$PYTHON_BIN" - "$root" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);plan=json.loads((p/'run-plan.json').read_text())
g=json.loads(Path(plan['geometry_report']).read_text())
assert g.get('pass',g.get('pass_',False)), 'Geometry preflight failed'
(p/'PRODUCTION_READY').write_text('Frozen identity, held field-only smoke, binding and full-input geometry passed\n')
PY
binary=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["binary"])' "$root/run-plan.json")
echo "Production start $(date --iso-8601=seconds) binary=$binary"
"${step[@]}" --wait=300 "$binary" "$root/run.yaml" > aprl.log 2>&1
echo "Completion audit $(date --iso-8601=seconds)"
"$PYTHON_BIN" "$root/scripts/check_completion.py" "$root"
echo "Finished $(date --iso-8601=seconds)"
