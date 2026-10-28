#!/usr/bin/env bash
#SBATCH --job-name=aprl-moving-full
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=144:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen case directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 1 && ${SLURM_NTASKS:-0} == 16 && ${SLURM_CPUS_PER_TASK:-0} == 8 ]]
[[ ! -e "$root/JOB_STARTED_AT" && ! -e "$root/run-manifest.yaml" ]]
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > JOB_ENDED_AT; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
export HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
binary=$("$PYTHON_BIN" -c 'import json,sys;print(json.load(open(sys.argv[1]))["binary"])' run-plan.json)
export LD_LIBRARY_PATH="$(dirname "$binary")/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$PYTHON_BIN" - "$root" <<'PY' > identity-check.json
import hashlib,json,os,sys
from pathlib import Path
r=Path(sys.argv[1]);p=json.loads((r/'run-plan.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
for name,expected in json.loads((r/'frozen.json').read_text()).items():assert sha(name)==expected,name
passed=json.loads(Path(p['preflight_pass']).read_text())
assert passed['pass'] and passed['binary_sha256']==p['binary_sha256']
assert passed['card_sha256'][str(r/'run.yaml')]==sha(r/'run.yaml')
st=Path(p['input_file']).stat()
assert (st.st_size,st.st_mtime_ns)==(p['input_size_bytes'],p['input_mtime_ns'])
assert passed['input_sha256']==p['input_sha256']
print(json.dumps({'pass':True,'binary_sha256':p['binary_sha256'],'input_sha256_preflight':passed['input_sha256'],
 'unchanged_input_stat':True,'slurm_job_id':os.environ['SLURM_JOB_ID'],'nodes':os.environ['SLURM_JOB_NODELIST']}))
PY
step=(srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=16 --ntasks-per-node=16 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1)
timeout --signal=TERM --kill-after=30s 120s "${step[@]}" --wait=30 --chdir="$root" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 16 8 > binding.log 2>&1
grep -q '\[binding\] result=PASS failures=0' binding.log
date --iso-8601=seconds > PRODUCTION_STARTED_AT
"${step[@]}" --wait=660 --chdir="$root" "$binary" "$root/run.yaml" > aprl.log 2>&1
"$PYTHON_BIN" "$root/../scripts/check_completion.py" "$root" > completion.log 2>&1
