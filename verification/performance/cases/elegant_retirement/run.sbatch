#!/usr/bin/env bash
#SBATCH --job-name=elegant-retirement
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=144:00:00
set -Eeuo pipefail
root=$(realpath "$1")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ -n "$SLURM_JOB_ID" && "$SLURM_NTASKS" == 256 && "$SLURM_JOB_NUM_NODES" == 2 ]]
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0 HDF5_USE_FILE_LOCKING=FALSE
# Preserve the existing two-node ELEGANT transport choice; no system changes.
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
ulimit -s unlimited
cd "$root"
trap 'rc=$?; echo "$rc" > "$root/EXIT_CODE"' EXIT
sha256sum --check frozen.sha256
"$PYTHON_BIN" - "$root" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
p=Path(sys.argv[1]); plan=json.loads((p/'run-plan.json').read_text())
assert (Path(plan['preflight_result'])/'PASS').is_file()
f=Path(plan['input_file']); stat=f.stat()
assert (stat.st_size,stat.st_mtime_ns)==(plan['input_size_bytes'],plan['input_mtime_ns'])
h=hashlib.sha256()
with f.open('rb') as stream:
    for block in iter(lambda:stream.read(16*1024*1024),b''):h.update(block)
assert h.hexdigest()==plan['input_sha256']
keys=['SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_NTASKS','SLURM_CPUS_PER_TASK',
      'OMP_NUM_THREADS','OMP_PLACES','OMP_PROC_BIND','APRL_PROFILE_INTERVAL',
      'OMPI_MCA_btl','SRUN_MPI','HDF5_USE_FILE_LOCKING']
(p/'execution-environment.json').write_text(json.dumps(
    {'environment':{k:os.environ.get(k) for k in keys},
     'input_sha256_verified':h.hexdigest()},indent=2)+'\n')
PY
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
srun --mpi="$SRUN_MPI" --cpu-bind=cores --kill-on-bad-exit=1 --wait=60 \
  /public/home/wfollowm/APRL_Workspace/build/cluster-tools/binding_probe 256 1 \
  > "$root/binding.log" 2>&1
echo "Start $(date --iso-8601=seconds) job=$SLURM_JOB_ID"
binary=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["binary"])' "$root/run-plan.json")
srun --mpi="$SRUN_MPI" --cpu-bind=cores --kill-on-bad-exit=1 --wait=300 \
  "$binary" "$root/run.yaml" > "$root/aprl.log" 2>&1
"$PYTHON_BIN" "$root/scripts/check_completion.py" "$root"
echo "Finished $(date --iso-8601=seconds)"
