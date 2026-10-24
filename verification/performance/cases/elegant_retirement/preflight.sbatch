#!/usr/bin/env bash
#SBATCH --job-name=elegant-retirement-preflight
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -Eeuo pipefail
root=$(realpath "${1:?Pass a new preflight result directory}")
scripts=$(realpath "${2:?Pass the frozen scripts directory}")
geometry_source=${3:-}
workspace=/public/home/wfollowm/APRL_Workspace
source "$workspace/cluster.env"
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
# Same established transport as the historical two-node APRL run.
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,vader,tcp}"
export HDF5_USE_FILE_LOCKING=FALSE
export APRL_TEARDOWN_TRACE=all
mkdir -p "$root"
cd "$root"
if [[ -f "$root/frozen.sha256" ]]; then sha256sum --check "$root/frozen.sha256"; fi
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
timeout --signal=TERM --kill-after=30 120 "${step[@]}" "$workspace/build/cluster-tools/binding_probe" 256 1 > "$root/binding.log" 2>&1
if [[ -n "$geometry_source" ]]; then
  "$PYTHON_BIN" - "$geometry_source" "$root" <<'PY'
import hashlib,json,shutil,sys
from pathlib import Path
source,root=map(Path,sys.argv[1:])
g=json.loads((source/'geometry.json').read_text())
a=json.loads((source/'arrival-projection.json').read_text())
assert g['pass'] and a['pass'] and g['retirement_length_m']==0.3
assert g['arrival_policy']['quiet_tail_fs']==0 and g['mesh_geometry']['guard_alignment_cells']==1
for record,prefix in ((g,'input'),(a,'source')):
    f=Path(record[prefix+'_file']); st=f.stat()
    assert (st.st_size,st.st_mtime_ns)==(record[prefix+'_size_bytes'],record[prefix+'_mtime_ns'])
names=('geometry.json','arrival-projection.json','geometry-policy.json')
identity={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names}
for name in names:
    assert not (root/name).exists()
    shutil.copy2(source/name,root/name)
(root/'REUSED_GEOMETRY.json').write_text(json.dumps({'source':str(source.resolve()),'sha256':identity,'reason':'Reuse completed geometry scan while checking production-matched MPI teardown settings.'},indent=2)+'\n')
PY
else
  srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 \
    "$PYTHON_BIN" "$scripts/compact_geometry.py" \
    --workspace "$workspace" --output-dir "$root" \
    --retirement-length 0.3 > "$root/geometry.log" 2>&1
fi
"$PYTHON_BIN" "$scripts/preflight.py" prepare --root "$root" --binary "$APRL_CANDIDATE" > "$root/prepare.log" 2>&1
timeout --signal=TERM --kill-after=30 600 "${step[@]}" --error="$root/functional/rank-%t.stderr" "$APRL_CANDIDATE" "$root/functional/run.yaml" > "$root/functional/aprl.log" 2>&1
"$PYTHON_BIN" "$scripts/preflight.py" check --root "$root" --binary "$APRL_CANDIDATE" > "$root/check.log" 2>&1
