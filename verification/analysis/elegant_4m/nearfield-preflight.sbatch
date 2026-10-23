#!/usr/bin/env bash
#SBATCH --job-name=nearfield-preflight
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -Eeuo pipefail
root=$(realpath "${1:?Result directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
cd "$root"
sha256sum --check preflight.sha256
"${step[@]}" "$PYTHON_BIN" -c 'import os; a=os.sched_getaffinity(0); print("Allowed CPUs:",sorted(a)); assert len(a)==1' > binding-preflight.log 2>&1
"${step[@]}" /public/home/wfollowm/APRL_Workspace/build/verification/postprocess/field_reconstruction/field_reconstruction "$root/small-reconstruction.yaml" > small-reconstruction.log 2>&1
"$PYTHON_BIN" - "$root/small-reconstructed.h5" <<'PY'
import h5py,sys,numpy as np
with h5py.File(sys.argv[1],'r') as f:
    g=f['reconstructed_field']; assert g['complete'][()]==1
    assert np.isfinite(g['electric_V_per_m'][:]).all()
    assert np.isfinite(g['magnetic_T'][:]).all()
    assert np.isfinite(g['cleaned_signed_power_W'][:]).all()
print('Small existing ballistic reconstruction fixture PASS')
PY
run=/public/home/wfollowm/APRL_Workspace/results/elegant_4m_analysis_20261024/original/runs/20261021T070618Z-861230
"${step[@]}" "$PYTHON_BIN" "$root/scripts/preflight_reconstruction.py" \
  --field "$run/detectors/exit-field.h5" \
  --particles "$run/detectors/exit-field-ballistic-reference.h5" \
  --output "$root/projection-preflight.json" > projection-preflight.log 2>&1
