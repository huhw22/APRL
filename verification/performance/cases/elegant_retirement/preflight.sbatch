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
workspace=/public/home/wfollowm/APRL_Workspace
source "$workspace/cluster.env"
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
# Same established transport as the historical two-node APRL run.
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,vader,tcp}"
export HDF5_USE_FILE_LOCKING=FALSE
mkdir -p "$root"
cd "$root"
if [[ -f "$root/frozen.sha256" ]]; then sha256sum --check "$root/frozen.sha256"; fi
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
timeout --signal=TERM --kill-after=30 120 "${step[@]}" "$workspace/build/cluster-tools/binding_probe" 256 1 > "$root/binding.log" 2>&1
srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 \
  "$PYTHON_BIN" "$scripts/geometry_preflight.py" \
  --input-file "$workspace/results/elegant_4m_analysis_20261024/original/particles.h5" \
  --historical-plan "$workspace/results/elegant_4m_analysis_20261024/original/runs/20261021T070618Z-861230/input-plan.json" \
  --output "$root/geometry.json" > "$root/geometry.log" 2>&1
"$PYTHON_BIN" "$scripts/preflight.py" prepare --root "$root" --binary "$APRL_CANDIDATE" > "$root/prepare.log" 2>&1
timeout --signal=TERM --kill-after=30 600 "${step[@]}" "$APRL_CANDIDATE" "$root/functional/run.yaml" > "$root/functional/aprl.log" 2>&1
"$PYTHON_BIN" "$scripts/preflight.py" check --root "$root" --binary "$APRL_CANDIDATE" > "$root/check.log" 2>&1
