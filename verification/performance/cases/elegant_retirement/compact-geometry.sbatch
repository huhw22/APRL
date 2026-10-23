#!/usr/bin/env bash
#SBATCH --job-name=elegant-compact-geometry
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -Eeuo pipefail
root=$(realpath "$1")
scripts=$(realpath "$2")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$root"
sha256sum --check frozen.sha256
scontrol show job "$SLURM_JOB_ID" > allocation.txt
srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 \
  --cpu-bind=cores --kill-on-bad-exit=1 \
  "$PYTHON_BIN" "$scripts/compact_geometry.py" \
  --workspace /public/home/wfollowm/APRL_Workspace --output-dir "$root" \
  --retirement-length 0.3 > geometry.log 2>&1
