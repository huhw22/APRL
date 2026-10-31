#!/bin/bash
#SBATCH --job-name=coupling-time-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -euo pipefail
root=${1:?output root required}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --mpi=none --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
  /public/home/wfollowm/micromamba/envs/physics/bin/python3 "$root/scripts/analyze_time_pairing.py" \
  --workspace /public/home/wfollowm/APRL_Workspace --out "$root/time_pairing" --workers 4
