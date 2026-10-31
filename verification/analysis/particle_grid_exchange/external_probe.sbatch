#!/bin/bash
#SBATCH --job-name=external-pusher-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --mpi=none --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
  /public/home/wfollowm/micromamba/envs/physics/bin/python3 "${1:?root}/scripts/run_external_probe.py" "$1"
