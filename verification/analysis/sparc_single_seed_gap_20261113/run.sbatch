#!/bin/bash
#SBATCH --job-name=single-seed-gap
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
root="$1"
timeout --signal=TERM --kill-after=20 840 "$PYTHON_BIN" -B "$root/scripts/analyze.py" /public/home/wfollowm/APRL_Workspace/results/sparc_random_startup_20261112_r2 "$root/analysis" > "$root/driver.log" 2>&1
