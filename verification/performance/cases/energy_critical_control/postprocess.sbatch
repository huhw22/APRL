#!/usr/bin/env bash
#SBATCH --job-name=critical-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:08:00
#SBATCH --no-requeue
set -Eeuo pipefail
review_root=$(realpath "${1:?Prepared postprocessing directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
ulimit -c 0
"$PYTHON_BIN" "$review_root/scripts/postprocess.py" "$review_root" > "$review_root/driver.log" 2>&1
