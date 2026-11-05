#!/usr/bin/env bash
#SBATCH --job-name=critical-power
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
power_root=$(realpath "${1:?Prepared power review directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR="$power_root/mpl-cache"
ulimit -c 0
"$PYTHON_BIN" "$power_root/scripts/power_review.py" "$power_root" > "$power_root/driver.log" 2>&1
