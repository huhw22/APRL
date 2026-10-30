#!/usr/bin/env bash
#SBATCH --job-name=aprl-charge-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen comparison root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/aprl-charge-scaling-mpl
"$PYTHON_BIN" "$root/postprocess/review_scaling.py" "$root" > "$root/postprocess/driver.log" 2>&1
