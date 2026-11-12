#!/bin/bash
#SBATCH --job-name=puffin-pic-source-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/source-propagation-${SLURM_JOB_ID}
ulimit -c 0
root="$1"
timeout --signal=TERM --kill-after=30 1140 "$PYTHON_BIN" -B "$root/scripts/analyze.py" "$root" > "$root/driver.log" 2>&1
