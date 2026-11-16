#!/bin/bash
#SBATCH --job-name=w355-MM-ledger-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/w355-ledger-review-${SLURM_JOB_ID}
review_root="$1"
timeout --signal=TERM --kill-after=30 1740 "$PYTHON_BIN" -B "$review_root/scripts/review.py" "$review_root/../six-aprl" "$review_root" > "$review_root/analysis.log" 2>&1
