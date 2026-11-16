#!/bin/bash
#SBATCH --job-name=bend-two-method-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/bend-two-method-review-${SLURM_JOB_ID}
root="$1"
timeout --signal=TERM --kill-after=30 3540 "$PYTHON_BIN" -B "$root/scripts/analyze.py" /public/home/wfollowm/APRL_Workspace/results/synchrotron_full_ledger_20261116_r1 "$root" /public/home/wfollowm/APRL_Workspace > "$root/analysis.log" 2>&1
