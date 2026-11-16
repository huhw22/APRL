#!/bin/bash
#SBATCH --job-name=w355-F3-energy-curve-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/w355-two-method-${SLURM_JOB_ID}
report_root="$1"
timeout --signal=TERM --kill-after=30 1140 "$PYTHON_BIN" -B "$report_root/three-point/scripts/review.py" "$report_root/../six-aprl" "$report_root/three-point" energy-binomial3 > "$report_root/three-point/analysis.log" 2>&1
