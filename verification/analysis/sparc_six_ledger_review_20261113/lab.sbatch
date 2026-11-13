#!/bin/bash
#SBATCH --job-name=sparc6-ledger-lab-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sparc6-lab-${SLURM_JOB_ID}
root="$1"
timeout --signal=TERM --kill-after=30 3540 "$PYTHON_BIN" -B "$root/lab-scripts/lab.py" "$root/.." "$root/lab" > "$root/lab.log" 2>&1
