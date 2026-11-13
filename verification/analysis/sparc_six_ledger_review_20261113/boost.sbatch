#!/bin/bash
#SBATCH --job-name=sparc6-ledger-boost-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sparc6-ledger-${SLURM_JOB_ID}
root="$1"
timeout --signal=TERM --kill-after=30 1140 "$PYTHON_BIN" -B "$root/scripts/boost.py" "$root/.." "$root/boost" > "$root/boost.log" 2>&1
