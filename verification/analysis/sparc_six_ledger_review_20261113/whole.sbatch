#!/bin/bash
#SBATCH --job-name=sparc6-ledger-lorentz-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sparc6-whole-${SLURM_JOB_ID}
root="$1"
timeout --signal=TERM --kill-after=30 1740 "$PYTHON_BIN" -B "$root/whole-scripts/whole.py" "$root/.." "$root/whole" > "$root/whole.log" 2>&1
