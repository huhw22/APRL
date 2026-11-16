#!/bin/bash
#SBATCH --job-name=w355-MM-existing-data-audit
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
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
audit_root="$1"
timeout --signal=TERM --kill-after=30 1140 "$PYTHON_BIN" -B "$audit_root/scripts/audit.py" "$audit_root/.." "$audit_root" > "$audit_root/analysis.log" 2>&1
