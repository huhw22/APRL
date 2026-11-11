#!/bin/bash
#SBATCH --job-name=SPARC6-ledger-preflight
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -euo pipefail
root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1
"$PYTHON_BIN" "$root/scripts/build_preflight.py" "$root" > "$root/preflight_driver.log" 2>&1
