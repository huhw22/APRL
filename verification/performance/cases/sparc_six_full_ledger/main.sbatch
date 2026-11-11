#!/bin/bash
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --no-requeue
set -euo pipefail
root="$1"
case_name="$2"
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1
"$PYTHON_BIN" "$root/scripts/run_case.py" "$root" "$case_name" > "$root/driver-$case_name.log" 2>&1
