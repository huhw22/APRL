#!/bin/bash
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue
set -Eeuo pipefail
root="$1"
mode="$2"
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
trap 'status=$?; printf "%s\n" "$status" > "$root/EXIT_$mode"' EXIT
if [ "$mode" = preflight ]; then
  "$PYTHON_BIN" "$root/scripts/preflight.py" "$root" > "$root/preflight.log" 2>&1
else
  "$PYTHON_BIN" "$root/scripts/run.py" "$root" "$mode" > "$root/driver-$mode.log" 2>&1
fi
