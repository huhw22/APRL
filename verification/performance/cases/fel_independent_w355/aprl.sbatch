#!/bin/bash
#SBATCH --job-name=aprl-6-independent-w355
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --array=0-1%1
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
ulimit -s unlimited
ulimit -c 0
root="$1"
case_names=(energy-binomial3 energy-momentum)
case_name="${case_names[$SLURM_ARRAY_TASK_ID]}"
"$PYTHON_BIN" -B "$root/scripts/run_pic.py" "$root/six-aprl" "$case_name" > "$root/six-aprl/driver-$case_name.log" 2>&1
