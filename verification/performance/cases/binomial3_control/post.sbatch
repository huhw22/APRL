#!/usr/bin/env bash
#SBATCH --job-name=fel-binomial3-growth
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
fel_root=$(realpath "${1:?Prepared FEL root}")
"$PYTHON_BIN" "$fel_root/scripts/analyze_fel.py" "$fel_root/one4one" >> "$fel_root/post.log" 2>&1
