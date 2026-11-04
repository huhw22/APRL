#!/usr/bin/env bash
#SBATCH --job-name=SPARC-entrance-analysis
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
analysis_root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$analysis_root"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$analysis_root/mpl-cache"
scontrol show job "$SLURM_JOB_ID" > allocation.txt
timeout --signal=TERM --kill-after=20 500 srun --mpi=none --exact --nodes=1 --ntasks=1 \
  --cpus-per-task=1 --cpu-bind=cores "$PYTHON_BIN" scripts/entrance_ballistic.py \
  /public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103 \
  "$analysis_root/analysis" > analysis.log 2>&1
