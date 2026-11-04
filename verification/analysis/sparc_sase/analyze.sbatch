#!/usr/bin/env bash
#SBATCH --job-name=SPARC-explicit-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --no-requeue
set -Eeuo pipefail
analysis_root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$analysis_root"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR="$analysis_root/mpl-cache"
scontrol show job "$SLURM_JOB_ID" > allocation.txt
timeout --signal=TERM --kill-after=20 2400 srun --mpi=none --exact --nodes=1 --ntasks=1 \
  --cpus-per-task=8 --cpu-bind=cores "$PYTHON_BIN" scripts/analyze_explicit.py \
  /public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103/compact-production-r1 \
  "$analysis_root" /public/home/wfollowm/APRL_Workspace/source/aprl > analysis.log 2>&1
