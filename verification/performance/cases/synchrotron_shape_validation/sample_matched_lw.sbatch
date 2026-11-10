#!/bin/bash
#SBATCH --job-name=synch-sample-lw
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --no-requeue

set -euo pipefail

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/aprl-sample-lw-mpl"

exec /public/home/wfollowm/micromamba/envs/physics/bin/python \
  /public/home/wfollowm/APRL_Workspace/source/aprl/verification/performance/cases/synchrotron_shape_validation/analyze_sample_matched_lw.py \
  "$@"
