#!/usr/bin/env bash
#SBATCH --job-name=aprl-FEL-3methods
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=05:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$root"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$root/mpl-cache"
ulimit -s unlimited
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > ENDED_AT; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
date --iso-8601=seconds > STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
timeout --signal=TERM --kill-after=20 17500 "$PYTHON_BIN" scripts/run_fel.py > driver.log 2>&1
