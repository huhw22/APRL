#!/usr/bin/env bash
#SBATCH --job-name=GENESIS-SPARC-SASE
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$root"
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
finish(){ code=$?; trap - EXIT; printf '%s\n' "$code" > GENESIS_EXIT_CODE; date --iso-8601=seconds > GENESIS_ENDED_AT; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
date --iso-8601=seconds > GENESIS_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > genesis-allocation.txt
timeout --signal=TERM --kill-after=30 14900 "$PYTHON_BIN" scripts/run_reference.py > genesis-driver.log 2>&1
