#!/usr/bin/env bash
#SBATCH --job-name=aprl-discrete-mpi
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen MPI campaign}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/STARTED" ]]
date --iso-8601=seconds > "$root/STARTED"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/EXIT_CODE"; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
timeout --signal=TERM --kill-after=20 1420 "$PYTHON_BIN" "$root/scripts/mpi.py" "$root" > "$root/driver.log" 2>&1
