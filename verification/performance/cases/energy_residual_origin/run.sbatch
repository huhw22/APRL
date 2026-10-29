#!/usr/bin/env bash
#SBATCH --job-name=aprl-free-work-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/STARTED" ]]
date --iso-8601=seconds > "$root/STARTED"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/EXIT_CODE"; date --iso-8601=seconds > "$root/ENDED"; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
timeout --signal=TERM --kill-after=20 1120 "$PYTHON_BIN" "$root/scripts/run.py" "$root" > "$root/driver.log" 2>&1
