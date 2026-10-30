#!/usr/bin/env bash
#SBATCH --job-name=aprl-charge-scale
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:28:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen comparison root}")
case_name=${2:?Case name}
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/$case_name/STARTED" ]]
date --iso-8601=seconds > "$root/$case_name/STARTED"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/$case_name/EXIT_CODE"; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
scontrol show job "$SLURM_JOB_ID" > "$root/$case_name/allocation.txt"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
timeout --signal=TERM --kill-after=20 1620 "$PYTHON_BIN" "$root/scripts/run_case.py" "$root" "$case_name" > "$root/$case_name/driver.log" 2>&1
