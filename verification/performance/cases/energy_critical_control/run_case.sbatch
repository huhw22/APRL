#!/usr/bin/env bash
#SBATCH --job-name=aprl-critical-control
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
case_root=$(realpath "${1:?Prepared case directory}")
shift
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
export APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
ulimit -c 0
scontrol show job "$SLURM_JOB_ID" > "$case_root/allocation.txt"
timeout --signal=TERM --kill-after=20 3480 "$PYTHON_BIN" "$case_root/../scripts/run_case.py" \
  "$case_root" "$@" > "$case_root/driver.log" 2>&1
