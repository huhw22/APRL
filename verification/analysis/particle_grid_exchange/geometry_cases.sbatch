#!/bin/bash
#SBATCH --job-name=particle-work-geometry
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
timeout --signal=TERM --kill-after=20 1440 "$PYTHON_BIN" "$root/scripts/run_geometry_cases.py" "$root" > "$root/driver.log" 2>&1
