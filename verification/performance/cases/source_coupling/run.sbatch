#!/bin/bash
#SBATCH --job-name=source-coupling-prototype
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?frozen campaign}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export HDF5_USE_FILE_LOCKING=FALSE
timeout --signal=TERM --kill-after=20 1440 "$PYTHON_BIN" "$root/scripts/run.py" "$root" > "$root/driver.log" 2>&1
