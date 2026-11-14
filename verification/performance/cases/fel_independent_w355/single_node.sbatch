#!/bin/bash
#SBATCH --job-name=puffin-11-4slip-w355-1node
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_io=romio321 HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
ulimit -c 0
root="$1"
timeout --signal=TERM --kill-after=60 258000 "$PYTHON_BIN" -B "$root/scripts/run_puffin.py" "$root" > "$root/driver.log" 2>&1
