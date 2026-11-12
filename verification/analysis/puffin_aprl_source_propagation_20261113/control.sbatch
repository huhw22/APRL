#!/bin/bash
#SBATCH --job-name=puffin-explicit-route
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_io=romio321 HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=/tmp/puffin-route-${SLURM_JOB_ID}
ulimit -s unlimited
ulimit -c 0
root="$1"
timeout --signal=TERM --kill-after=30 1140 "$PYTHON_BIN" -B "$root/scripts/control.py" "$root" > "$root/driver.log" 2>&1
