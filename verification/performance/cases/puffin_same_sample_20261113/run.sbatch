#!/bin/bash
#SBATCH --job-name=puffin-w64-u07
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=/tmp/puffin-matplotlib-${SLURM_JOB_ID}
ulimit -s unlimited
ulimit -c 0
root="$1"
timeout --signal=TERM --kill-after=30 3540 "$PYTHON_BIN" -B "$root/scripts/run.py" "$root" > "$root/driver.log" 2>&1
