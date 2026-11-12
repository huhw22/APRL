#!/bin/bash
#SBATCH --job-name=FEL-gain-origin
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:40:00
#SBATCH --no-requeue
set -euo pipefail
root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_btl=self,vader,tcp
export HDF5_USE_FILE_LOCKING=FALSE
"$PYTHON_BIN" "$root/scripts/run.py" "$root" > "$root/driver.log" 2>&1
