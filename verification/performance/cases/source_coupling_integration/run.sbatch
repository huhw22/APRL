#!/bin/bash
#SBATCH --job-name=aprl-source-preflight
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen campaign}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
ulimit -c 0
timeout --signal=TERM --kill-after=20 4000 "$PYTHON_BIN" "$root/scripts/run.py" "$root" > "$root/driver.log" 2>&1
