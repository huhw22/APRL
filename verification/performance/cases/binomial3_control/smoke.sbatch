#!/usr/bin/env bash
#SBATCH --job-name=aprl-binomial3-smoke
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
smoke_root=$(realpath "${1:?Prepared smoke root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
ulimit -c 0
timeout --signal=TERM --kill-after=20 1740 "$PYTHON_BIN" "$smoke_root/scripts/run_smoke.py" "$smoke_root" > "$smoke_root/driver.log" 2>&1
