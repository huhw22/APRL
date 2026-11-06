#!/usr/bin/env bash
#SBATCH --job-name=fel-one4one-binomial3
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=73:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
ulimit -c 0
fel_root=$(realpath "${1:?Prepared FEL root}")
timeout --signal=TERM --kill-after=30 260400 "$PYTHON_BIN" "$fel_root/scripts/run_fel.py" "$fel_root" > "$fel_root/driver.log" 2>&1
