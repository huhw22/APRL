#!/bin/bash
#SBATCH --job-name=random6D-FEL
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=05:10:00
#SBATCH --array=0-2
#SBATCH --no-requeue
set -Eeuo pipefail
root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
ulimit -c 0
timeout --signal=TERM --kill-after=30 18400 "$PYTHON_BIN" -B "$root/scripts/run_case.py" "$root" "$SLURM_ARRAY_TASK_ID" > "$root/driver-$SLURM_ARRAY_TASK_ID.log" 2>&1
