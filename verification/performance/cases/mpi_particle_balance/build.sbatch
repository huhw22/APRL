#!/usr/bin/env bash
#SBATCH --job-name=aprl-mpi-balance-build
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:45:00
#SBATCH --no-requeue
set -Eeuo pipefail
balance_build=$(realpath "${1:?Frozen build root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
ulimit -s unlimited
ulimit -c 0
timeout --signal=TERM --kill-after=20 2640 "$PYTHON_BIN" "$balance_build/scripts/build.py" "$balance_build" > "$balance_build/build.log" 2>&1
