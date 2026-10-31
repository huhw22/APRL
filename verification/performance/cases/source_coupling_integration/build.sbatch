#!/bin/bash
#SBATCH --job-name=aprl-source-build
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen build root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
ulimit -s unlimited
ulimit -c 0
timeout --signal=TERM --kill-after=20 1740 "$PYTHON_BIN" "$root/scripts/build.py" "$root" > "$root/build.log" 2>&1
