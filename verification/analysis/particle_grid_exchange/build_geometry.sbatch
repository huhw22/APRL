#!/bin/bash
#SBATCH --job-name=work-geometry-build
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
root=$(realpath "${1:?Frozen root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
export MPLCONFIGDIR="$root/mplconfig"
ulimit -s unlimited
timeout --signal=TERM --kill-after=20 1740 "$PYTHON_BIN" "$root/scripts/build_geometry.py" "$root" > "$root/build.log" 2>&1
