#!/bin/bash
#SBATCH --job-name=puffin-boundary-flow
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/puffin-boundary-flow-${SLURM_JOB_ID}
export OMPI_MCA_btl=self,vader,tcp
python=/public/home/wfollowm/micromamba/envs/physics/bin/python3
"$python" -B "$1/scripts/run.py" "$1" > "$1/run.log" 2>&1
