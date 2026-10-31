#!/bin/bash
#SBATCH --job-name=work-geometry-postcheck
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:06:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?postcheck root}")
campaign=$(realpath "${2:?closed campaign}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export HDF5_USE_FILE_LOCKING=FALSE
timeout --signal=TERM --kill-after=10 330 "$PYTHON_BIN" "$root/scripts/complete_geometry_campaign.py" "$root" "$campaign" > "$root/driver.log" 2>&1
