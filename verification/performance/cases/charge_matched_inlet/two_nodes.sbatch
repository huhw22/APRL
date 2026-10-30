#!/usr/bin/env bash
#SBATCH --job-name=aprl-inlet-two
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Result root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
scontrol show job "$SLURM_JOB_ID" > "$root/two_nodes_200/allocation.txt"
timeout --signal=TERM --kill-after=20 800 "$PYTHON_BIN" "$root/scripts/run_two_nodes.py" "$root" > "$root/two_nodes_200/driver.log" 2>&1
