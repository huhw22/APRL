#!/usr/bin/env bash
#SBATCH --job-name=aprl-SPARC-SASE-ABC
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=144:00:00
#SBATCH --no-requeue
#SBATCH --array=0-2%1
set -Eeuo pipefail
root="$1"
index="$SLURM_ARRAY_TASK_ID"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$root"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$root/mpl-cache"
ulimit -s unlimited
finish(){ code=$?; trap - EXIT; printf '%s\n' "$code" > "CASE_${index}_EXIT_CODE"; date --iso-8601=seconds > "CASE_${index}_ENDED_AT"; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
date --iso-8601=seconds > "CASE_${index}_STARTED_AT"
scontrol show job "$SLURM_JOB_ID" > "allocation-${index}.txt"
timeout --signal=TERM --kill-after=30 514800 "$PYTHON_BIN" scripts/run_main.py "$index" > "driver-${index}.log" 2>&1
