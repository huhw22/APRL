#!/usr/bin/env bash
#SBATCH --job-name=aprl-pending-test
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:40:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen validation campaign}")
stage=${2:?science or full}
[[ "$stage" == science || "$stage" == full ]]
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 1 && ${SLURM_NTASKS:-0} == 16 && ${SLURM_CPUS_PER_TASK:-0} == 8 ]]
mkdir "$root/$stage/execution"
cd "$root/$stage/execution"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > ENDED_AT; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
date --iso-8601=seconds > STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
# Driver stops actively at1800s; this outer limit preserves time for SIGTERM/KILL
# cleanup and remains below Slurm's40min hard ceiling. Cancellation is failure.
timeout --signal=TERM --kill-after=20 1840 "$PYTHON_BIN" "$root/scripts/run.py" --root "$root" --stage "$stage" > driver.log 2>&1
