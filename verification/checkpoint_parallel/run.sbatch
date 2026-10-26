#!/usr/bin/env bash
#SBATCH --job-name=aprl-checkpoint-io
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$1
cd "$root"
source /public/home/wfollowm/APRL_Workspace/cluster.env
# Reuse the accepted two-node job environment; no transport tuning is tested.
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > JOB_EXIT_CODE; date --iso-8601=seconds > JOB_ENDED_AT; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
timeout --signal=TERM --kill-after=15s 90s srun --mpi=pmix_v3 --exact \
  --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 \
  --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
  /public/home/wfollowm/APRL_Workspace/build/cluster-tools/binding_probe 256 1 \
  > binding.log 2>&1
timeout --signal=TERM --kill-after=30s 1000s "$PYTHON_BIN" scripts/run.py run --root "$root"
