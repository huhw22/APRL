#!/usr/bin/env bash
#SBATCH --job-name=aprl-held-cut-small
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Pass a fresh frozen campaign directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 1 && ${SLURM_CPUS_PER_TASK:-0} == 128 ]] || exit 2
export LD_LIBRARY_PATH="$APRL_WORKSPACE/build/cowan-accepted-r2-20261025/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
cd "$root"
trap 'rc=$?; echo "$rc" > JOB_EXIT_CODE; date --iso-8601=seconds > JOB_ENDED_AT' EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
timeout --signal=TERM --kill-after=10 1500 "$PYTHON_BIN" "$root/run.py" run --root "$root"
