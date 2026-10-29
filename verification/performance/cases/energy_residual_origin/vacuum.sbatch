#!/usr/bin/env bash
#SBATCH --job-name=aprl-vacuum-ledger
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:05:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?New output directory}")
probe=$(realpath "${2:?Frozen probe binary}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/STARTED" ]]
date --iso-8601=seconds > "$root/STARTED"
trap 'code=$?; printf "%s\n" "$code" > "$root/EXIT_CODE"' EXIT
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp
export LD_LIBRARY_PATH="$(dirname "$probe")/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
timeout --signal=TERM --kill-after=10 180 srun --mpi="${SRUN_MPI:-pmix_v3}" --exact --nodes=1 --ntasks=1 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=15 "$probe" "$root" > "$root/probe.jsonl" 2>&1
