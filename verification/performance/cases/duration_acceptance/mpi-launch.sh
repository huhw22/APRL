#!/usr/bin/env bash
# CTest steps inside the already reserved two-node acceptance allocation.
set -Eeuo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_JOB_NUM_NODES:-0} == 2 ]] || exit 2
[[ $# -ge 3 && $1 == -n && $2 =~ ^[1-9][0-9]*$ ]] || exit 2
ranks=$2
shift 2
threads=${OMP_NUM_THREADS:-1}
[[ $threads =~ ^[1-9][0-9]*$ && $1 == /* && -x $1 ]] || exit 2
(( ranks * threads <= 128 )) || exit 2
exec srun --mpi="${SRUN_MPI:-pmix_v3}" --exact --nodes=1 \
  --ntasks="$ranks" --ntasks-per-node="$ranks" --cpus-per-task="$threads" \
  --cpu-bind=cores --kill-on-bad-exit=1 --wait=10 "$@"
