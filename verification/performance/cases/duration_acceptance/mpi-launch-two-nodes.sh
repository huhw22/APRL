#!/usr/bin/env bash
# Two-rank fault injection, one APRL process per reserved compute node.
set -Eeuo pipefail
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_JOB_NUM_NODES:-0} == 2 ]] || exit 2
[[ $# -ge 3 && $1 == -n && $2 == 2 ]] || exit 2
shift 2
[[ $1 == /* && -x $1 && ${OMP_NUM_THREADS:-1} == 1 ]] || exit 2
exec srun --mpi="${SRUN_MPI:-pmix_v3}" --exact --nodes=2 \
  --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 --cpu-bind=cores \
  --kill-on-bad-exit=0 --wait=10 "$@"
