#!/usr/bin/env bash
# CTest MPI launcher. This script only runs inside an existing Slurm job.
set -Eeuo pipefail
fail() { printf 'MPI launcher: %s\n' "$*" >&2; exit 2; }
[[ -n "${SLURM_JOB_ID:-}" ]] || fail 'a Slurm allocation is required'
[[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}}" == 1 ]] || fail 'this acceptance launcher requires one node'
[[ $# -ge 3 && "$1" == -n ]] || fail 'usage: mpi-launch.sh -n RANKS EXECUTABLE [ARG ...]'
ranks="$2"
shift 2
threads="${OMP_NUM_THREADS:-1}"
[[ "$ranks" =~ ^[1-9][0-9]*$ ]] || fail 'RANKS must be a positive integer'
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || fail 'OMP_NUM_THREADS must be a positive integer'
[[ "${SLURM_NTASKS:-}" =~ ^[1-9][0-9]*$ ]] || fail 'SLURM_NTASKS is missing or invalid'
[[ "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]] || fail 'SLURM_CPUS_PER_TASK is missing or invalid'
(( ranks <= SLURM_NTASKS )) || fail 'requested ranks exceed the allocation'
(( threads <= SLURM_CPUS_PER_TASK )) || fail 'requested threads exceed allocated CPUs per task'
[[ "$1" == /* && -x "$1" ]] || fail 'the executable must be an absolute, executable path'
exec srun --mpi="${SRUN_MPI:-pmix_v3}" --nodes=1 --ntasks="$ranks" \
  --cpus-per-task="$threads" --exact --cpu-bind=verbose,cores \
  --kill-on-bad-exit=1 "$@"
