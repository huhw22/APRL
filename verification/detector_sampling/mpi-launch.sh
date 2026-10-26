#!/usr/bin/env bash
# Correctness tests only: a real Slurm allocation or the authorized small host.
set -Eeuo pipefail
[[ $# -ge 3 && $1 == -n && $2 =~ ^[1-9][0-9]*$ ]] || exit 2
ranks=$2
shift 2
threads=${OMP_NUM_THREADS:-1}
[[ $threads =~ ^[1-9][0-9]*$ && $1 == /* && -x $1 ]] || exit 2
if [[ -n ${SLURM_JOB_ID:-} ]]; then
  (( ranks * threads <= 128 )) || exit 2
  exec srun --mpi="${SRUN_MPI:-pmix_v3}" --exact --nodes=1 \
    --ntasks="$ranks" --ntasks-per-node="$ranks" --cpus-per-task="$threads" \
    --cpu-bind=cores --kill-on-bad-exit=1 --wait=10 "$@"
fi
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || exit 2
(( ranks * threads <= 16 )) || exit 2
exec mpiexec -n "$ranks" --map-by "slot:PE=$threads" --bind-to core \
  --nooversubscribe "$@"
