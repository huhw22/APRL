#!/usr/bin/env bash
# Capture read-only process/runtime evidence while teardown stragglers are
# still alive.  This is started only for APRL_TEARDOWN_TRACE=all runs.
set -u

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 APRL_LOG OUTPUT NODES DELAY_SECONDS MAX_WAIT_SECONDS ACTIVE_SENTINEL" >&2
  exit 2
fi

aprl_log=$1
output=$2
nodes=$3
delay_seconds=$4
max_wait_seconds=$5
active_sentinel=$6
deadline=$((SECONDS + max_wait_seconds))
marker_seen=0

while ((SECONDS < deadline)); do
  if [[ -f "${aprl_log}" ]] &&
     grep -q 'phase=final-stop-sync-enter' "${aprl_log}" 2>/dev/null; then
    marker_seen=1
    break
  fi
  sleep 5
done

temporary="${output}.tmp-$$"
if [[ "${marker_seen}" -ne 1 ]]; then
  {
    echo "snapshot_status=no-solve-tail-marker-before-deadline"
    echo "time=$(date --iso-8601=seconds)"
  } > "${temporary}"
  mv -- "${temporary}" "${output}"
  exit 0
fi

sleep "${delay_seconds}"
primary_step_observed=0
if [[ -e "${active_sentinel}" ]]; then
  primary_step_observed=1
fi
mpi_finalize_complete_rank_count=$(
  awk '
    /\[teardown\]/ && /phase=mpi-finalize-complete/ {
      for (field = 1; field <= NF; ++field) {
        if ($field ~ /^rank=/) {
          split($field, parts, "=")
          seen[parts[2]] = 1
        }
      }
    }
    END {
      count = 0
      for (rank in seen) ++count
      print count
    }
  ' "${aprl_log}" 2>/dev/null
)
{
  echo "snapshot_status=attempted"
  echo "time=$(date --iso-8601=seconds)"
  echo "delay_after_first_solve_tail_marker_s=${delay_seconds}"
  echo "primary_step_observed_active=${primary_step_observed}"
  echo "mpi_finalize_complete_rank_count=${mpi_finalize_complete_rank_count:-0}"
  timeout 120 srun --overlap \
    --nodes="${nodes}" --ntasks="${nodes}" --ntasks-per-node=1 \
    --cpus-per-task=1 --mpi=none \
    bash -lc '
      echo "===== host=$(hostname) ====="
      mapfile -t pids < <(pgrep -x aprl 2>/dev/null || true)
      echo "aprl_pid_count=${#pids[@]}"
      for pid in "${pids[@]}"; do
        echo "--- pid=${pid} ---"
        readlink -f "/proc/${pid}/exe" 2>&1 || true
        tr "\0" " " < "/proc/${pid}/cmdline" 2>/dev/null || true
        echo
        ps -p "${pid}" -o pid,ppid,stat,etime,time,%cpu,%mem,nlwp,wchan:32,cmd
        ps -L -p "${pid}" -o pid,tid,psr,stat,%cpu,wchan:32,comm
        grep -E "^(Name|State|Threads|VmRSS|VmPeak|voluntary_ctxt_switches|nonvoluntary_ctxt_switches):" \
          "/proc/${pid}/status" 2>/dev/null || true
        echo "mapped MPI/PMIx/HDF5 libraries:"
        grep -Ei "/(libmpi|libopen-rte|libopen-pal|libpmix|libpmi|libhdf5|mca_)" \
          "/proc/${pid}/maps" 2>/dev/null | awk "!seen[\$0]++" || true
      done
    '
  echo "snapshot_command_status=$?"
} > "${temporary}" 2>&1
mv -- "${temporary}" "${output}"
