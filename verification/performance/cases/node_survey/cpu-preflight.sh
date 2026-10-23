#!/usr/bin/env bash
# Standalone CPU binding diagnostics, only after explicit user approval of the
# gnode1 exception to this workspace's usual Slurm-only computation policy.
set -Eeuo pipefail

fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
mode="${1:---plan}"
[[ "$#" -le 1 && ( "$mode" == --plan || "$mode" == --run ) ]] || {
  printf 'Usage: %s [--plan|--run]\n' "$0" >&2
  exit 2
}
workspace="${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}"
[[ "$workspace" == /* ]] || fail 'APRL_WORKSPACE must be an absolute path.'
probe="$workspace/build/node-survey/binding_probe"
probe_source="$workspace/source/aprl/verification/performance/cases/node_survey/binding_probe.cpp"
manifest="$workspace/tests/node_survey/prepared-manifest.sha256"
load_limit=2
ranks_list=(16 8 4 2)
threads_list=(1 2 4 8)

print_plan() {
  printf '%s\n' \
    'Mode: plan only; no MPI launch, scientific simulation, SSH, or Slurm submission.' \
    'Run gate: explicit user approval, APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1, and hostname exactly gnode1.' \
    'Budget: 16 distinct physical CPU cores; SMT siblings are not additional workers.' \
    'Layouts: 16x1, 8x2, 4x4, 2x8, sequentially; binding diagnostics only.' \
    'Guard: workspace flock; reject one-minute load above 2 before each launch.' \
    'Timeout: 60 seconds per probe, then at most 10 seconds before killing this probe.' \
    'Topology: probe permits an affinity mask containing SMT siblings of one core, but rejects physical-core sharing between workers.'
  printf 'Frozen manifest: %s\nProbe: %s\nSource: %s\n' "$manifest" "$probe" "$probe_source"
  printf 'New results: %s/results/node_survey/cpu-preflight-gnode1-UTC-XXXXXX\n' "$workspace"
  printf '%s\n' 'Record SHA256, ldd -r, MPI/toolchain environment, topology, launch commands, binding logs, EXIT_CODE, and PASS only after all probes succeed.'
  local index
  for index in "${!ranks_list[@]}"; do
    printf 'OMP_NUM_THREADS=%s OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close ' "${threads_list[$index]}"
    printf 'timeout --kill-after=10 60 mpirun --host gnode1:16 -np %s --map-by slot:PE=%s --bind-to core --nooversubscribe --report-bindings %q %s %s\n' \
      "${ranks_list[$index]}" "${threads_list[$index]}" "$probe" "${ranks_list[$index]}" "${threads_list[$index]}"
  done
}

if [[ "$mode" == --plan ]]; then
  print_plan
  exit 0
fi

[[ "${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-}" == 1 ]] || \
  fail 'Execution requires explicit user approval of the no-Slurm exception and APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1.'
[[ "$(hostname)" == gnode1 ]] || fail 'This standalone preflight can run only on hostname exactly gnode1.'
[[ -f "$workspace/cluster.env" ]] || fail 'Prepared workspace cluster.env is missing.'
source "$workspace/cluster.env"
for utility in flock timeout mpirun ldd lscpu sha256sum awk; do
  command -v "$utility" >/dev/null || fail "Required utility is missing: $utility"
done
[[ -x "$probe" && -f "$probe_source" && -f "$manifest" ]] || \
  fail 'The standalone probe and frozen manifest must be prepared and reviewed first.'
[[ -x "$APRL_BASELINE" && -x "$APRL_CANDIDATE" && -x "$PYTHON_BIN" ]] || \
  fail 'Prepared baseline/candidate/Python paths are unavailable.'

mkdir -p "$workspace/results/node_survey"
exec 9>>"$workspace/results/node_survey/.cpu-preflight.lock"
flock -n 9 || fail 'Another CPU preflight holds this workspace lock; nothing will be stopped.'
run_dir="$(mktemp -d "$workspace/results/node_survey/cpu-preflight-gnode1-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")"
trap 'preflight_exit_code=$?; printf "%s\n" "$preflight_exit_code" > "$run_dir/EXIT_CODE"' EXIT
exec > >(tee "$run_dir/driver.log") 2>&1
printf 'Result directory: %s\n' "$run_dir"
printf 'node-survey-cpu-preflight\n' > "$run_dir/STAGE"
print_plan > "$run_dir/plan.txt"

check_load() {
  local load_one load_five load_fifteen other
  read -r load_one load_five load_fifteen other < /proc/loadavg
  printf '%s load1=%s load5=%s load15=%s limit=%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$load_one" "$load_five" "$load_fifteen" "$load_limit" >> "$run_dir/load-checks.log"
  awk -v observed="$load_one" -v maximum="$load_limit" 'BEGIN { exit !(observed <= maximum) }' || \
    fail "One-minute load $load_one exceeds $load_limit; refusing to start another probe."
}
check_load
cp "$manifest" "$run_dir/prepared-manifest.sha256"
sha256sum --check --strict --quiet "$run_dir/prepared-manifest.sha256" || fail 'Prepared artifacts changed.'
sha256sum "$probe" "$probe_source" "$APRL_BASELINE" "$APRL_CANDIDATE" \
  "$workspace/cluster.env" "${BASH_SOURCE[0]}" > "$run_dir/executables-and-sources.sha256"
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
{
  date -u --iso-8601=seconds
  hostname
  uname -a
  mpirun --version
  for variable in PATH LD_LIBRARY_PATH APRL_WORKSPACE APRL_BASELINE APRL_CANDIDATE PYTHON_BIN \
      APRL_ALLOW_UNSCHEDULED_CPU_TESTS OMP_DYNAMIC OMP_PLACES OMP_PROC_BIND OMP_THREAD_LIMIT \
      OPENBLAS_NUM_THREADS MKL_NUM_THREADS NUMEXPR_NUM_THREADS APRL_PROFILE_INTERVAL \
      SLURM_JOB_ID OMPI_MCA_btl OMPI_MCA_pml UCX_TLS UCX_NET_DEVICES; do
    printf '%s=%q\n' "$variable" "${!variable-}"
  done
} > "$run_dir/environment.log" 2>&1
lscpu > "$run_dir/lscpu.txt"
lscpu --parse=CPU,CORE,SOCKET,NODE,ONLINE > "$run_dir/cpu-topology.csv"
cp /proc/meminfo "$run_dir/meminfo.txt"
"$PYTHON_BIN" - "$run_dir/cpu-topology.csv" > "$run_dir/python-topology-check.json" <<'PY'
import csv, json, os, sys
import h5py, numpy, yaml
with open(sys.argv[1]) as stream:
    rows = list(csv.reader(line for line in stream if not line.startswith('#')))
online = [row for row in rows if row[4].strip().lower() in ('y', 'yes', '1')]
physical = {(row[2], row[1]) for row in online}
record = {'python': sys.executable, 'numpy': numpy.__version__, 'h5py': h5py.__version__,
          'hdf5': h5py.version.hdf5_version, 'yaml': yaml.__version__,
          'online_logical_cpus': len(online), 'online_physical_cores': len(physical),
          'shell_cpu_affinity': sorted(os.sched_getaffinity(0)),
          'required_distinct_physical_cores': 16}
print(json.dumps(record, indent=2))
if len(physical) != 16:
    raise SystemExit('gnode1 physical-core topology differs from the reviewed 16-core configuration')
PY
labels=(probe baseline candidate)
executables=("$probe" "$APRL_BASELINE" "$APRL_CANDIDATE")
for index in "${!labels[@]}"; do
  dependency_log="$run_dir/ldd-${labels[$index]}.log"
  ldd -r "${executables[$index]}" > "$dependency_log" 2>&1 || fail "Dependency audit failed: ${labels[$index]}"
  if grep -Eq 'not found|undefined symbol:' "$dependency_log"; then
    fail "Unresolved dynamic dependency: ${labels[$index]}"
  fi
done

for index in "${!ranks_list[@]}"; do
  ranks="${ranks_list[$index]}"
  threads="${threads_list[$index]}"
  export OMP_NUM_THREADS="$threads"
  check_load
  command_line=(timeout --kill-after=10 60 mpirun --host gnode1:16 -np "$ranks"
    --map-by "slot:PE=$threads" --bind-to core --nooversubscribe --report-bindings
    -x APRL_ALLOW_UNSCHEDULED_CPU_TESTS -x OMP_NUM_THREADS -x OMP_DYNAMIC -x OMP_PLACES -x OMP_PROC_BIND
    "$probe" "$ranks" "$threads")
  prefix="$run_dir/binding-r${ranks}-t${threads}"
  printf '%q ' "${command_line[@]}" > "$prefix.command"
  printf '\n' >> "$prefix.command"
  printf 'Checking %s MPI ranks x %s threads on 16 physical cores\n' "$ranks" "$threads"
  if "${command_line[@]}" > "$prefix.log" 2>&1; then
    probe_exit_code=0
  else
    probe_exit_code=$?
  fi
  printf '%s\n' "$probe_exit_code" > "$prefix.EXIT_CODE"
  [[ "$probe_exit_code" == 0 ]] || fail "Binding probe ${ranks}x${threads} failed with exit $probe_exit_code; see $prefix.log"
  grep -Fq '[binding] result=PASS failures=0' "$prefix.log" || fail 'Binding probe did not report a complete pass.'
done
sha256sum --check --strict --quiet "$run_dir/prepared-manifest.sha256" || fail 'Prepared artifacts changed during preflight.'
sha256sum --check --strict --quiet "$run_dir/executables-and-sources.sha256" || fail 'Recorded binaries or sources changed during preflight.'
printf 'All four standalone CPU binding layouts passed; no scientific simulation was run.\n' > "$run_dir/PASS"
printf 'Results: %s\n' "$run_dir"
