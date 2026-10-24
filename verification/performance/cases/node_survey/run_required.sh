#!/usr/bin/env bash
# Required regression on the explicitly authorized, unmanaged CPU test host.
set -Eeuo pipefail
[[ $# == 2 ]] || { echo 'Usage: run_required.sh BUILD_DIRECTORY NEW_RESULT_DIRECTORY' >&2; exit 2; }
build_dir=$(realpath "$1")
result_dir=$2
workspace=/public/home/wfollowm/APRL_Workspace
[[ $(hostname) == gnode1 && ${APRL_ALLOW_UNSCHEDULED_CPU_TESTS:-} == 1 ]] || {
  echo 'Requires gnode1 and explicit unscheduled CPU-test authorization' >&2; exit 2;
}
[[ -z ${SLURM_JOB_ID:-} && "$result_dir" == /* && ! -e "$result_dir" ]] || {
  echo 'Use a fresh absolute result directory; do not fabricate Slurm variables' >&2; exit 2;
}
[[ "$build_dir" == "$workspace/build/"* && -x "$build_dir/aprl" ]] || {
  echo 'Requires an independent verification build inside the workspace' >&2; exit 2;
}
grep -Fxq 'MPIEXEC_EXECUTABLE:FILEPATH=/usr/bin/mpiexec' "$build_dir/CMakeCache.txt"
grep -Fxq 'MPIEXEC_PREFLAGS:STRING=--host;gnode1:16;--map-by;slot:PE=2;--bind-to;core;--nooversubscribe' "$build_dir/CMakeCache.txt"
source "$workspace/cluster.env"
exec 9>>"$workspace/results/node_survey/.cpu-calibration.lock"
flock -n 9 || { echo 'Another CPU test holds the calibration lock' >&2; exit 2; }
required_resource_snapshot=$("$PYTHON_BIN" - <<'PY'
import json
import os
from pathlib import Path
import time
active = []
for entry in Path('/proc').iterdir():
    if not entry.name.isdigit():
        continue
    try:
        name = Path(os.fsdecode((entry / 'cmdline').read_bytes().split(b'\0')[0])).name
        if name in ('aprl', 'mpirun', 'mpiexec', 'orted'):
            active.append({'pid': int(entry.name), 'executable': name})
    except (OSError, ValueError):
        continue
def counters():
    values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
    return sum(values), values[3]
before, idle_before = counters()
time.sleep(2)
after, idle_after = counters()
fraction = (idle_after - idle_before) / (after - before) if after > before else 0
snapshot = {'loadavg': os.getloadavg(), 'cpu_idle_fraction': fraction,
            'sample_seconds': 2, 'compute_processes': active}
if active or fraction < .95:
    raise SystemExit('gnode1 has active work: ' + json.dumps(snapshot))
print(json.dumps(snapshot))
PY
)
mkdir "$result_dir"
mkdir "$build_dir/.gnode1-required-started"
finish() {
  code=$?
  trap - EXIT
  if [[ -d "$build_dir/Testing" ]]; then
    cp -a "$build_dir/Testing" "$result_dir/Testing" || code=1
  fi
  printf '%s\n' "$code" > "$result_dir/EXIT_CODE"
  if (( code == 0 )); then
    printf 'All registered required tests passed on gnode1\n' > "$result_dir/PASS"
  else
    printf 'Required regression failed with exit %s\n' "$code" > "$result_dir/FAIL"
  fi
  exit "$code"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
exec >"$result_dir/driver.log" 2>&1
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
date --iso-8601=seconds
hostname
lscpu --parse=CPU,CORE,SOCKET,NODE,ONLINE > "$result_dir/cpu-topology.csv"
sha256sum "$build_dir/aprl" "$build_dir/CMakeCache.txt" > "$result_dir/build.sha256"
cp "$build_dir/CMakeCache.txt" "$result_dir/CMakeCache.txt"
printf '%s\n' "$required_resource_snapshot" > "$result_dir/initial-resources.json"
if [[ -d "$build_dir/runtime-libs" ]]; then
  # Optional build-local copies of libraries absent from gnode1; no system
  # installation or shared production runtime is changed.
  export LD_LIBRARY_PATH="$build_dir/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  sha256sum "$build_dir/runtime-libs/"*.so* > "$result_dir/runtime-libs.sha256"
fi
# A loader failure must not accidentally satisfy a WILL_FAIL negative test.
"$PYTHON_BIN" - "$build_dir" "$result_dir/runtime-links.txt" <<'PY'
from pathlib import Path
import os
import subprocess
import sys
root = Path(sys.argv[1])
executables = [root / name for name in
    ('aprl', 'particle_text_to_hdf5', 'uniform_lab_plane_bunch_generator',
     'sase_bunch_generator')]
executables += list((root / 'verification').glob('aprl_*'))
executables += [root / 'postprocess' / name / name for name in
    ('field_reconstruction', 'field_plane_analysis', 'field_power_compare',
     'trajectory_radiation', 'energy_closure')]
failed = []
with open(sys.argv[2], 'w') as record:
    for executable in executables:
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise SystemExit('Missing required executable: ' + str(executable))
        result = subprocess.run(['ldd', '-r', str(executable)],
                                text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        record.write(str(executable) + '\n' + result.stdout + '\n')
        if result.returncode or 'not found' in result.stdout or 'undefined symbol' in result.stdout:
            failed.append(str(executable))
if failed:
    raise SystemExit('Runtime dependency preflight failed: ' + ', '.join(failed))
PY
OMP_NUM_THREADS=4 timeout --signal=TERM --kill-after=10 60 \
  mpirun --host gnode1:16 -np 4 --map-by slot:PE=4 --bind-to core \
  --nooversubscribe --report-bindings \
  "$workspace/build/node-survey/binding_probe" 4 4 > "$result_dir/binding.log" 2>&1
grep -q '^\[binding\] result=PASS failures=0' "$result_dir/binding.log"
# A new build owns its CTest output paths. Single-process tests may create up
# to four threads; MPI tests use the explicit CPU mapping recorded in its cache.
timeout --signal=TERM --kill-after=20 1200 taskset -c 0-15 \
  ctest --test-dir "$build_dir" -L required -j 1 --stop-on-failure --no-tests=error --output-on-failure
date --iso-8601=seconds
