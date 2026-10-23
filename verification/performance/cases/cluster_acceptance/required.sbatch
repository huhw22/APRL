#!/usr/bin/env bash
#SBATCH --job-name=aprl-required
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=aprl-required-%j.out
#SBATCH --error=aprl-required-%j.err
set -Eeuo pipefail
workspace="${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}"
source "$workspace/tests/cluster_acceptance/common.sh"
init_run required
: "${APRL_SMOKE_RESULT:?Set APRL_SMOKE_RESULT to the successful smoke result directory}"
require_pass "$APRL_SMOKE_RESULT" smoke
[[ "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}}" == 1 ]] || { echo 'Required regression needs one node' >&2; exit 2; }
[[ "${SLURM_NTASKS:-0}" -ge 4 && "${SLURM_CPUS_PER_TASK:-0}" -ge 2 ]] || { echo 'Required regression needs 4 tasks and 2 CPUs per task' >&2; exit 2; }
: "${APRL_VERIFICATION_BUILD_DIR:?Missing verification build directory}"
verification="$APRL_VERIFICATION_BUILD_DIR/verification"
held="$verification/work/hybrid-held"
card="$verification/hybrid-held.yaml"
launcher="$workspace/tests/cluster_acceptance/mpi-launch.sh"
comparator="$verification/aprl_hybrid_output_checks"
[[ -x "$comparator" && -f "$card" ]] || { echo 'Verification build is incomplete' >&2; exit 2; }

# CTest cards contain absolute scratch paths. Reusing this build could overwrite
# evidence, so every verification build may be executed only once.
if ! mkdir "$APRL_VERIFICATION_BUILD_DIR/.cluster-run-started"; then
  echo "This verification build was already used. Prepare a fresh build before retrying: $APRL_VERIFICATION_BUILD_DIR" >&2
  exit 2
fi
printf '%s\n' "$run_dir" > "$APRL_VERIFICATION_BUILD_DIR/.cluster-run-started/result-directory"
printf '%s\n' "$APRL_SMOKE_RESULT" > "$run_dir/SMOKE_RESULT"
snapshot_and_finish() {
  status=$?
  trap - EXIT
  set +e
  mkdir -p "$run_dir/verification"
  snapshot_status=$?
  for item in "$verification/work" "$APRL_VERIFICATION_BUILD_DIR/Testing" \
      "$APRL_VERIFICATION_BUILD_DIR/CMakeCache.txt" "$card"; do
    if [[ -e "$item" ]]; then
      cp -a -- "$item" "$run_dir/verification/" || snapshot_status=1
    fi
  done
  if (( status == 0 && snapshot_status != 0 )); then status=1; fi
  printf '%s\n' "$status" > "$run_dir/EXIT_CODE"
  if (( status == 0 )); then
    printf 'required and baseline-held checks passed\n' > "$run_dir/PASS"
  else
    printf 'required checks or evidence snapshot failed (exit %s)\n' "$status" > "$run_dir/FAIL"
  fi
  printf 'Regression results: %s\n' "$run_dir"
  exit "$status"
}
trap snapshot_and_finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0
# Slurm task/affinity does not guarantee that the batch shell has exactly the
# allocated CPUs. Discover them through a completed 4-rank, 2-thread step.
OMP_NUM_THREADS=2 srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=4 \
  --ntasks-per-node=4 --cpus-per-task=2 --cpu-bind=verbose,cores \
  --kill-on-bad-exit=1 "$workspace/build/cluster-tools/binding_probe" 4 2 \
  2>&1 | tee "$run_dir/binding-r4-t2.log"
ctest_cpus="$("$PYTHON_BIN" - "$run_dir/binding-r4-t2.log" <<'PY'
from pathlib import Path
import re
import sys

text = Path(sys.argv[1]).read_text(errors='replace')
if '[binding] result=PASS failures=0' not in text:
    raise SystemExit('Binding probe did not report a successful result')
pattern = re.compile(r'^\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=([0-9,]+) status=PASS reason=', re.MULTILINE)
records = pattern.findall(text)
if len(records) != 8 or {(int(r[0]), int(r[1])) for r in records} != {(rank, thread) for rank in range(4) for thread in range(2)}:
    raise SystemExit('Expected exactly eight successful 4-rank, 2-thread binding records')
if len({r[2] for r in records}) != 1 or len({(r[2], r[4], r[5]) for r in records}) != 8:
    raise SystemExit('Binding records must identify eight distinct physical cores on one node')
cpus = {int(cpu) for row in records for cpu in row[6].split(',')}
if len(cpus) != 8:
    raise SystemExit('Expected eight CPU IDs for this cluster with one hardware thread per core')
print(','.join(map(str, sorted(cpus))))
PY
)"
printf '%s\n' "$ctest_cpus" > "$run_dir/ctest-cpus.txt"
printf 'CTest CPU affinity, derived from the successful Slurm step: %s\n' "$ctest_cpus"
# The field check requests four threads. Restrict the CTest process and its
# direct children to the allocated eight CPUs without holding an enclosing
# Slurm step, so nested MPI test steps remain free to acquire their resources.
taskset --cpu-list "$ctest_cpus" ctest --test-dir "$APRL_VERIFICATION_BUILD_DIR" -L required -VV \
  --output-on-failure --output-junit "$run_dir/required-junit.xml" -j 1 \
  2>&1 | tee "$run_dir/ctest.log"

# The CTest held case compares candidate 1-thread and 2-thread results. Extend
# it with a same-rank baseline run using exactly the same card and input.
mkdir "$run_dir/held-baseline-r4-t1"
sha256sum "$card" "$held/particles.h5" > "$run_dir/held-input.sha256"
(
  cd "$run_dir/held-baseline-r4-t1"
  "$launcher" -n 4 "$APRL_BASELINE" "$card"
) > "$run_dir/held-baseline-r4-t1/run.log" 2>&1
sha256sum -c "$run_dir/held-input.sha256" | tee "$run_dir/held-input-check.log"
for threads in 1 2; do
  "$comparator" "$run_dir/held-baseline-r4-t1" "$held/threads-$threads" \
    2>&1 | tee "$run_dir/held-baseline-vs-threads-$threads.log"
done
"$PYTHON_BIN" - "$run_dir/held-baseline-r4-t1" "$held/threads-1" "$held/threads-2" \
  <<'PY' | tee "$run_dir/held-endpoint-check.log"
import math
from pathlib import Path
import re
import sys
import h5py
import numpy as np

for item in sys.argv[1:]:
    directory = Path(item)
    log = (directory / 'run.log').read_text(errors='replace')
    match = re.search(r'Configured stop reached: boost reference centre reached lab z=([0-9.eE+-]+) m \(target ([0-9.eE+-]+) m\)', log)
    if not match:
        raise RuntimeError(f'{directory}: missing physical endpoint marker')
    final_z, target_z = map(float, match.groups())
    if not all(map(math.isfinite, (final_z, target_z))) or target_z != 40e-6 or final_z < target_z:
        raise RuntimeError(f'{directory}: unexpected physical endpoint {match.group(0)}')
    if 'duration was exhausted' in log or 'stopped after a complete field step' in log:
        raise RuntimeError(f'{directory}: interrupted before the physical endpoint')
    ids = None
    for plane in ('activation', 'downstream'):
        with h5py.File(directory / 'detectors' / f'{plane}.h5', 'r') as file:
            group = file['particle_plane']
            records = group['records'][:]
            current_ids = np.sort(records['source_id'])
            if int(group['complete'][()]) != 1 or int(group['committed_records'][()]) != 10410 or len(current_ids) != 10410 or len(np.unique(current_ids)) != 10410:
                raise RuntimeError(f'{directory}/{plane}: each particle must be committed exactly once')
            if ids is not None and not np.array_equal(ids, current_ids):
                raise RuntimeError(f'{directory}: activation/downstream particle identities differ')
            ids = current_ids
    with h5py.File(directory / 'ledger' / 'energy-ledger.h5', 'r') as file:
        group = file['energy_ledger']
        records = group['records'][:]
        if int(group['complete'][()]) != 1 or int(group['committed_records'][()]) != len(records) or not len(records):
            raise RuntimeError(f'{directory}: incomplete energy ledger')
        if np.any(records['removed_macroparticles'] != 0) or np.any(records['active_macroparticles'] != 10410):
            raise RuntimeError(f'{directory}: particle conservation failed')
    print(f'{directory}: physical_endpoint_pass=True unique_particle_endpoints=10410 particle_conservation_pass=True')
PY
