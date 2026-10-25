#!/usr/bin/env bash
# Reuse a verified required suite without changing its build or result directory.
#SBATCH --job-name=aprl-cowan-acceptance
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Pass new prepared acceptance directory}")
build=$(realpath "${2:?Pass unchanged previously verified build directory}")
prior=$(realpath "${3:?Pass previous acceptance result directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
[[ "$build" == "$APRL_WORKSPACE/build/"* && -x "$build/aprl" && "$root" != "$prior" ]] || exit 2
[[ ! -e "$root/JOB_STARTED_AT" && ! -e "$root/ACCEPTANCE_PASS" ]] || exit 2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0 HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
ulimit -s unlimited
cd "$root"
mkdir .continuation-started
finish() {
  rc=$?
  trap - EXIT
  echo "$rc" > "$root/JOB_EXIT_CODE"
  date --iso-8601=seconds > "$root/JOB_ENDED_AT"
  if (( rc == 0 )); then
    echo 'Verified prior required suite and old race; new two-node fault and 256-rank scientific comparisons accepted.' > "$root/ACCEPTANCE_PASS"
  fi
  exit "$rc"
}
trap finish EXIT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
scontrol show hostnames "$SLURM_JOB_NODELIST" > allocation-nodes.txt
sha256sum --check frozen.sha256 > frozen-verification.log
sha256sum --check build-frozen.sha256 > build-verification.log
(cd "$prior" && sha256sum --check frozen.sha256) > prior-frozen-verification.log
"$PYTHON_BIN" - "$root" "$build" "$prior" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import yaml

root, build, prior = map(Path, sys.argv[1:])

def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()

def read_json(path):
    return json.loads(path.read_text())

plan, previous = (read_json(p / 'plan.json') for p in (root, prior))
assert plan == previous, 'Continuation must preserve the full frozen acceptance plan'
record, old_record = (read_json(p / 'build-record.json') for p in (root, prior))
fixed = previous['variants']['normal-fixed']
revision = fixed['source_revision']
assert re.fullmatch(r'[0-9a-f]{40}', revision), 'Exact compiled revision required'
assert record['compiled_revision'] == old_record['compiled_revision'] == revision
assert old_record['source_revision'] == revision
# The sole allowed changed prior source is the reviewed Slurm-status checker.
# Its original bytes remain recoverable from the exact compiled Git revision;
# binaries, caches and the original configured build driver must still match.
source = Path('/public/home/wfollowm/APRL_Workspace/source/aprl')
allowed_source = source / 'verification/check_duration_barrier_fault.cmake.in'
original_check = (prior / 'build-verification.log').read_text().splitlines()
frozen_entries = (prior / 'build-frozen.sha256').read_text().splitlines()
assert len(original_check) == len(frozen_entries) and all(line.endswith(': OK') for line in original_check)
verified = []
for line in frozen_entries:
    expected, filename = line.split('  ', 1)
    path = Path(filename)
    path = path if path.is_absolute() else prior / path
    assert filename + ': OK' in original_check, 'Prior startup did not verify frozen artifact'
    if digest(path) == expected:
        verified.append(filename + ': OK')
        continue
    assert path.resolve() == allowed_source, 'Prior build artifact changed: ' + filename
    original = subprocess.run(['git', '-C', str(source), 'show',
                               revision + ':verification/check_duration_barrier_fault.cmake.in'],
                              check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    assert hashlib.sha256(original).hexdigest() == expected, 'Original checker source identity changed'
    verified.append(filename + ': ORIGINAL SOURCE VERIFIED AT ' + revision)
(root / 'prior-build-verification.log').write_text('\n'.join(verified) + '\n')
assert Path(fixed['binary']).resolve() == build / 'aprl'
assert Path(old_record['binary']).resolve() == Path(record['binary']).resolve() == build / 'aprl'
assert record['binary_sha256'] == old_record['binary_sha256'] == fixed['binary_sha256']
for name in ('normal-control', 'normal-fixed', 'bounded-1', 'bounded-2'):
    assert plan['variants'][name] == previous['variants'][name], 'Changed program/termination contract: ' + name
    item = plan['variants'][name]
    assert digest(item['binary']) == item['binary_sha256'], 'Changed binary: ' + name
    assert not (root / name / 'aprl.log').exists() and not (root / name / 'EXIT_CODE').exists(), 'Fresh result required'
for name in ('control', 'fixed'):
    current, old = plan['toolchains'][name], previous['toolchains'][name]
    assert current == old and digest(current['cache']) == current['cache_sha256'], 'Changed CMake provenance'
assert plan['toolchains']['control']['selected'] == plan['toolchains']['fixed']['selected']
assert Path(plan['toolchains']['fixed']['cache']).resolve() == build / 'CMakeCache.txt'
assert record['production_sha256'] == old_record['production_sha256'] == plan['variants']['normal-control']['binary_sha256']
assert Path(record['production_binary']).resolve() == Path(plan['variants']['normal-control']['binary']).resolve()
flags = build / 'CMakeFiles/direct_eb_core.dir/flags.make'
assert re.search(r'APRL_SOURCE_REVISION=\\"' + re.escape(revision) + r'\\"', flags.read_text()), 'Compile flags revision differs'
manifest = yaml.safe_load((prior / 'two-node-fault/run-manifest.yaml').read_text())
assert manifest['source_revision'] == revision, 'Prior actual runtime revision differs'
assert (prior / 'required/PASS').read_text().strip() == 'All registered required checks passed'
assert (prior / 'required/EXIT_CODE').read_text().strip() == '0'
suite = ET.parse(prior / 'required/junit.xml').getroot()
assert suite.tag == 'testsuite' and int(suite.get('tests', '-1')) == 72
assert all(int(suite.get(k, '-1')) == 0 for k in ('failures', 'disabled', 'skipped'))
assert int(suite.get('errors', '0')) == 0
tests = suite.findall('testcase')
assert len(tests) == len({t.get('name') for t in tests}) == 72
assert all(t.get('name', '').startswith('required.') and t.get('status') == 'run' and
           all(t.find(k) is None for k in ('failure', 'error', 'skipped')) for t in tests)
assert {'required.runtime.duration-close', 'required.runtime.duration-barrier-watchdog'} <= {t.get('name') for t in tests}
assert '100% tests passed, 0 tests failed out of 72' in (prior / 'required/ctest.log').read_text()
reproduced = read_json(prior / 'old-race/REPRODUCED.json')
assert reproduced == {'expected_failure_reproduced': True, 'root_delayed_close_started': True,
                      'root_delayed_close_completed': False}
assert (prior / 'old-race/DRIVER_EXIT_CODE').read_text().strip() == '1'
old_log = (prior / 'old-race/run.log').read_text()
assert 'DURATION_CLOSE_DELAY_BEGIN' in old_log and 'mesh.duration was exhausted' in old_log
assert 'DURATION_CLOSE_DELAY_COMPLETE' not in old_log

copies = []
def copy_evidence(source, destination):
    assert not destination.exists(), 'Refusing to overwrite copied evidence'
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected = digest(source)
    shutil.copy2(source, destination)
    assert digest(destination) == expected == digest(source), 'Evidence changed during copy'
    copies.append({'source': str(source), 'copy': str(destination), 'sha256': expected})

for directory in ('required', 'old-race'):
    assert not (root / directory).exists(), 'Fresh evidence destination required'
    for source in sorted((prior / directory).rglob('*')):
        if source.is_file():
            copy_evidence(source, root / source.relative_to(prior))
for filename in ('plan.json', 'build-record.json', 'frozen.sha256', 'build-frozen.sha256',
                 'frozen-verification.log', 'build-verification.log',
                 'allocation.txt', 'allocation-nodes.txt', 'JOB_STARTED_AT', 'JOB_ENDED_AT', 'JOB_EXIT_CODE',
                 'two-node-fault/run-manifest.yaml'):
    copy_evidence(prior / filename, root / 'reused-evidence' / filename)
copy_evidence(flags, root / 'reused-evidence/compile-flags.make')
(root / 'REUSED_EVIDENCE.json').write_text(json.dumps({
    'pass': True, 'prior_root': str(prior), 'unchanged_build': str(build),
    'compiled_revision': revision, 'binary_sha256': fixed['binary_sha256'],
    'required_tests': 72, 'required_failures': 0, 'required_skipped': 0,
    'old_race_reproduced': True, 'copied_files': copies,
    'scope': 'Reuse only passed required suite and old-race reproduction; prior whole campaign remains failed.'
}, indent=2) + '\n')
print('Verified and copied 72 required checks and the original race reproduction.')
PY
srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=10 \
  "$PYTHON_BIN" "$root/scripts/check-runtime.py" "$build" "$root" \
  > runtime-links.log 2>&1
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30)
timeout --signal=TERM --kill-after=10 90 "${step[@]}" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 256 1 > binding-256.log 2>&1
"$PYTHON_BIN" - <<'PY'
import re
from pathlib import Path
text = Path('binding-256.log').read_text()
assert '[binding] result=PASS failures=0' in text
rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+).*?status=PASS', text)
assert len(rows) == 256 and {int(r[0]) for r in rows} == set(range(256))
assert {r[1] for r in rows} == {'0'} and len({r[2] for r in rows}) == 2
assert len({(r[2], r[3]) for r in rows}) == 256
assert all(sum(r[2] == host for r in rows) == 128 for host in {r[2] for r in rows})
PY
echo "Verify bounded exit with an absent peer $(date --iso-8601=seconds)"
cmake -DAPRL_EXECUTABLE="$build/aprl" \
  -DMPI_FAULT_LIBRARY="$build/verification/libaprl_duration_barrier_fault.so" \
  -DOUTPUT_CHECKER="$build/verification/aprl_output_checks" \
  -DDURATION_CARD="$root/scripts/duration-close.yaml" \
  -DMPI_FAULT_WORK="$root/two-node-fault" \
  -DMPI_FAULT_LAUNCHER="$root/scripts/mpi-launch-two-nodes.sh" \
  -DMPI_FAULT_RANKS=2 \
  -P "$root/scripts/check_duration_barrier_fault.cmake" \
  > two-node-fault-driver.log 2>&1
for variant in normal-control normal-fixed bounded-1 bounded-2; do
  binary="$build/aprl"
  if [[ $variant == normal-control ]]; then binary="$APRL_WORKSPACE/build/candidate/aprl"; fi
  extra=()
  if [[ $variant == bounded-* ]]; then
    extra=(env "LD_PRELOAD=$build/verification/libaprl_duration_close_delay.so" APRL_TEST_DELAY_H5_CLOSE=1)
  fi
  echo "256-rank $variant $(date --iso-8601=seconds)"
  set +e
  timeout --signal=TERM --kill-after=20 240 "${step[@]}" "${extra[@]}" "$binary" "$root/$variant/run.yaml" > "$root/$variant/aprl.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$root/$variant/EXIT_CODE"
  [[ $rc != 124 && $rc != 137 ]] || exit "$rc"
  "$PYTHON_BIN" "$root/scripts/cases.py" check --root "$root" --variant "$variant" > "$root/$variant/check.log" 2>&1
done
"$PYTHON_BIN" "$root/scripts/cases.py" check --root "$root" > scientific-comparison.log 2>&1
echo "Accepted continuation $(date --iso-8601=seconds)"
