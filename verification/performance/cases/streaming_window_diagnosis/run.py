#!/usr/bin/env python3
"""Freeze and run bounded residual-CPML comparisons on gnode1.

Successful execution establishes input/charge/finite-output gates. Gauss errors
are recorded, never used as acceptance thresholds. No production speed is claimed.
"""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import time


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
SOURCE = WORKSPACE / 'source/aprl'
CSV_NAMES = ('initial-particles.csv', 'release-particles.csv', 'detector-crossings.csv',
             'history.csv', 'field-plane.csv')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value, exclusive=False):
    with path.open('x' if exclusive else 'w') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def prepare(args):
    root, binary = args.root.resolve(), args.binary.resolve()
    binding = WORKSPACE / 'build/node-survey/binding_probe'
    checks = args.checks_binary.resolve() if args.checks_binary else None
    oracle_checks = args.oracle_checks_binary.resolve() if args.oracle_checks_binary else None
    if oracle_checks:require(oracle_checks.is_file(), 'Missing oracle checks')
    require(not root.exists(), 'Use a fresh result directory')
    for executable in ([binary, binding] + ([checks] if checks else [])):
        require(executable.is_file() and os.access(executable, os.X_OK),
                'Missing executable: ' + str(executable))
    cases = []
    specifications = args.case or ['b0-full-g24:0:24:full', 'b0-cpml-g24:0:24:cpml', 'b0-residual-g24:0:24:residual']
    for specification in specifications:
        parts = specification.split(':')
        require(len(parts) in (4,5,6), 'Use name:beta:guard:mode')
        name, beta, guard, mode = parts[:4]
        diagnostic_steps = int(parts[4]) if len(parts)>=5 else 0
        require(0<=diagnostic_steps<=176, 'Invalid diagnostic length')
        require(name and name.replace('-', '').replace('_', '').isalnum(), 'Invalid case name')
        require(mode in ('full', 'full-record', 'cpml', 'residual', 'residual-refresh', 'residual-oracle-feed', 'residual-oracle-cpml', 'residual-oracle-all', 'residual-native-b', 'held-record', 'residual-held-reference', 'residual-held-zero'), 'Unknown mode')
        require(float(beta) in (0., .23) and int(guard) in (16, 24, 40), 'Outside bounded fixture')
        require(name not in [case['name'] for case in cases], 'Duplicate case name')
        require(mode.startswith('residual') or mode=='held-record' or int(guard) == 24, 'Legacy replay uses only guard24')
        cases.append({'name': name, 'beta': beta, 'guard': int(guard), 'mode': mode,
                      'arguments': [beta, str(int(guard)), mode, str(diagnostic_steps)] + ([parts[5]] if len(parts)==6 else [])})
    previous_review = args.previous_review.resolve()
    previous = json.loads(previous_review.read_text())
    require(previous['format'] == 'aprl-streaming-cpml-review-v1' and
            previous['cohort_charge_and_input_consistency_status'] == 'PASS',
            'Use the independently reviewed original total-field CPML results')
    legacy_replay, legacy_files = {}, [previous_review]
    for case in cases:
        if case['mode'].startswith('residual') or case['mode']=='held-record':
            continue
        matches = [old for old in previous['cases'] if old['metrics']['beta_box'] == float(case['beta'])
                   and old['metrics']['mode'] == ('full' if case['mode']=='full-record' else case['mode']) and old['metrics']['guard_cells'] == case['guard']]
        require(len(matches) == 1, 'Missing unique original case for legacy replay: ' + case['name'])
        directory = Path(matches[0]['path']).resolve()
        legacy_replay[case['name']] = {'directory': str(directory),
                                      'csv_sha256': {name: digest(directory / name) for name in CSV_NAMES}}
        legacy_files.extend(directory / name for name in CSV_NAMES)
    snapshots = [SOURCE / 'verification' / name for name in
                 ('streaming_cpml_bunch_probe.cpp', 'streaming_oracle_checks.cpp', 'streaming_cpml_residual_checks.cpp', 'streaming_cpml_boundary.h',
                  'streaming_cpml_boundary.cpp', 'streaming_field_oracle.h', 'streaming_reference_boundary.h',
                  'streaming_reference_boundary.cpp', 'moving_window_remap.h',
                  'moving_window_remap.cpp', 'CMakeLists.txt')]
    snapshots += [SOURCE / 'CMakeLists.txt']
    for path in snapshots:
        require(path.is_file(), 'Missing source: ' + str(path))
    source_files = sorted(p for p in (SOURCE / 'src').rglob('*') if p.is_file())
    require(source_files, 'No core source files found')
    analysis_files = sorted(Path(__file__).resolve().parent.glob('*.py'))
    analysis_files += [Path(__file__).resolve().parent / 'README.md']
    source_files += snapshots + analysis_files
    revision = subprocess.check_output(
        ['git', '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip()
    # HEAD includes staged as well as unstaged tracked changes. Untracked probe
    # files are preserved explicitly above and represented in source-files.json.
    diff = subprocess.check_output(['git', '-C', str(SOURCE), 'diff', '--binary', 'HEAD'])
    status = subprocess.check_output(
        ['git', '-C', str(SOURCE), 'status', '--porcelain=v1'], text=True)
    root.mkdir(parents=True)
    shutil.copy2(__file__, root / 'run.py')
    (root / 'analysis').mkdir()
    for path in analysis_files:
        shutil.copy2(path, root / 'analysis' / path.name)
    for path in snapshots:
        destination = root / 'source' / path.relative_to(SOURCE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    save(root / 'source-files.json', {str(p): digest(p) for p in source_files})
    (root / 'source.diff').write_bytes(diff)
    (root / 'source-status.txt').write_text(status)
    runtime_files = []
    for executable in [binary] + ([checks] if checks else []):
        for directory in (executable.parent / 'runtime-libs', executable.parent.parent / 'runtime-libs'):
            if directory.is_dir():
                runtime_files.extend(p for p in directory.rglob('*') if p.is_file())
    plan = {'format': 'aprl-streaming-window-diagnosis-campaign-v1', 'cases': cases,
            'binary': str(binary), 'binary_sha256': digest(binary),
            'checks_binary': str(checks) if checks else None,
            'checks_binary_sha256': digest(checks) if checks else None,
            'oracle_checks_binary': str(oracle_checks) if oracle_checks else None,
            'oracle_checks_sha256': digest(oracle_checks) if oracle_checks else None,
            'previous_review': str(previous_review), 'previous_review_sha256': digest(previous_review),
            'legacy_binary_sha256': previous['binary_sha256'], 'legacy_replay': legacy_replay,
            'binding_probe': str(binding), 'source_base': revision,
            'ranks': 1, 'threads': 1, 'case_timeout_s': 180,
            'total_timeout_s': 7200,
            'scope': 'Bounded verification-only field-interface diagnosis; unchanged src kernels and original geometry. '
                     'Explicit modes select frozen, native-B reconstructed, physical-history oracle, or independent prescribed-source references. '
                     'Zero-source controls are not physical PIC runs; full held reference cost is included and not a production speedup. '
                     'Legacy full/cpml modes must reproduce all scientific CSVs byte for byte'}
    save(root / 'plan.json', plan)
    frozen = [p for p in root.rglob('*') if p.is_file()]
    frozen += ([oracle_checks] if oracle_checks else []) + [binary, binding] + runtime_files + legacy_files + ([checks] if checks else [])
    save(root / 'frozen.json', {str(p): digest(p) for p in sorted(set(frozen))})
    print(root / 'plan.json')


def finite_metrics(value):
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(finite_metrics(v) for v in value.values())
    if isinstance(value, list):
        return all(finite_metrics(v) for v in value)
    return True


def idle_check():
    other = []
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            raw = (entry / 'cmdline').read_bytes().split(b'\0')[0]
            name = Path(os.fsdecode(raw)).name
            if name in ('aprl', 'mpirun', 'mpiexec', 'orted') or name.startswith('aprl_'):
                other.append({'pid': int(entry.name), 'executable': name})
        except (OSError, ValueError):
            continue
    require(not other, 'Existing compute process: ' + json.dumps(other))

    def counters():
        values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(values), values[3]

    before, idle_before = counters()
    time.sleep(2)
    after, idle_after = counters()
    require(after > before, 'Invalid CPU-idle sample')
    fraction = (idle_after - idle_before) / (after - before)
    require(fraction >= .95, 'Test host is busy: idle fraction ' + str(fraction))
    return fraction


def terminate_group(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


def launch(command, directory, env, timeout, log_path, record_path):
    started = time.monotonic()
    watchdog = shutil.which('timeout')
    require(watchdog is not None, 'Bounded-runtime watchdog is unavailable')
    bounded_command = [watchdog, '--signal=TERM', '--kill-after=5', str(timeout)] + command
    record = {'command': bounded_command, 'timeout_s': timeout, 'unix_time': time.time()}
    process = None
    try:
        with log_path.open('x') as log:
            process = subprocess.Popen(bounded_command, cwd=directory, env=env,
                                       stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            record['pid'] = process.pid
            try:
                code = process.wait(timeout=timeout + 10)
            except subprocess.TimeoutExpired:
                record['timed_out'] = True
                raise RuntimeError('Process exceeded bounded runtime: ' + str(directory))
            record['exit_code'] = code
            if code in (124, 137):
                record['timed_out'] = True
            require(code == 0, 'Probe exited ' + str(code) + ': ' + str(log_path))
    except BaseException as error:
        record['error'] = str(error) or type(error).__name__
        raise
    finally:
        if process is not None:
            terminate_group(process)
            record['exit_code'] = process.returncode
        record['seconds'] = time.monotonic() - started
        save(record_path, record)


def run(args):
    root = args.root.resolve()
    require(socket.gethostname().split('.')[0] == 'gnode1' and
            os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1' and
            not os.environ.get('SLURM_JOB_ID'),
            'Requires authorized gnode1 CPU session')
    require(not (root / 'STARTED.json').exists(), 'Campaign already started; use a fresh root')
    for name, sha in json.loads((root / 'frozen.json').read_text()).items():
        require(digest(Path(name)) == sha, 'Frozen input changed: ' + name)
    plan = json.loads((root / 'plan.json').read_text())
    require(plan['ranks'] == plan['threads'] == 1 and plan['case_timeout_s'] == 180
            and plan['total_timeout_s'] == 7200, 'Outside the bounded CPU-only plan')
    with (WORKSPACE / 'results/node_survey/.cpu-calibration.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        start_record = {'host': socket.gethostname(), 'unix_time': time.time(),
                        'phase': 'preflight'}
        save(root / 'STARTED.json', start_record, exclusive=True)
        args.claimed_run = True
        deadline = time.monotonic() + plan['total_timeout_s']
        start_record.update(cpu_idle_fraction=idle_check(), loadavg=os.getloadavg())
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
                   OMP_PROC_BIND='close', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   NUMEXPR_NUM_THREADS='1', APRL_PROFILE_INTERVAL='0')
        executable = shutil.which('mpiexec')
        require(executable is not None, 'MPI launcher is unavailable')
        launcher = [executable, '-n', '1', '--map-by', 'slot:PE=1',
                    '--bind-to', 'core', '--nooversubscribe']
        launch(launcher + [plan['binding_probe'], '1', '1'], root, env, 45,
               root / 'binding.log', root / 'binding-launch.json')
        require('[binding] result=PASS failures=0' in (root / 'binding.log').read_text(),
                'Physical-core binding check failed')
        if plan.get('checks_binary'):
            checks_directory = root / 'checks'
            checks_directory.mkdir()
            launch(launcher + [plan['checks_binary'], str(checks_directory)], root, env, plan['case_timeout_s'],
                   root / 'check.log', root / 'check-launch.json')
            require('Residual CPML checks PASS' in (root / 'check.log').read_text(),
                    'Independent residual-CPML oracle failed')
            check_metrics = json.loads((checks_directory / 'metrics.json').read_text())
            require(check_metrics.get('completed') is True and check_metrics.get('bitwise_failures') == 0
                    and finite_metrics(check_metrics), 'Residual oracle numerical gates failed')
        if plan.get('oracle_checks_binary'):
            checks_directory=root/'oracle-checks'
            checks_directory.mkdir()
            launch(launcher+[plan['oracle_checks_binary'],str(checks_directory)],root,env,plan['case_timeout_s'],root/'oracle-check.log',root/'oracle-check-launch.json')
            require('Native field oracle checks PASS' in (root/'oracle-check.log').read_text(), 'Native replay checks failed')
            checked=json.loads((checks_directory/'metrics.json').read_text())
            require(checked.get('completed') is True and checked.get('bitwise_failures')==0 and finite_metrics(checked), 'Invalid native oracle checks')
        start_record.update(phase='cases', launcher=launcher)
        save(root / 'STARTED.json', start_record)
        results = []
        for case in plan['cases']:
            require(time.monotonic() + plan['case_timeout_s'] + 10 < deadline,
                    'Campaign deadline: no time for another bounded case')
            directory = root / case['name']
            directory.mkdir()
            command = launcher + [plan['binary'], str(directory)] + case['arguments']
            launch(command, directory, env, plan['case_timeout_s'],
                   directory / 'probe.log', directory / 'launch.json')
            metrics_path = directory / 'metrics.json'
            metrics = json.loads(metrics_path.read_text())
            require(isinstance(metrics, dict) and finite_metrics(metrics),
                    'Invalid or nonfinite metrics: ' + case['name'])
            require(metrics.get('completed') is True, 'Incomplete probe: ' + case['name'])
            require(metrics.get('mode') == case['mode'] and
                    metrics.get('beta_box') == float(case['beta']) and
                    metrics.get('guard_cells') == case['guard'], 'Requested fixture was not executed')
            require(metrics.get('cpml_input') == ('residual' if case['mode'].startswith('residual') else 'total'),
                    'Unexpected CPML derivative input')
            for key in ('residual_reference_work_bytes', 'residual_reference_input_values',
                        'residual_reference_input_wall_s'):
                require(metrics.get(key, -1) >= 0, 'Missing residual-reference accounting: ' + key)
                if not case['mode'].startswith('residual'):
                    require(metrics[key] == 0, 'Original mode unexpectedly performs residual work')
            if case['mode'].startswith('residual'):
                require(metrics['residual_reference_work_bytes'] > 0 and
                        metrics['residual_reference_input_values'] > 0,
                        'Residual reference subtraction was not exercised')
            require(metrics.get('history_policy') == 'world-position' and
                    metrics.get('cpml_cells') == [0, 0, 20], 'Unexpected CPML configuration')
            require(metrics.get('particles') == 1024 and metrics.get('admitted') == 1024 and
                    metrics.get('steps') == 176 and metrics.get('seed') == 710271,
                    'Incomplete or changed particle cohort')
            require(metrics.get('full_nz') == 424 and metrics.get('physical_full_nz') == 384 and
                    metrics.get('window_nz') == (424 if (case['mode'].startswith('full') or case['mode']=='held-record') else 56 + 2 * case['guard']),
                    'Window size does not include the separate CPML layers')
            require(metrics.get('inner_guard_cells') == case['guard'] and
                    metrics.get('cpml_z_cells') == 20 and
                    metrics.get('physical_window_nz') == metrics['window_nz'] - 40,
                    'Nonabsorbing guard and CPML thickness were conflated')
            require(metrics.get('pic_pushes_in_cpml') == 0 and
                    metrics.get('min_pic_inner_clearance_m', -1) >= 0,
                    'Particle PIC advance reached the absorbing layer')
            peak = metrics.get('cpml_history_peak_abs', -1)
            require(peak >= 0 and bool(metrics.get('cpml_history_nonzero_seen')) == (peak > 0),
                    'Inconsistent CPML history audit')
            if case['mode'] != 'full':
                require(metrics['cpml_history_nonzero_seen'], 'Candidate did not exercise CPML')
            for key in ('min_frame_release_inner_m', 'min_frame_detector_inner_m',
                        'min_event_release_inner_m', 'min_event_detector_inner_m'):
                require(metrics.get(key, -1) >= 0, 'A release/detector plane entered CPML: ' + key)
            charge = abs(metrics['charge_C'])
            require(charge > 0 and metrics['max_continuity_relative'] <= 1e-10,
                    'Interior continuity gate failed')
            require(metrics['max_open_charge_balance_C'] / charge <= 1e-10,
                    'Open-volume charge gate failed')
            require(metrics['max_remap_charge_balance_C'] / charge <= 1e-10,
                    'Moved-grid charge gate failed')
            # Initial and propagated Gauss residuals are measurements only.
            for name in ('initial-particles.csv', 'release-particles.csv', 'detector-crossings.csv'):
                with (directory / name).open() as stream:
                    require(sum(1 for _ in stream) == 1025, 'Truncated cohort CSV: ' + name)
            with (directory / 'history.csv').open() as stream:
                require(sum(1 for _ in stream) == 178, 'Truncated step history')
            with (directory / 'field-plane.csv').open() as stream:
                require(sum(1 for _ in stream) == metrics['field_plane_samples'] + 1,
                        'Truncated field-plane CSV')
            csv_hashes = {name: digest(directory / name) for name in CSV_NAMES}
            if case['name'] in plan['legacy_replay']:
                require(csv_hashes == plan['legacy_replay'][case['name']]['csv_sha256'],
                        'Original mode scientific CSV replay differs: ' + case['name'])
            results.append({'case': case['name'], 'metrics_sha256': digest(metrics_path),
                            'csv_sha256': csv_hashes,
                            'legacy_csv_byte_identical': True if case['name'] in plan['legacy_replay'] else None})
        save(root / 'COMPLETED.json', {'unix_time': time.time(), 'cases': results,
             'cohort_and_charge_gates': 'COMPLETED',
             'physical_acceptance': 'Not assessed: no complete moving-boundary or production acceptance',
             'performance_scope': 'Single-run different-domain workload observations, not production implementation speedup'})
    print(root / 'COMPLETED.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--binary', type=Path, required=True)
    prep.add_argument('--previous-review', type=Path, required=True)
    prep.add_argument('--checks-binary', type=Path,
                      help='Optional independently compiled residual oracle; provided for the planned comparison')
    prep.add_argument('--oracle-checks-binary', type=Path)
    prep.add_argument('--case', action='append',
                      help='name:beta:guard:full|cpml|residual; default is the three static guard24 preflight cases')
    run_parser = sub.add_parser('run')
    run_parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    args.claimed_run = False

    def stopped(signum, _frame):
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, stopped)
    try:
        globals()[args.action](args)
    except BaseException as error:
        # Never overwrite another invocation's active/successful result when
        # its shared lock or existing STARTED marker correctly rejects us.
        if args.action == 'run' and args.claimed_run:
            save(args.root.resolve() / 'FAILED.json',
                 {'message': str(error) or type(error).__name__, 'unix_time': time.time()})
        raise


if __name__ == '__main__':
    main()
