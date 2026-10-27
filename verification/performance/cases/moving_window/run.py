#!/usr/bin/env python3
"""Freeze and run bounded moving-window state/cost probes on gnode1.

Successful execution establishes only the probe's state-copy gates. It does
not accept a moving boundary's physics or establish production acceleration.
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
DEFAULT_CASE_NAMES = ('verify', 'bench64', 'bench128', 'bench256')
WAVE_CASE_NAMES = ('wave', 'wave-static', 'wave20', 'wave-static20')
CASE_NAMES = DEFAULT_CASE_NAMES + WAVE_CASE_NAMES


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
    require(not root.exists(), 'Use a fresh result directory')
    for executable in (binary, binding):
        require(executable.is_file() and os.access(executable, os.X_OK),
                'Missing executable: ' + str(executable))
    names = args.case or list(DEFAULT_CASE_NAMES)
    require(len(names) == len(set(names)), 'Duplicate case')
    cases = [{'name': name, 'arguments': [name] if name == 'verify' or name in WAVE_CASE_NAMES
              else ['bench', name[len('bench'):]]} for name in names]
    snapshots = [SOURCE / 'verification' / name for name in
                 ('moving_window_probe.cpp', 'moving_window_remap.h',
                  'moving_window_remap.cpp', 'moving_window_wave.h',
                  'moving_window_wave.cpp', 'CMakeLists.txt')]
    snapshots += [SOURCE / 'CMakeLists.txt']
    for path in snapshots:
        require(path.is_file(), 'Missing source: ' + str(path))
    source_files = sorted(p for p in (SOURCE / 'src').rglob('*') if p.is_file())
    require(source_files, 'No core source files found')
    source_files += snapshots + [Path(__file__).resolve()]
    revision = subprocess.check_output(
        ['git', '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip()
    # HEAD includes staged as well as unstaged tracked changes. Untracked probe
    # files are preserved explicitly above and represented in source-files.json.
    diff = subprocess.check_output(['git', '-C', str(SOURCE), 'diff', '--binary', 'HEAD'])
    status = subprocess.check_output(
        ['git', '-C', str(SOURCE), 'status', '--porcelain=v1'], text=True)
    root.mkdir(parents=True)
    shutil.copy2(__file__, root / 'run.py')
    for path in snapshots:
        destination = root / 'source' / path.relative_to(SOURCE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    save(root / 'source-files.json', {str(p): digest(p) for p in source_files})
    (root / 'source.diff').write_bytes(diff)
    (root / 'source-status.txt').write_text(status)
    runtime_files = []
    for directory in (binary.parent / 'runtime-libs', binary.parent.parent / 'runtime-libs'):
        if directory.is_dir():
            runtime_files.extend(p for p in directory.rglob('*') if p.is_file())
    plan = {'format': 'aprl-moving-window-state-cost-v1', 'cases': cases,
            'binary': str(binary), 'binary_sha256': digest(binary),
            'binding_probe': str(binding), 'source_base': revision,
            'ranks': 1, 'threads': 1, 'case_timeout_s': 180,
            'total_timeout_s': 7200,
            'scope': 'Same-width state remap/cost probes and optional independent vacuum wave; '
                     'not production acceleration or complete physical acceptance'}
    save(root / 'plan.json', plan)
    frozen = [p for p in root.rglob('*') if p.is_file()]
    frozen += [binary, binding] + runtime_files
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
            if case['arguments'][0] == 'verify':
                require(type(metrics.get('checks')) is int and metrics['checks'] > 0,
                        'No state checks executed')
                require(type(metrics.get('bitwise_failures')) is int and
                        metrics['bitwise_failures'] == 0, 'State-copy equality gate failed')
            elif case['arguments'][0] in WAVE_CASE_NAMES:
                require(metrics.get('mode') == case['arguments'][0] and metrics.get('finite') is True,
                        'Incomplete finite vacuum-wave comparison')
                require(metrics.get('initial_overlap_values_checked', 0) > 0 and
                        metrics.get('initial_overlap_bitwise_failures') == 0,
                        'Vacuum-wave initial overlap gate failed')
                require(metrics.get('reference_cpml_nonzero_seen') is True and
                        metrics.get('candidate_cpml_nonzero_seen') is True,
                        'Vacuum-wave comparison did not exercise CPML history')
                require(metrics.get('plane_frames') == metrics.get('steps', -1) + 1,
                        'Incomplete vacuum-wave history')
                with (directory / 'wave-plane.csv').open() as stream:
                    require(sum(1 for _ in stream) == metrics['plane_frames'] + 1,
                            'Truncated vacuum-wave plane CSV')
                z_cells = 20 if case['arguments'][0].endswith('20') else 8
                require(metrics.get('cpml_cells') == [6, 0, z_cells],
                        'Vacuum-wave CPML thickness control changed')
                require(metrics.get('full_nz') == 512 and
                        metrics.get('full_origin_global_k') == -192,
                        'Vacuum-wave large reference geometry changed')
                require(metrics.get('cpml_polynomial_order') == 3 and
                        metrics.get('cpml_kappa_max') == 8 and
                        metrics.get('cpml_alpha_fraction') == 0 and
                        metrics.get('cpml_target_reflection') == 1e-8,
                        'Vacuum-wave original CPML grading parameters changed')
                if case['arguments'][0].startswith('wave-static'):
                    require(metrics.get('window_nz') == 296 + z_cells - 8 and
                            metrics.get('window_initial_origin_global_k') == -192 and
                            metrics.get('window_final_origin_global_k') == -192 and
                            metrics.get('shifts') == 0 and metrics.get('window_beta') == 0,
                            'Fixed-CPML control geometry or immobility gate failed')
                    require(metrics.get('window_initial_cpml_inner_upper_global_k') == 96,
                            'Fixed-CPML upper entrance changed')
                    for key in ('field_moved_bytes', 'incoming_zero_bytes',
                                'history_moved_bytes', 'history_zeroed_bytes'):
                        require(metrics.get(key) == 0, 'Fixed control remapped state: ' + key)
                else:
                    require(metrics.get('window_nz') == 128 + z_cells - 8 and
                            metrics.get('window_initial_origin_global_k') == 0 and
                            metrics.get('window_initial_cpml_inner_upper_global_k') == 120 and
                            metrics.get('window_beta') == .8,
                            'Moving-wave upper entrance or geometry changed')
            else:
                require(metrics.get('numerical_status') == 'PASS',
                        'Benchmark same-width moving-state equality gate failed')
            results.append({'case': case['name'], 'metrics_sha256': digest(metrics_path)})
        save(root / 'COMPLETED.json', {'unix_time': time.time(), 'cases': results,
             'state_and_cost_probe_gates': 'COMPLETED',
             'physical_acceptance': 'Not assessed: no complete moving-boundary or production acceptance',
             'performance_scope': 'Same-width benchmark ratios and optional different-domain wave timings have separate scopes'})
    print(root / 'COMPLETED.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--binary', type=Path, default=WORKSPACE /
                      'build/moving-window-20261027/verification/aprl_moving_window_probe')
    prep.add_argument('--case', action='append', choices=CASE_NAMES,
                      help='Optional bounded subset, including wave/wave-static and their 20-cell controls; default: verify, bench64, bench128, bench256')
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
