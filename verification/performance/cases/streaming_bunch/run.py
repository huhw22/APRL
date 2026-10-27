#!/usr/bin/env python3
"""Freeze and run bounded charged-bunch window research cases on gnode1."""
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


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def prepare(args):
    root, binary = args.root.resolve(), args.binary.resolve()
    require(not root.exists(), 'Use a fresh result directory')
    require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable')
    cases = []
    specifications = args.case or [
        f'b{label}-{mode}-g{guard}:{beta}:{guard}:{mode}'
        for label, beta in [('0', '0'), ('023', '.23')]
        for mode, guard in [('full', 16), ('mur', 16), ('mur', 24),
                            ('mur', 40), ('clamp', 16), ('zero', 16)]]
    for specification in specifications:
        name, beta, guard, mode = specification.split(':')
        require(name.replace('-', '').isalnum(), 'Invalid case name')
        require(mode in ('full', 'mur', 'clamp', 'zero', 'mur-gauss'), 'Unknown mode')
        require(float(beta) in (0., .23) and 8 <= int(guard) <= 80, 'Outside small-test plan')
        require(name not in [c['name'] for c in cases], 'Duplicate case')
        cases.append({'name': name, 'beta': beta, 'guard': int(guard), 'mode': mode})
    root.mkdir(parents=True)
    shutil.copy2(__file__, root/'run.py')
    source = WORKSPACE/'source/aprl'
    files = [p for p in (source/'src').glob('*') if p.is_file()]
    files += [source/'verification'/n for n in
              ('streaming_bunch_probe.cpp', 'streaming_reference_boundary.cpp',
               'streaming_reference_boundary.h', 'CMakeLists.txt')]
    save(root/'source-files.json', {str(p): digest(p) for p in files})
    (root/'source.diff').write_bytes(subprocess.check_output(['git', '-C', str(source), 'diff']))
    for p in files:
        if p.parent.name == 'verification':
            shutil.copy2(p, root/p.name)
    plan = {'cases': cases, 'binary': str(binary), 'binary_sha256': digest(binary),
            'source_base': subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip(),
            'ranks': 1, 'threads': 1, 'case_timeout_s': 180, 'total_timeout_s': 7200,
            'scope': 'Experimental reference-fed moving window; not production CPML or speed acceptance'}
    save(root/'plan.json', plan)
    frozen = [root/'run.py', root/'plan.json', root/'source-files.json', root/'source.diff', binary]
    frozen += [root/p.name for p in files if p.parent.name == 'verification']
    save(root/'frozen.json', {str(p): digest(p) for p in frozen})
    print(root/'plan.json')


def run(args):
    root = args.root.resolve()
    require(socket.gethostname().split('.')[0] == 'gnode1' and
            os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1' and
            not os.environ.get('SLURM_JOB_ID'), 'Requires authorized gnode1 CPU session')
    require(not (root/'STARTED.json').exists(), 'Campaign already started')
    for name, sha in json.loads((root/'frozen.json').read_text()).items():
        require(digest(Path(name)) == sha, 'Frozen input changed: ' + name)
    plan = json.loads((root/'plan.json').read_text())
    lock = (WORKSPACE/'results/node_survey/.cpu-calibration.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    other = subprocess.run(['pgrep', '-a', '-x', 'aprl|mpirun|mpiexec|orted'], text=True, capture_output=True)
    require(other.returncode == 1, 'Existing compute process: ' + other.stdout)
    def counters():
        a = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(a), a[3]
    before, idle_before = counters()
    time.sleep(2)
    after, idle_after = counters()
    idle = (idle_after-idle_before)/(after-before)
    require(idle >= .95, 'Test host is busy')
    save(root/'STARTED.json', {'host': socket.gethostname(), 'unix_time': time.time(),
                               'cpu_idle_fraction': idle, 'loadavg': os.getloadavg()})
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
               OMP_PROC_BIND='close', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               APRL_PROFILE_INTERVAL='0')
    launcher = ['mpiexec', '-n', '1', '--map-by', 'slot:PE=1', '--bind-to', 'core', '--nooversubscribe']
    binding = subprocess.run(['timeout', '--signal=TERM', '--kill-after=5', '45'] + launcher +
                [str(WORKSPACE/'build/node-survey/binding_probe'), '1', '1'],
                env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=55)
    (root/'binding.log').write_text(binding.stdout)
    require(binding.returncode == 0 and '[binding] result=PASS failures=0' in binding.stdout,
            'Binding check failed')
    deadline = time.monotonic() + plan['total_timeout_s']
    for case in plan['cases']:
        require(time.monotonic() + plan['case_timeout_s'] < deadline, 'Campaign deadline')
        directory = root/case['name']
        directory.mkdir()
        command = launcher + [plan['binary'], str(directory), case['beta'], str(case['guard']), case['mode']]
        start = time.monotonic()
        with (directory/'probe.log').open('w') as log:
            process = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=plan['case_timeout_s'])
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                save(directory/'launch.json', {'command': command, 'timeout': True})
                raise RuntimeError('Case timed out: ' + case['name'])
        save(directory/'launch.json', {'command': command, 'exit_code': code,
                                      'seconds': time.monotonic()-start})
        require(code == 0, 'Probe failed: ' + case['name'])
        metrics = json.loads((directory/'metrics.json').read_text())
        require(all(math.isfinite(v) for v in metrics.values() if isinstance(v, (int, float))),
                'Nonfinite metrics: ' + case['name'])
        require(metrics['completed'] and metrics['admitted'] == metrics['particles'], 'Incomplete cohort')
        require(metrics['max_continuity_relative'] <= 1e-10, 'Continuity gate failed')
        require(metrics['max_open_charge_balance_C']/abs(metrics['charge_C']) <= 1e-10,
                'Open-volume charge gate failed')
        require(metrics['max_remap_charge_balance_C']/abs(metrics['charge_C']) <= 1e-10,
                'Moved-grid charge gate failed')
        require(metrics['initial_gauss_relative'] <= 4e-10, 'Cowan initialization gate failed')
    save(root/'COMPLETED.json', {'unix_time': time.time(), 'cases': len(plan['cases']),
                                 'physical_acceptance': 'Requires independent review of all measured differences'})
    print(root/'COMPLETED.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--binary', type=Path, default=WORKSPACE/'build/streaming-bunch-20261027/verification/aprl_streaming_bunch_probe')
    prep.add_argument('--case', action='append', help='name:beta:guard:mode; omit for the planned twelve cases')
    run_parser = sub.add_parser('run')
    run_parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    try:
        globals()[args.action](args)
    except Exception as error:
        if args.action == 'run':
            save(args.root/'FAILED.json', {'message': str(error), 'unix_time': time.time()})
        raise


if __name__ == '__main__':
    main()
