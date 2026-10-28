#!/usr/bin/env python3
"""Prepare/run two-node candidate replay against already accepted one-node data.

Prepare: mpi_review.py --prepare --parent PAIR_ROOT --root NEW_MPI_ROOT
Run inside a two-node Slurm allocation: mpi_review.py --root NEW_MPI_ROOT
Never submits jobs or changes the completed reference campaign.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

import h5py
import numpy as np
import yaml

RANKS, THREADS, COUNT = 16, 8, 31675
HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok: raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): value.update(block)
    return value.hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


def physical_card(card):
    result = copy.deepcopy(card)
    result['output']['manifest'] = '<OUTPUT>'
    result['mesh']['moving_window']['reference_directory'] = '<REFERENCE>'
    for section in ('detectors', 'energy_ledger', 'trajectory'):
        result[section]['directory'] = '<OUTPUT>'
    return result


def parent_gate(parent):
    plan = json.loads((parent / 'plan.json').read_text())
    path = parent / 'SCIENCE_PASS.json'
    gate = json.loads(path.read_text())
    require(gate['pass'] and gate['stage'] == 'science', 'One-node science has not passed')
    require(gate['binary_sha256'] == plan['binary_sha256'] and
            gate['frozen_manifest_sha256'] == sha(parent / 'frozen.json'),
            'One-node science identity differs')
    for magnet in ('undulator', 'free-drift'):
        name = magnet + '-16x8-candidate'
        require(gate['cases'][name]['pass'] and gate['comparisons'][name]['pass'],
                'Missing accepted one-node candidate: ' + name)
    return plan, gate


def prepare(parent, root):
    source_plan, gate = parent_gate(parent)
    require(not root.exists(), 'Use a fresh two-node directory')
    binary = Path(source_plan['binaries']['candidate'])
    require(sha(binary) == source_plan['binary_sha256']['candidate'], 'Candidate binary changed')
    root.mkdir(parents=True); scripts = root / 'scripts'; scripts.mkdir()
    for name in ('mpi_review.py', 'mpi.sbatch'):
        shutil.copy2(HERE / name, scripts / name)
    for name in ('smoke_review.py', 'ledger_review.py'):
        shutil.copy2(parent / 'scripts' / name, scripts / name)
    shutil.copy2(HERE / 'mpi.sbatch', root / 'mpi.sbatch')
    cases, reference_stats = [], {}
    for magnet in ('undulator', 'free-drift'):
        reference = parent / 'science' / (magnet + '-16x8-candidate')
        directory = root / magnet; directory.mkdir()
        original = yaml.safe_load((reference / 'run.yaml').read_text())
        card = copy.deepcopy(original)
        card['output']['manifest'] = str(directory / 'run-manifest.yaml')
        card['mesh']['moving_window']['reference_directory'] = str(directory / 'reference')
        for section, folder in (('detectors', 'detectors'), ('energy_ledger', 'ledger'), ('trajectory', 'trajectory')):
            card[section]['directory'] = str(directory / folder)
        require(physical_card(card) == physical_card(original), 'Physics/input changed while redirecting')
        (directory / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
        for path in sorted(reference.rglob('*.h5')):
            stat = path.stat()
            reference_stats[str(path)] = {'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
        cases.append({'name': magnet, 'directory': str(directory), 'reference': str(reference),
                      'steps': gate['cases'][magnet + '-16x8-candidate']['steps'], 'timeout_s': 240})
    plan = {'format': 'aprl-pending-two-node-v1', 'parent': str(parent),
        'parent_science_pass': str(parent / 'SCIENCE_PASS.json'),
        'parent_science_pass_sha256': sha(parent / 'SCIENCE_PASS.json'),
        'parent_frozen_manifest_sha256': sha(parent / 'frozen.json'),
        'binary': str(binary), 'binary_sha256': sha(binary),
        'binding_probe': source_plan['binding_probe'], 'nodes': 2, 'ranks': RANKS,
        'ranks_per_node': 8, 'threads_per_rank': THREADS, 'used_physical_cores': 128,
        'exclusive_allocated_physical_cores': 256, 'cases': cases,
        'small_input_file': str(parent / 'inputs' / 'held.h5'),
        'small_input_sha256': source_plan['small_input_sha256'], 'records': COUNT,
        'mpi_environment': source_plan['mpi_environment'], 'profile_interval': 0,
        'active_budget_s': 720, 'slurm_limit_minutes': 15,
        'reference_h5_stat': reference_stats,
        'acceptance': {'field_relative_l2': 2e-12, 'particle_roundoff_factor_eps64': 128,
                       'ledger': json.loads((parent / 'THRESHOLDS.json').read_text())['science_ledger']},
        'limits': ['Same candidate, same 16 MPI ranks and eight threads per rank; only host distribution changes.',
            'Two exclusive nodes reserve 256 cores but the steps use 128 cores; this is functional MPI acceptance, not a speedup benchmark.',
            'No baseline or full-input simulation is launched.',
            'Bitwise results are reported observations; existing strict scientific gates govern acceptance.',
            'Pending axial ledger statistics retain the declared affine diagnostic approximation.']}
    save(root / 'plan.json', plan)
    frozen = [root / 'plan.json', root / 'mpi.sbatch', binary, Path(plan['binding_probe']),
              parent / 'SCIENCE_PASS.json', parent / 'frozen.json', Path(plan['small_input_file'])]
    frozen += list(scripts.iterdir())
    frozen += [path for path in (binary.parent / 'runtime-libs').iterdir() if path.is_file()]
    for item in cases:
        frozen += [Path(item['directory']) / 'run.yaml', Path(item['reference']) / 'run.yaml',
                   Path(item['reference']) / 'run-manifest.yaml', Path(item['reference']) / 'REVIEW.json']
    save(root / 'frozen.json', {str(path): sha(path) for path in frozen})
    print(root / 'plan.json')


def check_frozen(root, plan):
    for name, expected in json.loads((root / 'frozen.json').read_text()).items():
        require(sha(name) == expected, 'Frozen input changed: ' + name)
    parent_gate(Path(plan['parent']))
    for name, expected in plan['reference_h5_stat'].items():
        stat = Path(name).stat()
        require((stat.st_size, stat.st_mtime_ns) == (expected['bytes'], expected['mtime_ns']),
                'Accepted reference output changed: ' + name)


def launch(command, directory, env, limit, name='aprl.log'):
    command = ['timeout', '--signal=TERM', '--kill-after=10', str(limit)] + command
    record = {'command': command, 'cwd': str(directory), 'limit_s': limit,
              'started_unix': time.time(), 'exit_code': None,
              'environment': {k: env.get(k) for k in ('SLURM_JOB_ID', 'OMP_NUM_THREADS', 'OMP_PLACES',
                  'OMP_PROC_BIND', 'OMP_DYNAMIC', 'OMPI_MCA_btl', 'APRL_PROFILE_INTERVAL', 'LD_LIBRARY_PATH')}}
    process = None; start = time.monotonic()
    try:
        with (directory / name).open('x') as stream:
            process = subprocess.Popen(command, cwd=directory, env=env, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            record['exit_code'] = process.wait(timeout=limit + 15)
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try: process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
        record['seconds'] = time.monotonic() - start
        save(directory / (name + '.launch.json'), record)
    require(record['exit_code'] == 0, 'MPI step failed/timed out: ' + str(directory / name))
    return record


def binding(path):
    text = path.read_text()
    rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)', text)
    require('[binding] result=PASS failures=0' in text and len(rows) == 128 and
            all(r[-1] == 'PASS' for r in rows), 'Incomplete or failed binding')
    require({(int(r[0]), int(r[1])) for r in rows} ==
            {(r, t) for r in range(RANKS) for t in range(THREADS)}, 'Rank/thread binding coverage differs')
    hosts = sorted({r[2] for r in rows})
    require(len(hosts) == 2 and len({(r[2], r[4], r[5]) for r in rows}) == 128,
            'Need two hosts with 128 distinct physical cores')
    host_ranks = {host: sorted({int(r[0]) for r in rows if r[2] == host}) for host in hosts}
    require(sorted(host_ranks.values()) == [list(range(8)), list(range(8, 16))], 'Expected block rank placement across the inter-node rank7/8 interface')
    require(all(sum(r[2] == host for r in rows) == 64 for host in hosts), 'Expected 64 used cores per host')
    return {'pass': True, 'hosts': hosts, 'host_ranks': host_ranks, 'used_physical_cores': 128,
            'reserved_physical_cores': 256, 'sha256': sha(path)}


def transactions_and_identity(directory):
    manifest = yaml.safe_load((directory / 'run-manifest.yaml').read_text())
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    require(manifest['mpi_size'] == RANKS and yaml.safe_load(manifest['configuration_yaml']) == card,
            'Manifest rank count/card differs')
    watchdog = manifest['resolved_runtime']['shutdown_watchdog']
    require(watchdog['enabled'] and watchdog['timeout_seconds'] == 600 and
            watchdog['finalize_timeout_seconds'] == 60, 'Shutdown watchdog changed')
    for path in sorted(directory.rglob('*.h5')):
        with h5py.File(path) as file:
            require(len(file) == 1, 'Unexpected diagnostic groups: ' + str(path))
            group = file[next(iter(file))]
            for key in ('run_id', 'source_revision', 'configuration_digest_fnv1a64'):
                value = group.attrs[key]
                value = value.decode() if isinstance(value, bytes) else str(value)
                require(value == str(manifest[key]), 'Diagnostic identity mismatch: ' + str(path) + '/' + key)
    return manifest


def trajectory_comparison(reference, candidate, helper):
    def ordered(directory):
        arrays = []
        for path in sorted((directory / 'trajectory').glob('*.h5')):
            with h5py.File(path) as file: arrays.append(file['trajectory/records'][:])
        records = np.concatenate(arrays)
        return records[np.lexsort((records['time_s'], records['event_type'], records['source_id']))]
    a, b = ordered(reference), ordered(candidate)
    require(a.dtype == b.dtype and a.shape == b.shape, 'Trajectory cohort/sample coverage differs')
    for key in ('particle_id', 'source_id', 'charge_C', 'weight', 'event_type', 'boundary_face'):
        require(np.array_equal(a[key], b[key]), 'Trajectory identity/event mismatch: ' + key)
    metrics = {key: helper.particle_error(a[key], b[key]) for key in ('position_m', 'proper_velocity', 'time_s')}
    return {'pass': all(value['pass'] for value in metrics.values()), 'records': len(a),
            'ordering': 'source_id, event_type, time-order within each event; independent of rank file ownership',
            'metrics': metrics, 'bitwise_equal': a.tobytes() == b.tobytes()}


def moving_csv(directory, expected_step):
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    log = (directory / 'aprl.log').read_text()
    geometry = re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+) ranks=16\b', log)
    require(geometry is not None, 'Missing moving-window geometry')
    reference_origin = card['mesh']['center'][2] - .5 * int(geometry[1]) * card['mesh']['cell_size'][2]
    path = Path(card['mesh']['moving_window']['reference_directory'])
    expected_files = {path / ('window-rank-%d.csv' % rank) for rank in range(RANKS)}
    require(set(path.glob('window-rank-*.csv')) == expected_files, 'Missing/unexpected rank window CSV')
    finals = []
    interval = card['energy_ledger']['sample_interval_steps']
    steps = list(range(interval, expected_step + 1, interval)) + [expected_step]
    for rank in range(RANKS):
        rows = np.atleast_1d(np.genfromtxt(path / ('window-rank-%d.csv' % rank), names=True, delimiter=','))
        require(all(np.isfinite(rows[key]).all() for key in rows.dtype.names), 'Nonfinite moving ledger')
        require(np.array_equal(rows['step'], steps) and np.all(rows['energy_measured'] == 1), 'Moving ledger step coverage differs')
        origin = np.longdouble(reference_origin) + (np.longdouble(int(geometry[3])) - rows['step'].astype(np.longdouble)) * np.longdouble(card['mesh']['cell_size'][2])
        require(np.max(abs(rows['origin_box_z_m'].astype(np.longdouble) - origin)) <= 128 * np.finfo(float).eps,
                'Window did not shift by exactly one cell per step')
        finals.append({key: float(rows[-1][key]) for key in rows.dtype.names})
    completed = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)', log)
    require(completed is not None and int(completed[1]) == expected_step and
            sum(row['particle_exits'] for row in finals) == int(completed[2]) and
            sum(row['cropped_carriers'] for row in finals) == int(completed[3]), 'Moving ledger/log counts differ')
    return {'pass': True, 'ranks': RANKS, 'step': expected_step, 'final_rank_rows': finals}


def run(root):
    plan = json.loads((root / 'plan.json').read_text())
    parent_gate(Path(plan['parent']))  # Must precede all srun/dependency work.
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == 2 and
            int(os.environ.get('SLURM_NTASKS', '0')) == RANKS and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == THREADS,
            'Require a real two-node, 16-rank, eight-thread Slurm allocation')
    require(plan['format'] == 'aprl-pending-two-node-v1', 'Unexpected MPI campaign')
    execution = root / 'execution'; execution.mkdir(exist_ok=True)
    require(not (execution / 'STARTED.json').exists(), 'MPI stage already started; use a fresh directory')
    check_frozen(root, plan)
    helper = module(root / 'scripts/smoke_review.py', 'pending_mpi_smoke')
    ledger = module(root / 'scripts/ledger_review.py', 'pending_mpi_ledger')
    require(ledger.thresholds(COUNT) == plan['acceptance']['ledger'], 'Frozen ledger precision contract differs')
    helper.RANKS = RANKS
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='8', OMP_DYNAMIC='FALSE', OMP_PLACES='cores', OMP_PROC_BIND='close',
               APRL_PROFILE_INTERVAL='0', HDF5_USE_FILE_LOCKING='FALSE', **plan['mpi_environment'])
    binary = Path(plan['binary'])
    env['LD_LIBRARY_PATH'] = str(binary.parent / 'runtime-libs') + ':' + env.get('LD_LIBRARY_PATH', '')
    start = time.monotonic(); deadline = start + plan['active_budget_s']
    save(execution / 'STARTED.json', {'slurm_job_id': os.environ['SLURM_JOB_ID'], 'started_unix': time.time(),
        'binary_sha256': plan['binary_sha256'], 'active_budget_s': plan['active_budget_s']})
    def interrupted(signum, frame): raise TimeoutError('MPI validation interrupted or active budget exhausted')
    signal.signal(signal.SIGALRM, interrupted); signal.signal(signal.SIGTERM, interrupted)
    signal.setitimer(signal.ITIMER_REAL, plan['active_budget_s'])
    def launcher(directory):
        return ['srun', '--mpi=' + env.get('SRUN_MPI', 'pmix_v3'), '--exact', '--nodes=2', '--ntasks=16',
                '--ntasks-per-node=8', '--cpus-per-task=8', '--distribution=block:block', '--cpu-bind=cores',
                '--kill-on-bad-exit=1', '--wait=660', '--chdir=' + str(directory)]
    reports = {}
    try:
        linked = subprocess.run(['ldd', '-r', str(binary)], env=env, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=30)
        save(execution / 'RUNTIME_LINKS.json', {'exit_code': linked.returncode, 'output': linked.stdout})
        require(linked.returncode == 0 and 'not found' not in linked.stdout and 'undefined symbol' not in linked.stdout,
                'Candidate runtime dependency failure')
        require(sha(plan['small_input_file']) == plan['small_input_sha256'], 'Small particle input changed')
        binddir = execution / 'binding'; binddir.mkdir()
        launch(launcher(binddir) + [plan['binding_probe'], '16', '8'], binddir, env, 45, 'binding.log')
        bound = binding(binddir / 'binding.log'); save(execution / 'BINDING.json', bound)
        for item in plan['cases']:
            directory, reference = Path(item['directory']), Path(item['reference'])
            require(deadline - time.monotonic() >= item['timeout_s'] + 30, 'Insufficient budget for complete next case')
            require(physical_card(yaml.safe_load((directory / 'run.yaml').read_text())) ==
                    physical_card(yaml.safe_load((reference / 'run.yaml').read_text())), 'Reference/input/physics differ')
            launched = launch(launcher(directory) + [str(binary), str(directory / 'run.yaml')], directory, env, item['timeout_s'])
            log = (directory / 'aprl.log').read_text(errors='replace')
            require(not re.search(r'CANCELLED|DUE TO TIME LIMIT|mesh.duration was exhausted|Simulation failure', log), 'Run did not finish normally')
            queue = re.search(r'Pending admission queue: waiting=(\d+) resident=(\d+)', log)
            require(queue is not None and int(queue[1]) > 0, 'Test did not exercise scheduled upstream particles')
            reference_summary, unused = helper.inspect_case(reference, free=item['name'] == 'free-drift')
            summary, planes = helper.inspect_case(directory, free=item['name'] == 'free-drift')
            require(summary['steps'] == reference_summary['steps'] == item['steps'], 'Physical stopping step changed')
            manifest = transactions_and_identity(directory)
            reference_manifest = transactions_and_identity(reference)
            require(manifest['source_revision'] == reference_manifest['source_revision'], 'Compiled candidate revision changed')
            comparison = helper.compare(reference, directory, planes)
            comparison['ledger'] = ledger.compare(reference, directory, COUNT)
            comparison['trajectory'] = trajectory_comparison(reference, directory, helper)
            comparison['pass'] &= comparison['ledger']['pass'] and comparison['trajectory']['pass']
            require(comparison['pass'], 'Same-rank cross-node scientific replay failed')
            window = moving_csv(directory, item['steps']); reference_window = moving_csv(reference, item['steps'])
            for a, b in zip(reference_window['final_rank_rows'], window['final_rank_rows']):
                for key in ('step', 'origin_box_z_m', 'particle_exits', 'cropped_carriers'):
                    require(a[key] == b[key], 'Window final geometry/count changed: ' + key)
            reports[item['name']] = {'pass': True, 'launcher': launched, 'initial_pending': int(queue[1]),
                'reference': str(reference), 'summary': summary, 'comparison': comparison, 'window': window}
            save(directory / 'MPI_REVIEW.json', reports[item['name']])
            save(execution / 'PROGRESS.json', {'completed': list(reports), 'elapsed_s': time.monotonic() - start})
            print('REVIEWED ' + item['name'], flush=True)
        check_frozen(root, plan)
        result = {'pass': True, 'format': 'aprl-pending-two-node-completion-v1', 'binary_sha256': plan['binary_sha256'],
            'parent_science_pass_sha256': plan['parent_science_pass_sha256'], 'binding': bound, 'cases': reports,
            'elapsed_s': time.monotonic() - start, 'frozen_manifest_sha256': sha(root / 'frozen.json'), 'limits': plan['limits']}
        save(root / 'MPI_PASS.json', result)
        print('MPI_REPLAY_PASS', flush=True)
    except BaseException as error:
        save(execution / 'FAILED_OR_INCOMPLETE.json', {'pass': False, 'error': str(error),
            'completed_cases': list(reports), 'elapsed_s': time.monotonic() - start,
            'scope': 'No incomplete case promoted to success; no automatic resubmission or budget extension.'})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--parent', type=Path)
    parser.add_argument('--prepare', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        require(args.parent is not None, '--prepare requires --parent')
        prepare(args.parent.resolve(), args.root.resolve())
    else:
        run(args.root.resolve())


if __name__ == '__main__': main()
