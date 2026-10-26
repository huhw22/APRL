#!/usr/bin/env python3
"""Prepare or explicitly run bounded checkpoint failures; never overwrite results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import struct
import subprocess
import time

import yaml

HERE = Path(__file__).resolve().parent


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def save_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def card(directory, ranks, bad_parent=False):
    checkpoint = directory / ('not-a-directory/child' if bad_parent else 'snapshots')
    return {
        'units': {'length': 'micrometer', 'time': 'femtosecond'},
        'runtime': {
            'mode': 'throughput', 'resource_monitor': {'enabled': False},
            'shutdown_watchdog': {'enabled': True, 'timeout_seconds': 600,
                                  'finalize_timeout_seconds': 60},
        },
        'output': {'overwrite': False, 'manifest': str(directory / 'manifest.yaml')},
        'mesh': {'field_solver': 'cowan-z', 'current_filter': 'none',
                 'cells': [8, 8, max(64, 2 * ranks)],
                 'cell_size': [0.5, 0.5, 0.1], 'center': [0., 0., 0.],
                 'duration': 5., 'boost_gamma': 2.,
                 'particle_steps_per_undulator_period': 1},
        'boundary': {'type': 'pec'}, 'initial_self_field': {'enabled': False},
        'beam': {
            'reference': {'initial_center_z': -0.5},
            'input': {'type': 'generated-gaussian', 'electrons': 8.,
                      'macroparticles': 8, 'gamma': 4., 'direction': [0., 0., 1.],
                      'center': [0., 0., 0.], 'sigma_position': [0., 0., 0.],
                      'sigma_proper_velocity': [0., 0., 0.], 'random_seed': 11},
        },
        'sources': {'incident_waves': [], 'magnetic_elements': []},
        'stop': {'mode': 'reference-center-z', 'z': 0.5},
        'trajectory': {'enabled': False},
        'checkpoint': {'directory': str(checkpoint), 'steps': [1],
                       'io_timeout_seconds': 3},
    }


def prepare(args):
    output = args.directory.resolve()
    require(not output.exists(), 'Result directory already exists; choose a new directory')
    binary, library = args.binary.resolve(), args.library.resolve()
    require(binary.is_file() and os.access(binary, os.X_OK), 'APRL binary is not executable')
    require(library.is_file(), 'Build the test preload library before preparing')
    launcher = json.loads(args.launcher_json)
    require(isinstance(launcher, list) and launcher and
            all(isinstance(item, str) and item for item in launcher),
            'launcher-json must be a nonempty argument array, not a shell command')
    require(args.ranks >= 2 and args.ranks <= 256, 'Fault test needs 2..256 ranks')
    require(0 <= args.target_rank < args.ranks, 'Target rank is outside the allocation')
    output.mkdir(parents=True)
    cases = []
    for mode in ('create', 'write', 'bad-parent'):
        directory = output / mode
        directory.mkdir()
        (directory / 'markers').mkdir()
        if mode == 'bad-parent':
            (directory / 'not-a-directory').write_text('Intentional regular-file parent.\n')
        configuration = directory / 'config.yaml'
        configuration.write_text(yaml.safe_dump(card(directory, args.ranks,
                                                    mode == 'bad-parent'), sort_keys=False))
        cases.append({'mode': mode, 'directory': str(directory),
                      'configuration': str(configuration), 'sha256': digest(configuration)})
    plan = {'schema': 'aprl-checkpoint-fault-plan-v1', 'prepared_at': time.time(),
            'binary': str(binary), 'binary_sha256': digest(binary),
            'library': str(library), 'library_sha256': digest(library),
            'source_sha256': digest(HERE / 'fault.c'),
            'launcher': launcher, 'ranks': args.ranks, 'target_rank': args.target_rank,
            'checkpoint_timeout_seconds': 3, 'injected_block_seconds': 20,
            'launcher_timeout_seconds': 15, 'cases': cases}
    save_new(output / 'PLAN.json', plan)
    print(json.dumps({'status': 'PREPARED_NOT_RUN', 'plan': str(output / 'PLAN.json')}))


def marker_evidence(directory, ranks):
    records = []
    for rank in range(ranks):
        stem = directory / 'markers' / ('rank-%06d' % rank)
        entry = Path(str(stem) + '.enter')
        exit_file = Path(str(stem) + '.exit')
        item = {'rank': rank, 'entered': False, 'local_exit_124': False}
        if entry.exists() and entry.stat().st_size == 8:
            item['entered'] = True
            item['entered_monotonic_ns'] = struct.unpack('=Q', entry.read_bytes())[0]
        if exit_file.exists() and exit_file.stat().st_size == 16:
            magic, status, when = struct.unpack('=IiQ', exit_file.read_bytes())
            item.update({'exit_magic': magic, 'exit_status': status,
                         'exited_monotonic_ns': when,
                         'local_exit_124': magic == 0x4150524c and status == 124})
            if item['entered']:
                item['after_entry_seconds'] = (when - item['entered_monotonic_ns']) / 1e9
        records.append(item)
    return records


def launch_case(plan, case, timeout_program):
    directory = Path(case['directory'])
    mode = case['mode']
    environment = os.environ.copy()
    for key in list(environment):
        if key.startswith('APRL_CHECKPOINT_FAULT_'):
            del environment[key]
    # Existing preloads could obscure the evidence; this is an isolated test process.
    require(not environment.get('LD_PRELOAD'), 'Unset existing LD_PRELOAD for this isolated test')
    environment['OMP_NUM_THREADS'] = '1'
    if mode != 'bad-parent':
        environment.update({
            'LD_PRELOAD': plan['library'], 'APRL_CHECKPOINT_FAULT_ENABLE': '1',
            'APRL_CHECKPOINT_FAULT_PHASE': mode,
            'APRL_CHECKPOINT_FAULT_BINARY': plan['binary'],
            'APRL_CHECKPOINT_FAULT_RANK': str(plan['target_rank']),
            'APRL_CHECKPOINT_FAULT_MARKERS': str(directory / 'markers'),
        })
    command = [timeout_program, '--kill-after=3s', '15s'] + plan['launcher'] + [
        plan['binary'], case['configuration']]
    start = time.monotonic()
    outer_timeout = False
    with (directory / 'aprl.log').open('xb') as log:
        process = subprocess.Popen(command, cwd=directory, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            returncode = process.wait(timeout=25)
        except subprocess.TimeoutExpired:
            outer_timeout = True
            os.killpg(process.pid, signal.SIGKILL)
            returncode = process.wait(timeout=5)
    elapsed = time.monotonic() - start
    log_text = (directory / 'aprl.log').read_text(errors='replace')
    published = [str(path) for path in directory.rglob('step-*.h5')]
    partial = [str(path) for path in directory.rglob('step-*.h5.partial')]
    result = {'mode': mode, 'command': command, 'returncode': returncode,
              'elapsed_seconds': elapsed, 'outer_timeout': outer_timeout,
              'published_files': published, 'partial_files': partial}
    if mode != 'bad-parent':
        markers = marker_evidence(directory, plan['ranks'])
        injected_file = directory / 'markers/injected.json'
        injection = json.loads(injected_file.read_text()) if injected_file.exists() else None
        all_local = all(item['entered'] and item['local_exit_124'] and
                        0 <= item.get('after_entry_seconds', 100) < 10 for item in markers)
        result.update({'injection': injection, 'rank_evidence': markers,
                       'all_ranks_local_watchdog_exit': all_local})
        result['pass'] = bool(not outer_timeout and not published and returncode != 0 and
                              elapsed < 14 and all_local and injection and
                              injection.get('rank') == plan['target_rank'] and
                              injection.get('phase') == mode and
                              injection.get('block_seconds') == 20)
    else:
        explicit_error = ('Cannot create directory' in log_text or
                          'Output directory creation failed' in log_text or
                          'Not a directory' in log_text)
        result['explicit_directory_error'] = explicit_error
        result['pass'] = bool(not outer_timeout and elapsed < 14 and
                              returncode not in (0, 124, 137) and explicit_error and
                              not published and not partial)
    save_new(directory / 'RESULT.json', result)
    return result


def run(args):
    output = args.directory.resolve()
    plan = json.loads((output / 'PLAN.json').read_text())
    host = socket.gethostname().split('.')[0]
    require(os.environ.get('SLURM_JOB_ID') or host == 'gnode1',
            'Run computation only inside Slurm or the explicitly authorized gnode1 exception')
    require(host != 'gnode1' or plan['ranks'] <= 16, 'gnode1 budget is 16 physical CPU cores')
    require(digest(plan['binary']) == plan['binary_sha256'], 'APRL binary changed since preparation')
    require(digest(plan['library']) == plan['library_sha256'], 'Preload library changed since preparation')
    for case in plan['cases']:
        require(digest(case['configuration']) == case['sha256'], 'Prepared configuration changed')
    timeout_program = shutil.which('timeout')
    require(timeout_program, 'GNU timeout is required')
    save_new(output / 'RUN_STARTED.json', {'time': time.time(), 'host': host,
                                         'slurm_job_id': os.environ.get('SLURM_JOB_ID')})
    results = []
    try:
        for case in plan['cases']:
            results.append(launch_case(plan, case, timeout_program))
        result = {'status': 'PASS' if all(item['pass'] for item in results) else 'FAIL',
                  'cases': results, 'binary_sha256': plan['binary_sha256'],
                  'note': 'Strict PASS requires actual _exit(124) evidence from every rank; '
                          'launcher cleanup alone is not accepted as local watchdog proof.'}
    except Exception as error:
        result = {'status': 'FAIL', 'error': str(error), 'cases': results}
    save_new(output / 'RESULT.json', result)
    (output / result['status']).write_text('\n')
    print(json.dumps({'status': result['status'], 'result': str(output / 'RESULT.json')}))
    return 0 if result['status'] == 'PASS' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--directory', required=True, type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--library', type=Path)
    parser.add_argument('--launcher-json')
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--target-rank', type=int, default=1)
    args = parser.parse_args()
    if args.action == 'prepare':
        require(args.binary and args.library and args.launcher_json,
                'prepare requires binary, library and launcher-json')
        prepare(args)
        return 0
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
