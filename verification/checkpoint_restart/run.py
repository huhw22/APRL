#!/usr/bin/env python3
"""Bounded checkpoint correctness campaign; prepare never launches computation."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import socket
import subprocess
import time

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
C = 299792458.0
IDENTITY_ATTRIBUTES = {
    'run_id', 'run_manifest', 'source_revision', 'configuration_digest_fnv1a64',
    'restart_parent', 'restart_step', 'restart_time_box_s', 'output_scope',
}


def require(value, message):
    if not value:
        raise RuntimeError(message)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def bitwise_equal(first, second, label):
    a, b = np.asarray(first), np.asarray(second)
    require(a.shape == b.shape and a.dtype == b.dtype,
            label + ': shape/type mismatch')
    if a.dtype.names:
        # Compound padding is not scientific state and can be uninitialized.
        for name in a.dtype.names:
            bitwise_equal(a[name], b[name], label + '.' + name)
    elif a.dtype.kind in 'SUO':
        require(np.array_equal(a, b), label + ': value mismatch')
    else:
        require(a.tobytes() == b.tobytes(), label + ': bitwise mismatch')


def finite(data, label):
    if data.dtype.names:
        for name in data.dtype.names:
            finite(data[name], label + '.' + name)
    elif data.dtype.kind in 'fc':
        require(np.isfinite(data).all(), label + ': nonfinite output')


def inspect_file(path):
    with h5py.File(path) as file:
        groups = [item for item in file.values() if isinstance(item, h5py.Group)]
        require(len(groups) == 1, str(path) + ': unexpected output group count')
        group = groups[0]
        arrays = {}
        for name, item in group.items():
            require(isinstance(item, h5py.Dataset), 'Unexpected nested scientific group')
            arrays[name] = item[()]
            finite(np.asarray(arrays[name]), str(path) + '/' + name)
        require(int(arrays.get('complete', -1)) == 1, str(path) + ': incomplete output')
        for commit, dataset in [('committed_records', 'records'),
                                ('committed_samples', 'time_s'),
                                ('committed_validation_records', 'validation_records')]:
            if commit in arrays:
                require(dataset in arrays and int(arrays[commit]) == len(arrays[dataset]),
                        str(path) + ': committed length mismatch for ' + dataset)
        if 'fields' in arrays:
            require(len(arrays['fields']) == len(arrays['time_s']),
                    str(path) + ': field/time length mismatch')
        return group.name, arrays, dict(group.attrs)


def output_files(directory):
    return {str(path.relative_to(directory)): path
            for folder in ('detectors', 'ledger', 'trajectory')
            for path in (directory / folder).glob('*.h5')}


def snapshot_metadata(path):
    prefixes, time_box = {}, None
    with h5py.File(path) as file:
        require(int(file['header'][2]) == 1, 'Snapshot was not committed')
        block = file['blocks/metadata/diagnostic_prefixes']
        counts, data = block['rank_counts'][:], block['data'][:]
        offset = 0
        for count in counts:
            count = int(count)
            item = yaml.safe_load(data[offset:offset + count].tobytes().decode())
            offset += count
            time_box = item['time_box_s'] if time_box is None else time_box
            require(item['time_box_s'] == time_box, 'Checkpoint rank time mismatch')
            for output in item['outputs']:
                filename = str(Path(output['file']).resolve())
                require(filename not in prefixes, 'Duplicate prefix output path')
                prefixes[filename] = output
        require(offset == len(data), 'Checkpoint prefix rank counts differ from data')
    return prefixes, time_box


def compare_outputs(reference, candidate, snapshot_step=None, snapshot=None):
    expected, actual = output_files(reference), output_files(candidate)
    require(expected.keys() == actual.keys(), 'Scientific output file set differs')
    prefixes, snapshot_time = snapshot_metadata(snapshot) if snapshot else ({}, None)
    checked = []
    for name in sorted(expected):
        group_a, arrays_a, attrs_a = inspect_file(expected[name])
        group_b, arrays_b, attrs_b = inspect_file(actual[name])
        require(group_a == group_b and arrays_a.keys() == arrays_b.keys(),
                name + ': scientific layout differs')
        prefix = prefixes.get(str(expected[name].resolve()))
        if snapshot:
            require(prefix is not None, name + ': no committed checkpoint prefix')
            require(int(attrs_b.get('restart_step', -1)) == snapshot_step and
                    float(attrs_b.get('restart_time_box_s', -1)) == snapshot_time,
                    name + ': invalid restart step/time provenance')
            parent = attrs_b.get('restart_parent', '')
            parent = parent.decode() if isinstance(parent, bytes) else str(parent)
            require(Path(parent).resolve() == snapshot.resolve(), name + ': wrong restart parent')
            scope = attrs_b.get('output_scope', '')
            scope = scope.decode() if isinstance(scope, bytes) else str(scope)
            require(scope.startswith('restart-segment'), name + ': missing segment scope')
        for attr in set(attrs_a) | set(attrs_b):
            if attr in IDENTITY_ATTRIBUTES:
                continue
            require(attr in attrs_a and attr in attrs_b, name + ': missing attribute ' + attr)
            bitwise_equal(attrs_a[attr], attrs_b[attr], name + '@' + attr)
        for dataset, data in arrays_b.items():
            original = arrays_a[dataset]
            label = name + '/' + dataset
            if snapshot_step is None or np.ndim(data) == 0:
                if snapshot_step is not None and dataset.startswith('committed_'):
                    continue  # Already checked against this segment's actual length.
                bitwise_equal(original, data, label)
            elif dataset == 'records' and data.dtype.names and 'step' in data.dtype.names:
                bitwise_equal(original[original['step'] > snapshot_step], data, label)
            else:
                require(prefix is not None, label + ': restart lacks prefix length')
                count = int(prefix['validation_records'] if dataset == 'validation_records'
                            else prefix['records'])
                require(0 <= count <= len(original), label + ': invalid committed prefix')
                tail = original[count:]
                bitwise_equal(tail, data, label)
            checked.append(label)
    require(checked, 'No scientific datasets compared')
    return checked


def inspect_snapshot(path):
    with h5py.File(path) as file:
        particles = file['blocks/particles/data'][:]
        cpml_nonzero = 0
        if 'cpml' in file['blocks']:
            cpml_nonzero = sum(np.count_nonzero(block['data'][:])
                               for block in file['blocks/cpml'].values())
        return {'particles': len(particles),
                'held': int(np.count_nonzero(particles['frozen'])),
                'cpml_nonzero': int(cpml_nonzero),
                'cpml_carriers': len(file['blocks/carriers/cpml/data']),
                'retirement_carriers': len(file['blocks/carriers/retirement/data'])}


def gaussian_card():
    return {
        'units': {'length': 'micrometer', 'time': 'femtosecond'},
        'runtime': {'mode': 'throughput', 'resource_monitor': {'enabled': False}},
        'output': {'overwrite': False, 'manifest': 'run-manifest.yaml'},
        'mesh': {'field_solver': 'cowan-z', 'current_filter': 'none',
                 'cells': [16, 16, 256], 'cell_size': [1, 1, .5],
                 'center': [0, 0, 0], 'duration': 600, 'boost_gamma': 2,
                 'particle_steps_per_undulator_period': 4,
                 'self_field_feedback': True, 'ballistic_outside_transverse': False},
        'boundary': {'type': 'cpml', 'cells': [2, 2, 4], 'polynomial_order': 3,
                     'target_reflection': 1e-8, 'kappa_max': 8, 'alpha_fraction': 0},
        'initial_self_field': {'enabled': True, 'model': 'relativistic-poisson',
                              'relative_tolerance': 1e-8, 'maximum_iterations': 2000},
        'beam': {'reference': {'initial_center_z': 1}, 'input': {
            'type': 'generated-gaussian', 'electrons': 64, 'macroparticles': 64,
            'gamma': 4, 'direction': [0, 0, 1], 'center': [0, 0, 0],
            'sigma_position': [.1, .1, .1], 'sigma_proper_velocity': [0, 0, 0],
            'random_seed': 261026}},
        'sources': {'incident_waves': [], 'magnetic_elements': []},
        'detectors': {'enabled': True, 'directory': 'detectors',
                      'field_planes': [{'name': 'field', 'z': 60, 'rhythm': 2.3,
                                        'buffer_samples': 4, 'compression': 0}],
                      'particle_planes': [{'name': 'entry', 'z': 18,
                                           'buffer_records': 16, 'compression': 0}]},
        'stop': {'mode': 'reference-center-z', 'z': 100},
        'trajectory': {'enabled': True, 'directory': 'trajectory',
                       'basename': 'trajectory', 'rhythm': 6.7, 'buffer_records': 64,
                       'compression': 0},
        'energy_ledger': {'enabled': True, 'directory': 'ledger',
                          'filename': 'energy-ledger.h5', 'sample_interval_steps': 5,
                          'buffer_records': 8, 'compression': 0,
                          'warning_relative_tolerance': .05}}


def held_input(path):
    dtype = np.dtype([('plane_position_m', '<f8', (2,)),
                      ('arrival_time_offset_s', '<f8'),
                      ('proper_velocity', '<f8', (3,)),
                      ('source_id', '<u8'), ('macro_weight', '<f8')])
    records = np.zeros(64, dtype=dtype)
    rng = np.random.default_rng(261026)
    records['plane_position_m'] = rng.normal(0, .05e-6, (64, 2))
    records['arrival_time_offset_s'] = np.linspace(-2e-15, 2e-15, 64)
    records['proper_velocity'][:, 2] = math.sqrt(4 * 4 - 1)
    records['source_id'] = np.arange(1, 65)
    records['macro_weight'] = 1
    with h5py.File(path, 'x') as file:
        group = file.create_group('particles')
        group.attrs['format_version'] = 4
        group.create_dataset('records', data=records)


def cases(root):
    simple = gaussian_card()
    held = copy.deepcopy(simple)
    held['beam'] = {'reference': {'initial_center_z': 1, 'input_plane_z': 0,
                                 'frozen_lab_plane_injection': True,
                                 'activation_plane_z': 0,
                                 'hold_transverse_until_activation': True},
                    'input': {'type': 'hdf5', 'file': str(root / 'inputs/held.h5'),
                              'electrons': 64, 'position_offset': [0, 0, 0]}}
    held['particle_retirement'] = {'enabled': True, 'entrance_z': 20, 'length': 10}
    reference = yaml.safe_load((HERE.parent / 'cards/postprocess-ballistic-run.yaml').read_text())
    reference['output']['overwrite'] = False
    reference['detectors']['field_planes'][0]['particle_background']['validation']['enabled'] = True
    reference['stop'] = {'mode': 'reference-center-z', 'z': 2}
    reference['energy_ledger'] = copy.deepcopy(simple['energy_ledger'])
    reference['trajectory'] = copy.deepcopy(simple['trajectory'])
    return {'gaussian-cpml': (simple, [0, 7, 17, 55], [0, 13.7]),
            'held-retirement': (held, [1, 7, 17, 55], []),
            'ballistic-validation': (reference, [1, 7, 17], [])}


def write_card(directory, base, checkpoint=None):
    directory.mkdir(parents=True, exist_ok=False)
    card = copy.deepcopy(base)
    card['output'] = {'overwrite': False, 'manifest': str(directory / 'run-manifest.yaml')}
    for section in ('trajectory', 'detectors', 'energy_ledger'):
        if section in card:
            target = 'ledger' if section == 'energy_ledger' else section
            card[section]['directory'] = str(directory / target)
    if checkpoint is not None:
        card['checkpoint'] = copy.deepcopy(checkpoint)
    (directory / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
    return card


def prepare(args):
    root = args.root.resolve()
    require(not root.exists(), 'Use a fresh result directory')
    for executable in (args.binary, args.baseline):
        require(executable.is_absolute() and executable.is_file() and os.access(executable, os.X_OK),
                'Supply absolute executable paths')
    launcher = json.loads(args.launcher_json)
    require(isinstance(launcher, list) and all(isinstance(x, str) for x in launcher),
            '--launcher-json must be a JSON argument array')
    require(any('{ranks}' in token for token in launcher), 'Launcher must include {ranks}')
    require(args.ranks in (1, 2, 4, 8) and args.threads > 0 and
            args.ranks * args.threads <= 16, 'Small test budget is at most 16 CPUs')
    require(args.timeout > 0, 'Per-run timeout must be positive')
    root.mkdir(parents=True)
    (root / 'inputs').mkdir()
    held_input(root / 'inputs/held.h5')
    plan = {'format': 'aprl-checkpoint-acceptance-v1', 'binary': str(args.binary),
            'baseline': str(args.baseline), 'binary_sha256': digest(args.binary),
            'baseline_sha256': digest(args.baseline), 'launcher': launcher,
            'ranks': args.ranks, 'threads': args.threads, 'timeout_seconds': args.timeout,
            'cases': {}, 'scope': 'Small same-rank correctness, not full-scale I/O throughput'}
    for name, (card, steps, times) in cases(root).items():
        checkpoint = {'directory': str(root / name / 'continuous/checkpoints'),
                      'steps': steps, 'times': times, 'io_timeout_seconds': 60}
        write_card(root / name / 'continuous', card, checkpoint)
        write_card(root / name / 'disabled', card)
        write_card(root / name / 'accepted', card)
        plan['cases'][name] = {'steps': steps, 'times': times, 'base': card}
    save(root / 'plan.json', plan)
    save(root / 'frozen.json', {str(p): digest(p) for p in sorted(root.rglob('*')) if p.is_file()})
    print(root / 'plan.json')


def launch(root, plan, directory, binary=None, ranks=None, rejection=None):
    binary = binary or plan['binary']
    ranks = plan['ranks'] if ranks is None else ranks
    command = [token.format(ranks=ranks, threads=plan['threads']) for token in plan['launcher']]
    command += [binary, str(directory / 'run.yaml')]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(plan['threads']), OMP_DYNAMIC='FALSE',
               OMP_PLACES='cores', OMP_PROC_BIND='close', APRL_PROFILE_INTERVAL='0',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    libraries = Path(binary).parent / 'runtime-libs'
    if libraries.is_dir():
        env['LD_LIBRARY_PATH'] = str(libraries) + ':' + env.get('LD_LIBRARY_PATH', '')
    start = time.monotonic()
    with (directory / 'aprl.log').open('w') as log:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=plan['timeout_seconds'])
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise RuntimeError('Run timed out; cannot count timeout as a passed rejection: ' + str(directory))
    text = (directory / 'aprl.log').read_text(errors='replace')
    result = {'command': command, 'exit_code': code, 'seconds': time.monotonic() - start,
              'expected_rejection': rejection}
    save(directory / 'launch.json', result)
    if rejection:
        require(code not in (0, 124, 137, 143) and re.search(rejection, text, re.I),
                'Expected explicit rejection was not observed: ' + str(directory))
        require('Configured stop reached:' not in text, 'Rejected run reached physical completion')
    else:
        require(code == 0 and 'Configured stop reached:' in text,
                'Small fixture must reach its physical stop: ' + str(directory))
        for path in output_files(directory).values():
            inspect_file(path)
    return result


def checkpoint_path(directory, step):
    return directory / 'checkpoints' / ('step-%012d.h5' % step)


def reject(root, plan, name, card, checkpoint, changes=None, ranks=None):
    directory = root / 'rejections' / name
    modified = copy.deepcopy(card)
    if changes:
        changes(modified)
    write_card(directory, modified, checkpoint)
    return launch(root, plan, directory, ranks=ranks,
                  rejection=r'checkpoint|restart|retire|schema|mesh|configuration')


def run(args):
    root = args.root.resolve()
    require(not (root / 'RUN_STARTED.json').exists(), 'Campaign already started; use a fresh directory')
    host = socket.gethostname().split('.')[0]
    require(os.environ.get('SLURM_JOB_ID') or (host == 'gnode1' and args.allow_gnode1),
            'Run inside a real Slurm allocation or explicitly authorized gnode1 session')
    for name, expected in json.loads((root / 'frozen.json').read_text()).items():
        require(digest(Path(name)) == expected, 'Prepared input changed: ' + name)
    plan = json.loads((root / 'plan.json').read_text())
    require(digest(Path(plan['binary'])) == plan['binary_sha256'], 'Candidate binary changed')
    require(digest(Path(plan['baseline'])) == plan['baseline_sha256'], 'Accepted binary changed')
    save(root / 'RUN_STARTED.json', {'host': host, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                                    'started_unix': time.time()})
    comparisons, timings, snapshots = [], [], {}
    for name, item in plan['cases'].items():
        case_root = root / name
        for variant in ('accepted', 'disabled', 'continuous'):
            launch(root, plan, case_root / variant,
                   binary=plan['baseline'] if variant == 'accepted' else plan['binary'])
        for first, second in [('accepted', 'disabled'), ('disabled', 'continuous')]:
            compared = compare_outputs(case_root / first, case_root / second)
            comparisons.append({'case': name, 'pair': [first, second], 'datasets': compared})
        for step in item['steps']:
            path = checkpoint_path(case_root / 'continuous', step)
            require(path.is_file(), 'Scheduled snapshot missing: ' + str(path))
            snapshots[name + ':' + str(step)] = inspect_snapshot(path)
            directory = case_root / ('restart-%d' % step)
            settings = {'restart_from': str(path), 'reset_diagnostics': False,
                        'io_timeout_seconds': 60}
            if name == 'gaussian-cpml' and step == 7:
                settings.update(directory=str(directory / 'checkpoints'), steps=[17])
            write_card(directory, item['base'], settings)
            original_input = root / 'inputs/held.h5'
            hidden_input = root / 'inputs/held.h5.hidden-during-restart'
            if name == 'held-retirement':
                original_input.rename(hidden_input)
            try:
                launch(root, plan, directory)
            finally:
                if hidden_input.exists():
                    hidden_input.rename(original_input)
            compared = compare_outputs(case_root / 'continuous', directory, step, path)
            comparisons.append({'case': name, 'restart_step': step, 'datasets': compared})
        for line in (case_root / 'continuous/aprl.log').read_text().splitlines():
            if '[checkpoint]' in line and 'phase=write' in line:
                timings.append({'case': name, 'line': line})
        if item['times']:
            # Cowan's time step is fixed by dz. The requested 13.7fs is
            # deliberately between grid times, away from rounding boundaries.
            dt_fs = item['base']['mesh']['cell_size'][2] * 1e-6 / C / 1e-15
            for requested in item['times']:
                actual_step = int(math.ceil(requested / dt_fs))
                time_snapshot = checkpoint_path(case_root / 'continuous', actual_step)
                require(time_snapshot.is_file(), 'Time-triggered snapshot missing')
                _, saved_time = snapshot_metadata(time_snapshot)
                require(saved_time >= requested * 1e-15 and
                        saved_time < requested * 1e-15 + dt_fs * 1e-15,
                        'Snapshot time trigger did not use the first complete field step')
    require(timings, 'No checkpoint write timing was reported')
    require(snapshots['gaussian-cpml:17']['cpml_nonzero'] > 0,
            'CPML snapshot did not exercise nonzero convolution history')
    require(snapshots['held-retirement:1']['held'] > 0,
            'Held fixture has no frozen particles at its early checkpoint')
    require(snapshots['held-retirement:17']['retirement_carriers'] > 0,
            'Held fixture has no in-flight retirement-current carriers at step17')
    require(snapshots['held-retirement:55']['particles'] == 0 and
            snapshots['held-retirement:55']['retirement_carriers'] == 0,
            'Retirement fixture was not fully empty before step55')
    validation_snapshot = checkpoint_path(root / 'ballistic-validation/continuous', 7)
    validation_prefixes, _ = snapshot_metadata(validation_snapshot)
    require(any('ballistic-reference' in filename and
                prefix['records'] > prefix.get('validation_records', 0)
                for filename, prefix in validation_prefixes.items()),
            'Ballistic fixture did not straddle reference and validation crossings')
    # The new archive's diagnostic prefixes belong to the first restart segment,
    # so compare the second restart against that immediate parent segment.
    second_parent = root / 'gaussian-cpml/restart-7'
    second_snapshot = checkpoint_path(second_parent, 17)
    second_restart = root / 'gaussian-cpml/second-generation-17'
    write_card(second_restart, plan['cases']['gaussian-cpml']['base'],
               {'restart_from': str(second_snapshot), 'io_timeout_seconds': 60})
    launch(root, plan, second_restart)
    comparisons.append({'case': 'gaussian-cpml', 'generation': 2,
        'datasets': compare_outputs(second_parent, second_restart, 17, second_snapshot)})
    base = plan['cases']['gaussian-cpml']['base']
    source = checkpoint_path(root / 'gaussian-cpml/continuous', 7)
    negative = root / 'negative-inputs'
    negative.mkdir()
    partial = negative / 'unfinished.h5'
    # Deliberately generate a valid HDF5 container with an incomplete archive
    # header, not a random file whose failure could hide a schema-check bug.
    with h5py.File(source) as original, h5py.File(partial, 'x') as broken:
        header = original['header'][:]
        header[2] = 0
        broken.create_dataset('header', data=header)
    reject(root, plan, 'partial', base, {'restart_from': str(partial)})
    reject(root, plan, 'different-mpi', base, {'restart_from': str(source)},
           ranks=1 if plan['ranks'] != 1 else 2)
    reject(root, plan, 'different-mesh', base, {'restart_from': str(source)},
           changes=lambda card: card['mesh']['cell_size'].__setitem__(0, 1.25))
    retired = checkpoint_path(root / 'held-retirement/continuous', 55)
    retirement_base = plan['cases']['held-retirement']['base']
    reject(root, plan, 'retired-layer-change', retirement_base,
           {'restart_from': str(retired), 'allow_retirement_change': True,
            'reset_diagnostics': True},
           changes=lambda card: card['particle_retirement'].__setitem__('length', 11))
    for name, invalid in [('negative-step', {'steps': [-1]}),
                          ('negative-time', {'times': [-1]}),
                          ('infinite-time', {'times': [float('inf')]}),
                          ('zero-timeout', {'io_timeout_seconds': 0}),
                          ('empty-path', {'restart_from': ''})]:
        reject(root, plan, name, base, invalid)
    ledger_file = root / 'held-retirement/continuous/ledger/energy-ledger.h5'
    _, arrays, _ = inspect_file(ledger_file)
    records = arrays['records']
    before = records[records['step'] <= 55]
    require(len(before) and before[-1]['active_macroparticles'] == 0 and
            before[-1]['removed_macroparticles'] == 64,
            'Retirement fixture did not clear all 64 particles before step55')
    result = {'status': 'PASS', 'comparisons': comparisons, 'write_timings': timings,
              'snapshot_coverage': snapshots,
              'limitations': ['Small deterministic same-rank cases only',
                              'Scientific post-snapshot suffixes; prior output remains in parent',
                              'No full-scale filesystem throughput or normal-loop speed claim']}
    save(root / 'RESULT.json', result)
    (root / 'PASS').write_text('All checkpoint correctness gates passed.\n')
    print(root / 'RESULT.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prepare_parser = sub.add_parser('prepare')
    prepare_parser.add_argument('--root', type=Path, required=True)
    prepare_parser.add_argument('--binary', type=Path, required=True)
    prepare_parser.add_argument('--baseline', type=Path, required=True)
    prepare_parser.add_argument('--ranks', type=int, default=4)
    prepare_parser.add_argument('--threads', type=int, default=1)
    prepare_parser.add_argument('--timeout', type=int, default=180)
    prepare_parser.add_argument('--launcher-json', required=True)
    run_parser = sub.add_parser('run')
    run_parser.add_argument('--root', type=Path, required=True)
    run_parser.add_argument('--allow-gnode1', action='store_true')
    args = parser.parse_args()
    try:
        prepare(args) if args.action == 'prepare' else run(args)
    except Exception as error:
        if args.root.exists() and args.action == 'run':
            save(args.root / 'FAIL.json', {'error': str(error), 'time_unix': time.time()})
        raise


if __name__ == '__main__':
    main()
