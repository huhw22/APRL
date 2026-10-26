#!/usr/bin/env python3
"""Bounded held field-only acceptance; no submission and no core changes.

The two four-rank runs differ only in particle-plane output. Both stop after
all surviving particles cross the final field plane, without a paired-plane
requirement. This is a small diagnostic regression, not radiation validation.
"""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import re
import socket
import time
import h5py

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    'field_only_checkpoint_helpers', HERE.parents[2] / 'checkpoint_restart/run.py')
common = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(common)


def prepare(args):
    root = args.root.resolve()
    common.require(not root.exists(), 'Use a fresh result directory')
    binary = args.binary.resolve()
    common.require(binary.is_file() and os.access(binary, os.X_OK),
                   'Supply an existing executable')
    launcher = json.loads(args.launcher_json)
    common.require(isinstance(launcher, list) and launcher and
                   all(isinstance(item, str) for item in launcher) and
                   any('{ranks}' in item for item in launcher),
                   'Launcher must be a JSON argument array containing {ranks}')
    root.mkdir(parents=True)
    (root / 'inputs').mkdir()
    common.held_input(root / 'inputs/held.h5')
    card = common.gaussian_card()
    card['beam'] = {'reference': {
        'initial_center_z': 1, 'input_plane_z': 0,
        'frozen_lab_plane_injection': True, 'activation_plane_z': 0,
        'hold_transverse_until_activation': True},
        'input': {'type': 'hdf5', 'file': str(root / 'inputs/held.h5'),
                  'electrons': 64, 'position_offset': [0, 0, 0]}}
    card['mesh']['duration'] = 300
    card['sources']['magnetic_elements'] = [{
        'type': 'uniform-dipole', 'entrance_z': 18, 'length': 2,
        'field_T': .1, 'polarization_angle_rad': .3}]
    card['particle_retirement'] = {'enabled': False}
    card['trajectory'] = {'enabled': False}
    card['stop'] = {'mode': 'after-last-element'}
    field = card['detectors']['field_planes'][0]
    field.update(z=60, rhythm=.01, particle_background={'enabled': False})
    card['energy_ledger']['sample_interval_steps'] = 1
    card['detectors']['particle_planes'] = [{
        'name': 'paired-particles', 'z': 60, 'buffer_records': 16,
        'compression': 0}]
    common.write_card(root / 'paired', card, {
        'directory': str(root / 'paired/checkpoints'),
        'steps': [10], 'io_timeout_seconds': 60})
    only = copy.deepcopy(card)
    only['detectors']['particle_planes'] = []
    common.write_card(root / 'field-only', only, {
        'directory': str(root / 'field-only/checkpoints'),
        'steps': [10], 'io_timeout_seconds': 60})
    plan = {
        'format': 'aprl-held-field-only-smoke-v1', 'binary': str(binary),
        'binary_sha256': common.digest(binary), 'launcher': launcher,
        'ranks': 4, 'threads': 1, 'timeout_seconds': 180,
        'expected_input_particles': 64,
        'scope': 'Identical held model; particle output present versus absent',
        'physical_stop': 'after-last-element; final field plane z=60 micrometres',
    }
    common.save(root / 'plan.json', plan)
    common.save(root / 'frozen.json', {
        str(path): common.digest(path) for path in sorted(root.rglob('*'))
        if path.is_file()})
    print(root / 'plan.json')


def compare(first, second):
    ga, a, _ = common.inspect_file(first)
    gb, b, _ = common.inspect_file(second)
    common.require(ga == gb and a.keys() == b.keys(),
                   'Output groups/datasets differ: ' + str(second))
    for name in a:
        common.bitwise_equal(a[name], b[name], str(second) + '/' + name)
    return list(a)


def run(args):
    root = args.root.resolve()
    common.require(not (root / 'RUN_STARTED.json').exists(),
                   'Smoke already started; use a fresh result directory')
    host = socket.gethostname().split('.')[0]
    common.require(os.environ.get('SLURM_JOB_ID') or
                   (host == 'gnode1' and args.allow_gnode1 and
                    os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1'),
                   'Requires a real Slurm allocation or authorized gnode1')
    for path, expected in json.loads((root / 'frozen.json').read_text()).items():
        common.require(common.digest(Path(path)) == expected,
                       'Prepared input changed: ' + path)
    plan = json.loads((root / 'plan.json').read_text())
    common.require(common.digest(Path(plan['binary'])) == plan['binary_sha256'],
                   'Executable changed after preparation')
    common.save(root / 'RUN_STARTED.json', {
        'host': host, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
        'started_unix': time.time()})
    launches, coverage, stops, snapshots = {}, {}, [], {}
    for name in ('paired', 'field-only'):
        directory = root / name
        launches[name] = common.launch(root, plan, directory)
        text = (directory / 'aprl.log').read_text()
        held = re.search(r'Frozen lab-plane injection: activations=(\d+), '
                         r'cross-rank activation pushes=(\d+), still frozen=(\d+)', text)
        common.require(held and int(held[1]) == 64 and int(held[3]) == 0,
                       'Small run did not activate all held particles: ' + name)
        common.require('held-transverse pseudo-trajectory' in text,
                       'Held transverse startup missing: ' + name)
        matches = re.findall(r'Configured stop reached:.*', text)
        common.require(len(matches) == 1, 'Missing unique physical stop: ' + name)
        stops.append(matches[0])
        expected = 1 if name == 'paired' else 0
        common.require(f'Laboratory detector planes active: fields=1, particles={expected};' in text,
                       'Unexpected detector counts: ' + name)
        names = {p.name for p in (directory / 'detectors').glob('*.h5')}
        common.require(names == ({'field.h5', 'paired-particles.h5'} if expected else {'field.h5'}),
                       'Unexpected field-only or paired output file set: ' + name)
        _, fields, _ = common.inspect_file(directory / 'detectors/field.h5')
        common.require(len(fields['time_s']) > 10, 'Too few field samples')
        coverage[name] = {
            'held_activated': int(held[1]), 'cross_rank_activation_pushes': int(held[2]),
            'still_held': int(held[3]), 'field_samples': len(fields['time_s']),
            'detector_files': sorted(names)}
        checkpoint = common.checkpoint_path(directory, 10)
        common.require(list((directory / 'checkpoints').glob('step-*.h5')) == [checkpoint],
                       'Expected exactly one complete step-10 checkpoint')
        with h5py.File(checkpoint) as file:
            common.require(int(file['header'][2]) == 1,
                           'Checkpoint was not committed')
        snapshots[name] = common.inspect_snapshot(checkpoint)
        common.require(snapshots[name]['particles'] == 64 and
                       snapshots[name]['retirement_carriers'] == 0,
                       'Full snapshot must retain all particles without retirement')
    common.require(stops[0] == stops[1], 'Removing particle output changed physical stop')
    checked = {}
    for path in ('detectors/field.h5', 'ledger/energy-ledger.h5'):
        checked[path] = compare(root / 'paired' / path, root / 'field-only' / path)
    snapshot_checked = []
    with h5py.File(common.checkpoint_path(root / 'paired', 10)) as first, \
            h5py.File(common.checkpoint_path(root / 'field-only', 10)) as second:
        for category in ('fields', 'particles', 'cpml', 'carriers'):
            paths = []
            first['blocks/' + category].visititems(
                lambda name, obj: paths.append(name) if isinstance(obj, h5py.Dataset) else None)
            for path in paths:
                label = 'blocks/' + category + '/' + path
                common.bitwise_equal(first[label][()], second[label][()], label)
                snapshot_checked.append(label)
    _, particles, _ = common.inspect_file(root / 'paired/detectors/paired-particles.h5')
    common.require(len(particles['records']) == 64,
                   'Paired control must record all 64 particles at final field plane')
    common.save(root / 'RESULT.json', {
        'status': 'PASS', 'launches': launches, 'coverage': coverage,
        'physical_stop': stops[0], 'exact_datasets': checked,
        'snapshots': snapshots, 'exact_snapshot_datasets': snapshot_checked,
        'paired_records': len(particles['records']),
        'limitations': ['Small same-four-rank held diagnostic test',
                        'No full-scale runtime or radiation-convergence claim']})
    (root / 'PASS').write_text('Held field-only outputs and physical stop independent of particle planes.\n')
    print(root / 'RESULT.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    setup = sub.add_parser('prepare')
    setup.add_argument('--root', type=Path, required=True)
    setup.add_argument('--binary', type=Path, required=True)
    setup.add_argument('--launcher-json', required=True)
    execute = sub.add_parser('run')
    execute.add_argument('--root', type=Path, required=True)
    execute.add_argument('--allow-gnode1', action='store_true')
    args = parser.parse_args()
    try:
        prepare(args) if args.action == 'prepare' else run(args)
    except Exception as error:
        if args.action == 'run' and args.root.is_dir():
            common.save(args.root / 'FAIL.json', {'error': str(error), 'time_unix': time.time()})
        raise


if __name__ == '__main__':
    main()
