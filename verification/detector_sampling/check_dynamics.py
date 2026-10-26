#!/usr/bin/env python3
"""Small same-layout detector-fix dynamics/restart regression; never submits jobs.

Prepare is read/write setup only. Run requires Slurm or explicitly authorized
gnode1, performs three bounded four-rank jobs, and checks actual physical stops.
"""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import struct
import time

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    'checkpoint_acceptance_common', HERE.parent / 'checkpoint_restart/run.py')
common = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(common)


def fixture():
    card = common.gaussian_card()
    card['mesh'].update(cells=[12, 12, 80], duration=150)
    card['detectors']['field_planes'][0].update(z=28, rhythm=.01)
    card['stop']['z'] = 36
    card['energy_ledger']['sample_interval_steps'] = 1
    return card


def prepare(args):
    root = args.root.resolve()
    common.require(not root.exists(), 'Use a fresh result directory')
    for binary in (args.old, args.new):
        common.require(binary.is_absolute() and binary.is_file() and
                       os.access(binary, os.X_OK), 'Supply absolute executable paths')
    launcher = json.loads(args.launcher_json)
    common.require(isinstance(launcher, list) and launcher and
                   all(isinstance(item, str) for item in launcher) and
                   any('{ranks}' in item for item in launcher),
                   'Launcher must be a JSON argument array containing {ranks}')
    root.mkdir(parents=True)
    card = fixture()
    for name in ('old', 'continuous'):
        directory = root / name
        common.write_card(directory, card, {
            'directory': str(directory / 'checkpoints'), 'steps': [7, 30],
            'io_timeout_seconds': 60})
    plan = {
        'format': 'aprl-detector-dynamics-v1', 'baseline': str(args.old),
        'binary': str(args.new), 'baseline_sha256': common.digest(args.old),
        'binary_sha256': common.digest(args.new), 'launcher': launcher,
        'ranks': 4, 'threads': 1, 'timeout_seconds': 180, 'base': card,
        'scope': 'Same-model same-layout correctness; no performance claim',
    }
    common.save(root / 'plan.json', plan)
    common.save(root / 'frozen.json', {
        str(path): common.digest(path) for path in sorted(root.rglob('*'))
        if path.is_file()})
    print(root / 'plan.json')


def compare_nonfield_outputs(first, second):
    a, b = common.output_files(first), common.output_files(second)
    common.require(a.keys() == b.keys(), 'Scientific output file sets differ')
    checked = []
    for name in sorted(a):
        ga, da, aa = common.inspect_file(a[name])
        gb, db, ab = common.inspect_file(b[name])
        common.require(ga == gb and da.keys() == db.keys(), name + ': output layout differs')
        # The detector fields intentionally change; their native timing does not.
        if 'fields' in da:
            common.bitwise_equal(da['time_s'], db['time_s'], name + '/time_s')
            checked.append(name + '/time_s')
            continue
        for attr in set(aa) | set(ab):
            if attr in common.IDENTITY_ATTRIBUTES:
                continue
            common.require(attr in aa and attr in ab, name + ': missing attribute ' + attr)
            common.bitwise_equal(aa[attr], ab[attr], name + '@' + attr)
        for dataset in da:
            common.bitwise_equal(da[dataset], db[dataset], name + '/' + dataset)
            checked.append(name + '/' + dataset)
    common.require(any(name.startswith('ledger/') for name in checked), 'No ledger comparison')
    common.require(any(name.startswith('trajectory/') for name in checked), 'No trajectory comparison')
    return checked


def normalized_scalars(file):
    """Exclude only six x87 padding bytes from each of 19 long-double values.

    transferCheckpointScalars schema v1 packs scalars individually. This prefix
    is uint32 version, double time, size_t step, double trajectory time, uint32
    trajectory flush count, then 19 native long doubles. Every meaningful byte,
    including all later ledger/counter state, remains subject to exact equality.
    Refuse other ABIs instead of guessing their padding layout.
    """
    header = file['header'][:]
    common.require(tuple(int(header[i]) for i in (4, 5, 6, 7, 8, 9)) ==
                   (8, 8, 16, 64, 16384, 1), 'Unsupported scalar ABI for padding-aware check')
    block = file['blocks/state/scalars']
    result = block['data'][:].copy()
    offset = 0
    for count in block['rank_counts'][:]:
        count = int(count)
        common.require(count >= 32 + 19 * 16 and
                       struct.unpack_from('<I', result, offset)[0] == 1,
                       'Unsupported or truncated scalar state schema')
        for index in range(19):
            first_padding = offset + 32 + index * 16 + 10
            result[first_padding:first_padding + 6] = 0
        offset += count
    common.require(offset == len(result), 'Scalar rank counts do not cover data')
    return result


def compare_snapshots(first, second):
    checked = []
    with h5py.File(first) as a, h5py.File(second) as b:
        common.bitwise_equal(a['header'][:], b['header'][:], 'checkpoint/header')
        common.require(int(a['header'][2]) == 1, 'Checkpoint incomplete')
        paths = []
        a['blocks'].visititems(lambda name, obj: paths.append(name)
                              if isinstance(obj, h5py.Dataset) and
                              not name.startswith('metadata/') else None)
        other_paths = []
        b['blocks'].visititems(lambda name, obj: other_paths.append(name)
                              if isinstance(obj, h5py.Dataset) and
                              not name.startswith('metadata/') else None)
        common.require(paths == other_paths, 'Checkpoint scientific block layout differs')
        for name in paths:
            if name == 'state/scalars/data':
                da, db = normalized_scalars(a), normalized_scalars(b)
            else:
                da, db = a['blocks/' + name][()], b['blocks/' + name][()]
            common.bitwise_equal(da, db, name)
            checked.append(name)
        common.require(any(name.startswith('fields/') for name in checked) and
                       'particles/data' in checked and 'state/scalars/data' in checked,
                       'Missing expected scientific checkpoint state')
    return checked


def detector_coverage(directory, card, snapshot):
    with h5py.File(snapshot) as file:
        block = file['blocks/state/geometry']
        data = block['data'][:int(block['rank_counts'][0])].tobytes()
        common.require(struct.unpack_from('<I', data)[0] == 1,
                       'Unsupported saved geometry version')
        # Version + three size_t cell counts + three double extents.
        origin_time, origin_z, origin_box_z = struct.unpack_from('<3d', data, 4 + 24 + 24)
    _, arrays, _ = common.inspect_file(directory / 'detectors/field.h5')
    gamma = card['mesh']['boost_gamma']
    beta = math.sqrt(1 - 1 / gamma ** 2)
    plane = card['detectors']['field_planes'][0]['z'] * 1e-6
    time_box = gamma * (arrays['time_s'] - origin_time - beta * (plane - origin_z) / common.C)
    dz = card['mesh']['cell_size'][2] * 1e-6
    nz = card['mesh']['cells'][2]
    box_z = origin_box_z + (plane - origin_z) / gamma - beta * common.C * time_box
    origin_grid = card['mesh']['center'][2] * 1e-6 - .5 * nz * dz
    coordinate = (box_z - origin_grid) / dz
    steps = np.rint(time_box / (dz / common.C)).astype(int)
    common.require(np.array_equal(steps, np.arange(len(steps))),
                   'Fixture must sample every complete field step including initialization')
    owners = np.minimum(np.maximum(np.floor(coordinate).astype(int), 0), nz - 1) // (nz // 4)
    changes = np.flatnonzero(owners[1:] != owners[:-1]) + 1
    common.require(len(changes) >= 2 and any(index > 7 for index in changes),
                   'Detector did not cross two rank boundaries with one after restart')
    common.require(any(abs(int(index) - 7) <= 4 for index in changes),
                   'Restart checkpoint is not near a detector ownership transition')
    common.require(steps[-1] >= 30 and np.isfinite(coordinate).all(), 'Insufficient dynamics coverage')
    return {'samples': len(steps), 'last_step': int(steps[-1]),
            'owner_transition_steps': changes.tolist(),
            'transition_coordinates': coordinate[changes].tolist(),
            'checkpoint_steps': [7, 30], 'restart_step': 7}


def field_change(first, second):
    _, old, _ = common.inspect_file(first / 'detectors/field.h5')
    _, new, _ = common.inspect_file(second / 'detectors/field.h5')
    fields_old, fields_new = old['fields'], new['fields']
    result = {}
    changed_frames = set()
    for group, prefix in [('electric_V_per_m', 'E'), ('magnetic_T', 'B')]:
        for component, axis in enumerate('xyz'):
            delta = np.abs(fields_new[group][..., component] - fields_old[group][..., component])
            result[prefix + axis] = float(np.max(delta))
            changed_frames.update(np.flatnonzero(np.any(delta != 0, axis=(1, 2))).tolist())
    common.require(any(value > 0 for value in result.values()),
                   'Fixture did not exercise any corrected detector samples')
    return {'maximum_absolute_difference_by_component_native_SI': result,
            'changed_sample_indices': sorted(changed_frames),
            'interpretation': 'Diagnostic output intentionally changes; no physical spectrum claim'}


def run(args):
    root = args.root.resolve()
    common.require(not (root / 'RUN_STARTED.json').exists(), 'Use a fresh prepared directory')
    host = socket.gethostname().split('.')[0]
    common.require(os.environ.get('SLURM_JOB_ID') or
                   (host == 'gnode1' and args.allow_gnode1),
                   'Run in Slurm or the explicitly authorized gnode1 session')
    for path, expected in json.loads((root / 'frozen.json').read_text()).items():
        common.require(common.digest(Path(path)) == expected, 'Prepared input changed: ' + path)
    plan = json.loads((root / 'plan.json').read_text())
    for key, hash_key in [('baseline', 'baseline_sha256'), ('binary', 'binary_sha256')]:
        common.require(common.digest(Path(plan[key])) == plan[hash_key], 'Binary changed: ' + key)
    common.save(root / 'RUN_STARTED.json', {'host': host, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                                           'started_unix': time.time()})
    launches = {}
    for name, binary in [('old', plan['baseline']), ('continuous', plan['binary'])]:
        launches[name] = common.launch(root, plan, root / name, binary=binary)
    snapshots = {}
    for step in (7, 30):
        snapshots[str(step)] = compare_snapshots(
            common.checkpoint_path(root / 'old', step),
            common.checkpoint_path(root / 'continuous', step))
    nonfield = compare_nonfield_outputs(root / 'old', root / 'continuous')
    parent = common.checkpoint_path(root / 'continuous', 7)
    coverage = detector_coverage(root / 'continuous', plan['base'], parent)
    changes = field_change(root / 'old', root / 'continuous')
    restart = root / 'restart-7'
    common.write_card(restart, plan['base'], {
        'restart_from': str(parent), 'directory': str(restart / 'checkpoints'),
        'steps': [30], 'reset_diagnostics': False, 'io_timeout_seconds': 60})
    launches['restart-7'] = common.launch(root, plan, restart)
    restart_outputs = common.compare_outputs(root / 'continuous', restart, 7, parent)
    restart_state = compare_snapshots(common.checkpoint_path(root / 'continuous', 30),
                                      common.checkpoint_path(restart, 30))
    result = {'status': 'PASS', 'launches': launches, 'same_step_state': snapshots,
              'nonfield_datasets_exact': nonfield, 'coverage': coverage,
              'intentional_field_change': changes, 'restart_outputs_exact': restart_outputs,
              'restart_step30_state_exact': restart_state,
              'excluded_checkpoint_data': ['metadata (identity, configuration paths, output prefixes)',
                                            'six padding bytes per x87 long double; all meaningful bits compared'],
              'limitations': ['Small generated-gaussian same-four-rank correctness test',
                              'No full-scale throughput or pure-radiation validation']}
    common.save(root / 'RESULT.json', result)
    (root / 'PASS').write_text('Dynamics unchanged and detector continuation exact.\n')
    print(root / 'RESULT.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    setup = sub.add_parser('prepare')
    setup.add_argument('--root', type=Path, required=True)
    setup.add_argument('--old', type=Path, required=True)
    setup.add_argument('--new', type=Path, required=True)
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
