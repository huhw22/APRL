#!/usr/bin/env python3
"""Bounded two-node checkpoint I/O, with exact continued-state verification."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

import h5py
import numpy as np
import yaml


def save(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def equal(a, b, name):
    assert a.shape == b.shape and a.dtype == b.dtype, name
    if a.dtype.names:
        for field in a.dtype.names:
            equal(a[field], b[field], name + '.' + field)
    else:
        assert a.tobytes() == b.tobytes(), name


def prepare(args):
    root = args.root.resolve()
    assert not root.exists() and args.binary.is_absolute() and args.binary.is_file()
    root.mkdir(parents=True)
    (root / 'scripts').mkdir()
    for path in Path(__file__).parent.iterdir():
        if path.is_file():
            (root / 'scripts' / path.name).write_bytes(path.read_bytes())
    step_fs = .5e-6 / 299792458.0 / 1e-15
    card = {
        'units': {'length': 'micrometer', 'time': 'femtosecond'},
        'runtime': {'mode': 'throughput', 'resource_monitor': {
            'enabled': True, 'calibration_steps': 1, 'progress_interval_steps': 3}},
        'output': {'overwrite': False, 'manifest': 'run-manifest.yaml'},
        'mesh': {'field_solver': 'cowan-z', 'cells': [128, 128, 2048],
                 'cell_size': [1, 1, .5], 'center': [0, 0, 0],
                 'duration': 11.5 * step_fs, 'boost_gamma': 2,
                 'particle_steps_per_undulator_period': 4,
                 'self_field_feedback': True},
        'boundary': {'type': 'cpml', 'cells': [4, 4, 4],
                     'polynomial_order': 3, 'target_reflection': 1e-8,
                     'kappa_max': 8, 'alpha_fraction': 0},
        'initial_self_field': {'enabled': False},
        'beam': {'reference': {'initial_center_z': 1}, 'input': {
            'type': 'generated-gaussian', 'macroparticles': 32768,
            'electrons': 32768, 'gamma': 4, 'direction': [0, 0, 1],
            'center': [0, 0, 0], 'sigma_position': [10, 10, 20],
            'sigma_proper_velocity': [0, 0, 0], 'random_seed': 261026}},
        'sources': {'incident_waves': [], 'magnetic_elements': []},
        'detectors': {'enabled': True, 'directory': 'detectors',
                      'field_planes': [{'name': 'field', 'z': 10,
                                        'rhythm': step_fs * 2.3,
                                        'buffer_samples': 2, 'compression': 0}],
                      'particle_planes': []},
        'stop': {'mode': 'reference-center-z', 'z': 1000},
        'trajectory': {'enabled': False},
        'energy_ledger': {'enabled': True, 'directory': 'ledger',
                          'filename': 'energy-ledger.h5', 'sample_interval_steps': 3,
                          'buffer_records': 4, 'compression': 0}}
    for name in ('continuous', 'restart'):
        folder = root / name
        folder.mkdir()
        current = json.loads(json.dumps(card))
        current['checkpoint'] = {'directory': str(folder / 'checkpoints'),
                                 'steps': [3, 9] if name == 'continuous' else [9],
                                 'io_timeout_seconds': 300}
        if name == 'restart':
            current['checkpoint']['restart_from'] = str(
                root / 'continuous/checkpoints/step-000000000003.h5')
        (folder / 'run.yaml').write_text(yaml.safe_dump(current, sort_keys=False))
    save(root / 'plan.json', {'binary': str(args.binary), 'binary_sha256': sha(args.binary),
                             'nodes': 2, 'ranks': 256, 'threads': 1,
                             'cells': [128, 128, 2048], 'particles': 32768,
                             'step_limit': 12, 'run_timeout_seconds': 420,
                             'scope': 'bounded generated I/O case; not full production'})
    save(root / 'frozen.json', {str(p): sha(p) for p in root.rglob('*') if p.is_file()})
    print(root)


def launch(root, plan, name):
    folder = root / name
    command = ['srun', '--mpi=pmix_v3', '--exact', '--nodes=2', '--ntasks=256',
               '--ntasks-per-node=128', '--cpus-per-task=1', '--cpu-bind=cores',
               '--kill-on-bad-exit=1', '--wait=30', plan['binary'], str(folder / 'run.yaml')]
    start = time.monotonic()
    with (folder / 'aprl.log').open('w') as log:
        result = subprocess.run(['timeout', '--signal=TERM', '--kill-after=30s',
                                 '420s'] + command, cwd=folder, stdout=log, stderr=log)
    log = (folder / 'aprl.log').read_text()
    save(folder / 'launch.json', {'command': command, 'exit_code': result.returncode,
                                 'wall_seconds': time.monotonic() - start})
    assert result.returncode == 1 and 'mesh.duration was exhausted' in log, name
    assert 'phase=duration-exhaustion-readiness-sync-complete' in log, name
    assert 'phase=final step=12 ' in log, name + ': wrong physical endpoint'
    return log


def run(args):
    root = args.root.resolve()
    assert os.environ.get('SLURM_JOB_ID'), 'Slurm allocation required'
    assert int(os.environ['SLURM_JOB_NUM_NODES']) == 2
    plan = json.loads((root / 'plan.json').read_text())
    for path, expected in json.loads((root / 'frozen.json').read_text()).items():
        assert sha(Path(path)) == expected, path
    assert sha(Path(plan['binary'])) == plan['binary_sha256']
    logs = {name: launch(root, plan, name) for name in ('continuous', 'restart')}
    compared = []
    paths = [root / name / 'checkpoints/step-000000000009.h5'
             for name in ('continuous', 'restart')]
    with h5py.File(paths[0]) as a, h5py.File(paths[1]) as b:
        assert a['header'][2] == b['header'][2] == 1
        for family in ('fields', 'cpml', 'particles', 'carriers', 'state/detectors'):
            def check(name, obj, family=family):
                if not isinstance(obj, h5py.Dataset):
                    return
                full = 'blocks/' + family + ('/' + name if name else '')
                other = b[full]
                assert obj.shape == other.shape and obj.dtype == other.dtype, full
                for start in range(0, len(obj), 131072):
                    equal(obj[start:start+131072], other[start:start+131072], full)
                compared.append(full)
            a['blocks/' + family].visititems(check)
    with h5py.File(root / 'continuous/ledger/energy-ledger.h5') as a, \
         h5py.File(root / 'restart/ledger/energy-ledger.h5') as b:
        assert a['energy_ledger/complete'][()] == b['energy_ledger/complete'][()] == 0
        x, y = a['energy_ledger/records'][:], b['energy_ledger/records'][:]
        assert int(a['energy_ledger/committed_records'][()]) == len(x)
        assert int(b['energy_ledger/committed_records'][()]) == len(y)
        equal(x[x['step'] > 3], y, 'ledger post-step3')
    writes = []
    for name, log in logs.items():
        for line in log.splitlines():
            if '[checkpoint] phase=write' in line:
                values = dict(re.findall(r'(step|bytes|seconds|gib_per_s|path)=([^ ]+)', line))
                writes.append({'run': name, **values})
    assert len(writes) == 3 and all(int(w['bytes']) > (1 << 30) for w in writes)
    result = {'status': 'PASS', 'job': os.environ['SLURM_JOB_ID'],
              'nodes': os.environ['SLURM_JOB_NODELIST'], 'ranks': 256, 'threads': 1,
              'state_datasets': compared, 'writes': writes,
              'limitations': 'One moderate 2-node generated case, warm filesystem/cache possible; '
                             'not original production snapshot bandwidth or durability under power loss'}
    save(root / 'RESULT.json', result)
    (root / 'PASS').write_text('Two-node restart state and ledger suffix match\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    for name in ('prepare', 'run'):
        s = sub.add_parser(name)
        s.add_argument('--root', type=Path, required=True)
        if name == 'prepare':
            s.add_argument('--binary', type=Path, required=True)
    args = p.parse_args()
    try:
        (prepare if args.action == 'prepare' else run)(args)
    except Exception as error:
        if args.action == 'run':
            save(args.root / 'FAIL.json', {'error': repr(error)})
        raise
