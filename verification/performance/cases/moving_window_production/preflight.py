#!/usr/bin/env python3
"""Run bounded 16x8 moving-window smoke/replay inside an existing Slurm allocation.

Never submits jobs. Full-input geometry validation is a separate preflight step.
The numerical gates below are fixed before execution; failure preserves evidence.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import time

import h5py
import numpy as np
import yaml

RANKS, THREADS, COUNT = 16, 8, 31675
CASE_LIMIT, TOTAL_LIMIT = 240, 900
BINARY_SHA = 'c81baefb2953ad04181621e31e0f60644872a562bd5bccbd1f45961e44844f11'
FIELD_TOL = 2e-12
ROUND_FACTOR = 128
IDENTITY = ('source_id', 'particle_id', 'charge_C', 'mass_kg', 'weight', 'crossing_direction')
LD = np.longdouble


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def finite(array):
    a = np.asarray(array)
    if a.dtype.names:
        return all(finite(a[name]) for name in a.dtype.names)
    return bool(np.isfinite(a).all()) if a.dtype.kind in 'fc' else True


def launch(command, directory, env, limit, log_name):
    require(limit > 0, 'Time budget exhausted before launch')
    full = ['timeout', '--signal=TERM', '--kill-after=5', str(limit)] + command
    record = {'command': full, 'cwd': str(directory), 'timeout_s': limit,
              'started_unix': time.time(), 'environment': {k: env.get(k) for k in (
                  'OMP_NUM_THREADS', 'OMP_PLACES', 'OMP_PROC_BIND', 'OMP_DYNAMIC',
                  'LD_LIBRARY_PATH', 'SRUN_MPI', 'SLURM_JOB_ID')}}
    start, process = time.monotonic(), None
    try:
        with (directory / log_name).open('x') as stream:
            process = subprocess.Popen(full, cwd=directory, env=env, stdin=subprocess.DEVNULL,
                                       stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            record['exit_code'] = process.wait(timeout=limit + 12)
        require(record['exit_code'] == 0, 'Launch failed: ' + str(directory / log_name))
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        record['seconds'] = time.monotonic() - start
        save(directory / (log_name + '.launch.json'), record)


def validate_binding(path):
    text = path.read_text()
    rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)', text)
    require('[binding] result=PASS failures=0' in text, 'Binding probe did not PASS')
    require(len(rows) == RANKS * THREADS, 'Binding probe must report exactly 128 threads')
    require(all(row[-1] == 'PASS' for row in rows), 'A binding thread failed')
    require({(int(r[0]), int(r[1])) for r in rows} ==
            {(r, t) for r in range(RANKS) for t in range(THREADS)}, 'Binding rank/thread coverage differs')
    require(len({r[2] for r in rows}) == 1, 'Preflight must use exactly one node')
    require(len({(r[2], r[4], r[5]) for r in rows}) == 128, 'Threads share physical cores')
    return {'pass': True, 'host': rows[0][2], 'ranks': RANKS, 'threads_per_rank': THREADS,
            'physical_cores': 128, 'log_sha256': digest(path)}


def inspect_case(directory, free=False):
    log = (directory / 'aprl.log').read_text()
    require('Configured stop reached:' in log, 'Missing physical configured stop')
    require('phase=mpi-finalize-complete' in log, 'MPI finalization did not complete')
    moving = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)', log)
    require(moving is not None, 'Missing moving-window completion')
    require(int(moving[2]) == COUNT, 'Not all particles exercised the moving right exit')
    activation = re.search(r'Frozen lab-plane injection: activations=(\d+),.*still frozen=(\d+)', log)
    require(activation is not None and int(activation[1]) == COUNT and int(activation[2]) == 0,
            'Source release coverage differs')
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    require(not card.get('checkpoint', {}).get('enabled', False), 'Simulation checkpoint enabled')
    # reference/checkpoints.bin is the blocked tridiagonal solver summary, not a simulation snapshot.
    require(not any(p.is_dir() and 'checkpoint' in p.name.lower() for p in directory.iterdir()),
            'Unexpected simulation-checkpoint directory')
    files = sorted(directory.glob('detectors/*.h5')) + sorted(directory.glob('ledger/*.h5')) + sorted(directory.glob('trajectory/*.h5'))
    require(len(files) == 5 + RANKS, 'Unexpected HDF5 file count')
    checks, planes, terminal_counts = {}, {}, np.zeros(COUNT, dtype=np.int64)
    trajectory_ranks = []
    prescribed = None
    for path in files:
        with h5py.File(path, 'r') as f:
            require(len(f) == 1, 'Unexpected HDF5 groups: ' + str(path))
            group = f[next(iter(f))]
            require(int(group['complete'][()]) == 1, 'Incomplete HDF5: ' + str(path))
            for commit, values in (('committed_records', 'records'), ('committed_samples', 'time_s')):
                if commit in group:
                    require(int(group[commit][()]) == len(group[values]), 'Truncated HDF5: ' + str(path))
            datasets = {}
            for name, dataset in group.items():
                if not isinstance(dataset, h5py.Dataset):
                    continue
                datasets[name] = {'shape': list(dataset.shape), 'dtype': str(dataset.dtype)}
                chunks = [dataset[()]] if not dataset.shape else (
                    dataset[i:i + (4 if name == 'fields' else 65536)]
                    for i in range(0, len(dataset), 4 if name == 'fields' else 65536))
                for chunk in chunks:
                    require(finite(chunk), 'Nonfinite HDF5: ' + str(path) + '/' + name)
                    if path.parent.name == 'trajectory' and name == 'records':
                        require(np.all((chunk['source_id'] >= 1) & (chunk['source_id'] <= COUNT)), 'Unexpected trajectory ID')
                        terminal = chunk[chunk['event_type'] != 0]
                        require(np.all(terminal['event_type'] == 2) and np.all(terminal['boundary_face'] == 5),
                                'Unexpected trajectory terminal event')
                        np.add.at(terminal_counts, terminal['source_id'].astype(np.int64) - 1, 1)
            if 'fields' in group:
                ts = group['time_s'][:]
                require(len(ts) > 0 and len(group['fields']) == len(ts) and np.all(np.diff(ts) > 0),
                        'Invalid field sample/time alignment')
            if path.parent.name == 'trajectory':
                require(int(group.attrs['mpi_size']) == RANKS, 'Wrong trajectory MPI size')
                trajectory_ranks.append(int(group.attrs['mpi_rank']))
            if path.parent.name == 'detectors' and path.stem != 'exit-field':
                a = group['records'][:]
                require(len(a) == COUNT, 'Particle plane count differs: ' + path.stem)
                a = a[np.argsort(a['source_id'])]
                require(np.array_equal(a['source_id'], np.arange(1, COUNT + 1)), 'Missing/duplicate plane IDs')
                require(np.all(a['crossing_direction'] == 1), 'Unexpected backward detector crossing')
                planes[path.stem] = a
            if path.parent.name == 'ledger':
                a = group['records'][:]
                require(len(a) > 0 and int(a['step'][-1]) == int(moving[1]), 'Ledger does not reach final step')
                prescribed = float(np.max(np.abs(a['prescribed_work_J'])))
                if free:
                    require(prescribed == 0., 'No-magnet prescribed work is not exactly zero')
            checks[str(path.relative_to(directory))] = {'complete': True, 'datasets': datasets,
                                                        'size_bytes': path.stat().st_size}
    require(set(planes) == {'release', 'aux', 'exit'}, 'Missing particle plane')
    for a in planes.values():
        for key in IDENTITY:
            require(np.array_equal(a[key], planes['release'][key]), 'Particle cohort identity changed: ' + key)
    require(sorted(trajectory_ranks) == list(range(RANKS)), 'Missing trajectory rank')
    require(np.all(terminal_counts == 1), 'Missing/duplicate final particle exit events')
    require(prescribed is not None, 'Missing energy ledger')
    return {'pass': True, 'steps': int(moving[1]), 'particle_exits': int(moving[2]),
            'cropped_carriers': int(moving[3]), 'plane_counts': {k: len(a) for k, a in planes.items()},
            'prescribed_work_max_abs_J': prescribed, 'terminal_events': int(terminal_counts.sum()),
            'simulation_checkpoint_enabled': False, 'hdf5_files': checks}, planes


def particle_error(first, second):
    require(first.shape == second.shape, 'Particle comparison shapes differ')
    a, b = np.asarray(first, dtype=LD), np.asarray(second, dtype=LD)
    delta = b - a
    scale = np.max(np.abs(a), axis=0)
    tolerance = ROUND_FACTOR * LD(np.finfo(np.float64).eps) * scale
    excess = abs(delta) > tolerance
    norm = np.sum(a * a)
    return {'pass': not bool(np.any(excess)), 'bitwise_equal': first.dtype == second.dtype and first.tobytes() == second.tobytes(),
            'maximum_abs_difference': float(np.max(abs(delta))),
            'rms_abs_difference': float(np.sqrt(np.mean(delta * delta))),
            'relative_l2': float(np.sqrt(np.sum(delta * delta) / norm)) if norm else None,
            'reference_component_max_abs': np.asarray(scale, dtype=float).tolist(),
            'per_component_absolute_roundoff_gate': np.asarray(tolerance, dtype=float).tolist(),
            'elements_over_roundoff_gate': int(np.count_nonzero(excess))}


def physical_card(path):
    card = yaml.safe_load(Path(path).read_text())
    card['beam']['input']['file'] = '<INPUT>'
    card['mesh']['moving_window']['reference_directory'] = '<REFERENCE>'
    card['output']['manifest'] = '<MANIFEST>'
    for section in ('detectors', 'trajectory', 'energy_ledger'):
        card[section]['directory'] = '<OUTPUT>'
    return card


def compare(reference, candidate, planes):
    require(physical_card(reference / 'run.yaml') == physical_card(candidate / 'run.yaml'),
            'Replay physical cards differ beyond output/input/reference paths')
    result = {'reference': str(reference), 'field_relative_l2_limit': FIELD_TOL,
              'particle_roundoff_factor_eps64': ROUND_FACTOR, 'planes': {}}
    for name, b in planes.items():
        path = reference / 'detectors' / (name + '.h5')
        with h5py.File(path) as f:
            g = f[next(iter(f))]
            require(int(g['complete'][()]) == 1 and int(g['committed_records'][()]) == COUNT,
                    'Reference particle plane incomplete')
            a = g['records'][:]
            a = a[np.argsort(a['source_id'])]
        require(finite(a), 'Reference plane nonfinite')
        for key in IDENTITY:
            require(np.array_equal(a[key], b[key]), 'Cross-layout identity changed: ' + key)
        metrics = {key: particle_error(a[key], b[key]) for key in ('position_m', 'proper_velocity', 'time_s')}
        ua, ub = np.asarray(a['proper_velocity'], dtype=LD), np.asarray(b['proper_velocity'], dtype=LD)
        dk = np.sum((ub - ua) * (ub + ua), axis=1) / (
            np.sqrt(1 + np.sum(ua * ua, axis=1)) + np.sqrt(1 + np.sum(ub * ub, axis=1))) * LD('510998.9499961642')
        metrics['delta_K'] = {'mean_eV': float(np.mean(dk)), 'rms_eV': float(np.sqrt(np.mean(dk * dk))),
                              'max_abs_eV': float(np.max(abs(dk))), 'scope': 'Reported, no independent energy gate added'}
        result['planes'][name] = metrics
    with h5py.File(reference / 'detectors/exit-field.h5') as fa, h5py.File(candidate / 'detectors/exit-field.h5') as fb:
        a, b = fa[next(iter(fa))], fb[next(iter(fb))]
        require(int(a['complete'][()]) == 1 and int(a['committed_samples'][()]) == len(a['time_s']),
                'Reference field incomplete')
        require(np.array_equal(a['time_s'][:], b['time_s'][:]), 'Cross-layout field times changed')
        require(a['fields'].shape == b['fields'].shape, 'Cross-layout field shapes changed')
        for key in ('x_first_m', 'y_first_m', 'dx_m', 'dy_m'):
            require(a.attrs[key] == b.attrs[key], 'Field plane geometry changed: ' + key)
        result['fields'] = {}
        for name in ('electric_V_per_m', 'magnetic_T'):
            norm, error, maxabs, equal = LD(0), LD(0), LD(0), True
            for i in range(0, len(a['fields']), 4):
                aa, bb = a['fields'][i:i+4][name], b['fields'][i:i+4][name]
                require(finite(aa) and finite(bb), 'Nonfinite compared field')
                equal = equal and aa.tobytes() == bb.tobytes()
                x, d = np.asarray(aa, dtype=LD), np.asarray(bb, dtype=LD) - np.asarray(aa, dtype=LD)
                norm += np.sum(x * x)
                error += np.sum(d * d)
                maxabs = max(maxabs, np.max(abs(d)))
            rel = float(np.sqrt(error / norm)) if norm else (0. if error == 0 else None)
            result['fields'][name] = {'relative_l2': rel, 'max_abs_difference': float(maxabs),
                                     'bitwise_equal': equal, 'pass': rel is not None and rel <= FIELD_TOL}
    result['pass'] = all(v['pass'] for v in result['fields'].values()) and all(
        v[key]['pass'] for v in result['planes'].values() for key in ('position_m', 'proper_velocity', 'time_s'))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    args = p.parse_args()
    campaign = args.root.resolve()
    plan = json.loads((campaign / 'plan.json').read_text())
    root = campaign / 'preflight'
    root.mkdir(exist_ok=True)
    require(not (root / 'SMOKE_STARTED.json').exists(), 'Smoke already started; preserve evidence and use a new campaign')
    require(os.environ.get('SLURM_JOB_ID'), 'Must run inside a real Slurm allocation')
    require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == 1, 'Require one allocated node')
    require(int(os.environ.get('SLURM_NTASKS', '0')) == RANKS and
            int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == THREADS, 'Require allocated16MPIx8CPU')
    require(plan['ranks'] == RANKS and plan['threads'] == THREADS, 'Campaign layout differs')
    binary, binding = Path(plan['binary']), Path(plan['binding_probe'])
    require(plan['binary_sha256'] == BINARY_SHA and digest(binary) == BINARY_SHA, 'R3 binary hash mismatch')
    require(os.access(binary, os.X_OK) and os.access(binding, os.X_OK), 'Missing executable')
    work = Path(plan['workspace'])
    old = work / 'results/moving_window_integration_20261028'
    source_card = old / 'dense-two-node-r3/moving/run.yaml'
    source_input = old / 'dense-two-node-r3/inputs/held.h5'
    reference = old / 'dense-r1/moving'
    save(root / 'THRESHOLDS.json', {'field_relative_l2': FIELD_TOL,
        'field_scope': 'Existing MPI module roundoff bound, not physical radiation accuracy',
        'particle_formula': 'Every element: abs(candidate-reference) <=128*eps64*max_abs(reference component over this plane)',
        'particle_roundoff_factor': ROUND_FACTOR, 'particle_absolute_floor': 0,
        'identities': 'Exact source_id/particle_id/charge/mass/weight/crossing_direction; 31675 unique IDs on all three planes',
        'no_magnet_prescribed_work_J': 0., 'case_timeout_s': CASE_LIMIT, 'total_timeout_s': TOTAL_LIMIT})
    save(root / 'SMOKE_STARTED.json', {'started_unix': time.time(), 'host': socket.gethostname(),
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'slurm_nodes': os.environ.get('SLURM_JOB_NODELIST'),
        'allocated_cpus': 128, 'ranks': RANKS, 'threads': THREADS})
    start = time.monotonic()
    deadline = start + TOTAL_LIMIT
    def deadline_expired(signum, frame):
        raise TimeoutError('Total900s smoke budget exceeded')
    signal.signal(signal.SIGALRM, deadline_expired)
    signal.setitimer(signal.ITIMER_REAL, TOTAL_LIMIT)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(THREADS), OMP_DYNAMIC='FALSE', OMP_PLACES='cores', OMP_PROC_BIND='close',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', APRL_PROFILE_INTERVAL='0')
    libs = binary.parent / 'runtime-libs'
    require(libs.is_dir(), 'R3 runtime-libs directory missing')
    env['LD_LIBRARY_PATH'] = str(libs) + ':' + env.get('LD_LIBRARY_PATH', '')
    launcher = ['srun', '--mpi=' + env.get('SRUN_MPI', 'pmix_v3'), '--exact', '--nodes=1', '--ntasks=16',
                '--ntasks-per-node=16', '--cpus-per-task=8', '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']
    try:
        inputs = root / 'inputs'
        inputs.mkdir()
        shutil.copy2(source_input, inputs / 'held.h5')
        shutil.copy2(source_card, inputs / 'source-card.yaml')
        with h5py.File(inputs / 'held.h5') as f:
            require(len(f['particles/records']) == COUNT, 'Wrong smoke input count')
        frozen = {str(path): digest(path) for path in (binary, binding, source_card, source_input,
                  inputs / 'held.h5', inputs / 'source-card.yaml', reference / 'run.yaml', Path(__file__).resolve())}
        frozen.update({str(path): digest(path) for path in libs.iterdir() if path.is_file()})
        require(digest(inputs / 'held.h5') == digest(old / 'dense-r1/inputs/held.h5'), 'Reference fixture input differs')
        save(root / 'FROZEN_INPUTS.json', frozen)
        launch(launcher + ['--chdir=' + str(root), str(binding), str(RANKS), str(THREADS)], root, env, 45, 'binding.log')
        binding_result = validate_binding(root / 'binding.log')
        summaries, comparison = {}, None
        for name, free in (('small-undulator', False), ('small-free-drift', True)):
            directory = root / name
            directory.mkdir()
            card = yaml.safe_load(source_card.read_text())
            card['beam']['input']['file'] = str(inputs / 'held.h5')
            card['mesh']['moving_window']['reference_directory'] = str(directory / 'reference')
            card['output']['manifest'] = str(directory / 'run-manifest.yaml')
            for section in ('detectors', 'trajectory', 'energy_ledger'):
                card[section]['directory'] = str(directory / ('ledger' if section == 'energy_ledger' else section))
            require('checkpoint' not in card, 'Frozen short fixture unexpectedly contains checkpoint')
            if free:
                card['sources']['magnetic_elements'] = []
            (directory / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
            frozen[str(directory / 'run.yaml')] = digest(directory / 'run.yaml')
            save(root / 'FROZEN_INPUTS.json', frozen)
            remaining = int(deadline - time.monotonic() - 20)
            require(remaining >= 60, 'Insufficient smoke time budget')
            launch(launcher + ['--chdir=' + str(directory), str(binary), str(directory / 'run.yaml')],
                   directory, env, min(CASE_LIMIT, remaining), 'aprl.log')
            summaries[name], planes = inspect_case(directory, free=free)
            save(directory / 'COMPLETE_REVIEW.json', summaries[name])
            if not free:
                comparison = compare(reference, directory, planes)
                save(root / 'UNDULATOR_REPLAY.json', comparison)
                require(comparison['pass'], 'Cross-layout numerical replay gate failed; see UNDULATOR_REPLAY.json')
            require(time.monotonic() < deadline, 'Total900s smoke budget exceeded')
        require(all(digest(path) == value for path, value in frozen.items()), 'Frozen input changed during smoke')
        save(root / 'SMOKE_PASS.json', {'pass': True, 'status': 'PASS', 'binding': binding_result,
            'binary': str(binary), 'binary_sha256': BINARY_SHA, 'cases': summaries,
            'undulator_replay': comparison, 'elapsed_s': time.monotonic() - start,
            'limits': ['Short16x8 functional and MPI roundoff replay only; not full-input throughput or weak-radiation acceptance.',
                       'Full-input geometry/MPI-count checks are a separate gate.']})
        print('SMOKE_PASS ' + str(root / 'SMOKE_PASS.json'), flush=True)
    except BaseException as error:
        save(root / 'SMOKE_FAILED.json', {'pass': False, 'error': str(error), 'elapsed_s': time.monotonic() - start})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
