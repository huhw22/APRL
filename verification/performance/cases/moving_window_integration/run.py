#!/usr/bin/env python3
"""Prepare/run bounded full-program moving-window fixtures; never submits jobs."""
import argparse
import copy
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

import h5py
import numpy as np
import yaml

WORK = Path('/public/home/wfollowm/APRL_Workspace')
SOURCE = WORK / 'source/aprl'
HERE = Path(__file__).resolve().parent
C, E = 299792458., 1.602176634e-19
GAMMA, DZ = 929.1305492858543, 1.2543262415359034e-6
BETA = math.sqrt(1 - GAMMA**-2)
DT = DZ / C
PERIOD, RELEASE, FIELD = .0234, -.05, .06


def require(test, message):
    if not test:
        raise RuntimeError(message)


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fixture(path, count):
    # Sample moments from the frozen 2026-10-28 source audit, not a full input read.
    rng = np.random.default_rng(281028)
    q = np.clip(rng.normal(size=(count, 6)), -4., 4.)
    q[0, 2], q[1, 2] = -4.5, 4.5  # Explicit inlet witnesses; not random-tail statistics.
    xy_sigma = np.array([6.393651865400727e-5, 6.126395728478917e-5])
    u_mean = np.array([3.526198929811919e-6, 2.294658317104036e-5, 1174.1852350615736])
    u_sigma = np.array([.02709500523622301, .02543628834980338, 2.568639437850893])
    corr = np.array([-.7334153049192772, -.6546423205475509])
    dtype = [('plane_position_m', '<f8', (2,)), ('arrival_time_offset_s', '<f8'),
             ('proper_velocity', '<f8', (3,)), ('source_id', '<u8'), ('macro_weight', '<f8')]
    data = np.zeros(count, dtype=dtype)
    data['plane_position_m'] = q[:, :2] * xy_sigma
    data['proper_velocity'][:, :2] = u_mean[:2] + u_sigma[:2] * (
        corr * q[:, :2] + np.sqrt(1-corr**2) * q[:, 3:5])
    data['proper_velocity'][:, 2] = u_mean[2] + u_sigma[2] * q[:, 5]
    beta_lab = u_mean[2] / math.hypot(1., u_mean[2])
    held_beta = (beta_lab-BETA)/(1-BETA*beta_lab)
    sigma_time = 8 * DZ / (GAMMA*C*(BETA+held_beta))
    data['arrival_time_offset_s'] = q[:, 2] * sigma_time
    data['source_id'], data['macro_weight'] = np.arange(1, count+1), 1.
    with h5py.File(path, 'x') as stream:
        group = stream.create_group('particles')
        group.attrs['format_version'] = 4
        group.attrs['fixture_scope'] = 'fixed-seed functional fixture with two explicit +/-4.5 sigma inlet witnesses'
        group.create_dataset('records', data=data)
    electrons = 41732775. * 8. / 10540.
    # Reconstruct the v4 common-time projection and the held frame's earliest
    # crossing origin. The reference-centre boost worldline does NOT start at
    # initial_center_z on the t'=0 plane; using that would underestimate duration.
    u = np.asarray(data['proper_velocity'], dtype=np.longdouble)
    gamma_lab = np.sqrt(1+np.sum(u*u, axis=1))
    beta_z = u[:, 2]/gamma_lab
    arrival = np.asarray(data['arrival_time_offset_s'], dtype=np.longdouble)
    center_offset = np.longdouble(float(RELEASE+1e-7)-RELEASE)
    reference_time = (center_offset/C+np.mean(beta_z*arrival))/np.mean(beta_z)
    input_time = arrival-reference_time
    # Stable light-front form used by BoostFrameTransform. Do not change the
    # already established arrival-time fixture above when improving estimates.
    gamma_beta = np.sqrt((np.longdouble(GAMMA)-1)*(np.longdouble(GAMMA)+1))
    plus = (np.sqrt(1+u[:, 2]**2)+u[:, 2])/(GAMMA+gamma_beta)
    held_u = .5*(plus-1/plus)
    held = held_u/np.sqrt(1+held_u**2)
    frame_t0 = np.min(input_time-2*DZ/(GAMMA*C*(BETA+held)))
    reference_z0 = RELEASE+np.longdouble(GAMMA)**2*(center_offset+BETA*C*frame_t0)
    source_z_cells = -GAMMA*C*(BETA+held)*(input_time-frame_t0)/DZ
    translation_cells = -.5*(source_z_cells.min()+source_z_cells.max())
    require(abs(float(np.mean(held))) <= .25, 'Fixture held beta exceeds the explicit reference bound')
    require(source_z_cells.min() < -64, 'Fixture does not exercise held particles behind the inlet')
    return {'records': count, 'seed': 281028, 'represented_electrons': electrons,
            'macro_electrons': electrons/count, 'sigma_arrival_time_s': sigma_time,
            'nominal_sigma_boost_z_cells': 8., 'witness_source_ids': [1, 2],
            'witness_arrival_sigma': [-4.5, 4.5],
            'approximate_mean_held_beta': held_beta,
            'reconstructed_mean_held_beta_estimate': float(np.mean(held)),
            'reference_initial_lab_z_estimate_m': float(reference_z0),
            'initial_source_z_cells_estimate': [float(source_z_cells.min()), float(source_z_cells.max())],
            'held_box_translation_cells_estimate': float(translation_cells),
            'held_box_translation_m_estimate': float(translation_cells*DZ),
            'witness_initial_z_cells_estimate': {
                str(identifier): float(source_z_cells[identifier-1]+translation_cells)
                for identifier in (1, 2)},
            'witness_held_beta_estimate': {str(identifier): float(held[identifier-1]) for identifier in (1, 2)},
            'moment_source': 'results/held_cut_20261028/source-moment-sample.json',
            'scope': 'Functional input, not production-density or N-convergence validation'}


def geometry(source):
    full_nz = 896
    origin = -.5 * full_nz * DZ
    translation = source['held_box_translation_m_estimate']
    release_box = translation
    offset = math.floor((release_box-origin)/DZ) - 64
    cells = math.ceil((FIELD-RELEASE)/(GAMMA*DZ)) + 64 + 128 + 2
    # Match the runtime's conservative ceil(duration/dt)+1 preflight bound.
    # Right-margin ceil is one for all of these short fixtures. Recheck after
    # deriving the duration from the actual held reference-centre worldline.
    right_margin = 5
    outer = RELEASE + GAMMA*(origin+(offset+cells)*DZ-translation)
    safe_exit = outer - right_margin*GAMMA*DZ
    stop = safe_exit + 6*GAMMA*DZ
    initial_center = RELEASE + 1e-7
    stop_step = math.ceil((stop-source['reference_initial_lab_z_estimate_m'])/(GAMMA*BETA*DZ))
    max_steps = stop_step+1
    duration = (max_steps-.5)*DT
    bound_steps = math.ceil(duration/DT)+1
    require(right_margin == 4+math.ceil((1-BETA)*bound_steps), 'Right-margin estimate changed')
    # beta=.25 bounds the fixture's small longitudinal momentum spread; exact
    # measured beta and reference range are validated again by the executable.
    minimum = offset-(1+.25)*bound_steps-1
    require(minimum >= 0 and offset+cells <= full_nz-1, 'Reference range is insufficient')
    require(stop_step < max_steps and stop > outer, 'Stop does not exercise outer cropping')
    initial_origin = origin+offset*DZ
    witness_relative = {key: value-initial_origin/DZ for key, value in source['witness_initial_z_cells_estimate'].items()}
    require(witness_relative['2'] < 0, 'Source ID2 must remain outside the initial moving window')
    crossings = {}
    for key, relative in witness_relative.items():
        if relative >= 0:
            continue
        held_beta = source['witness_held_beta_estimate'][key]
        for step in range(1, max_steps+1):
            before = relative+(step-1)*(1+held_beta)
            after = before+held_beta
            if before < 0 <= after or after < 0 <= after+1:
                crossings[key] = {'step': step, 'relative_before_drift_cells': before,
                    'relative_after_drift_cells': after,
                    'entry_type': 'physical_drift' if after >= 0 else 'window_shift'}
                break
    require(crossings.get('2', {}).get('entry_type') == 'physical_drift', 'ID2 must exercise physical left entry')
    return {'geometry_contract_version': 2,
            'reference_cells': full_nz, 'active_cells_estimate': cells,
            'initial_offset_estimate': offset, 'duration_box_s': duration,
            'held_box_translation_m_estimate': translation,
            'initial_origin_box_z_estimate_m': initial_origin,
            'witness_relative_to_initial_inlet_cells': witness_relative,
            'witness_entry_estimates': crossings,
            'max_steps': max_steps, 'reference_preflight_steps': bound_steps,
            'reference_minimum_query_bound_cells': minimum,
            'reference_maximum_query_bound_cells': offset+cells,
            'bound_mean_held_beta': .25, 'outer_lab_z_estimate_m': outer,
            'safe_exit_lab_z_estimate_m': safe_exit, 'stop_reference_lab_z_m': stop,
            'predicted_reference_stop_step': stop_step, 'initial_center_lab_z_m': initial_center,
            'fringe_start_lab_z_m': -math.sqrt(-2*math.log(1e-9))*PERIOD/(2*math.pi),
            'fringe_end_lab_z_m': PERIOD+math.sqrt(-2*math.log(1e-9))*PERIOD/(2*math.pi)}


def card(root, geo, source):
    plane = lambda name, z: {'name': name, 'z': z, 'buffer_records': source['records'], 'compression': 0}
    return {
        'units': {'length': 'm', 'time': 's'},
        'runtime': {'mode': 'throughput', 'resource_monitor': {
            'enabled': True, 'progress_interval_steps': 20, 'calibration_steps': 5},
            'shutdown_watchdog': {'enabled': True, 'timeout_seconds': 600, 'finalize_timeout_seconds': 60}},
        'output': {'overwrite': False, 'manifest': 'run-manifest.yaml'},
        'mesh': {'field_solver': 'cowan-z', 'current_filter': 'none',
            'ballistic_outside_transverse': False, 'self_field_feedback': True,
            'cells': [160, 160, geo['reference_cells']], 'cell_size': [5e-6, 5e-6, DZ],
            'center': [0., 0., 0.], 'duration': geo['duration_box_s'], 'boost_gamma': GAMMA,
            'particle_steps_per_undulator_period': 80, 'maximum_particle_substeps': 4096},
        'boundary': {'type': 'cpml', 'cells': [10, 10, 0], 'polynomial_order': 3.,
            'target_reflection': 1e-8, 'kappa_max': 8., 'alpha_fraction': 0.},
        'initial_self_field': {'enabled': True, 'model': 'cowan-relativistic-poisson',
            'relative_tolerance': 1e-8, 'maximum_iterations': 10000,
            'maximum_mean_transverse_beta': 3e-5, 'include_transverse_cpml_particles': False},
        'beam': {'reference': {'input_plane_z': RELEASE, 'initial_center_z': geo['initial_center_lab_z_m'],
            'frozen_lab_plane_injection': True, 'hold_transverse_until_activation': True,
            'activation_plane_z': RELEASE}, 'input': {'type': 'hdf5',
            'file': str(root/'inputs/held.h5'), 'electrons': source['represented_electrons'],
            'position_offset': [0., 0., 0.]}},
        'sources': {'incident_waves': [], 'magnetic_elements': [{'type': 'planar-undulator',
            'characteristic': True, 'strength_parameter': 1.092714657005, 'period': PERIOD,
            'periods': 1, 'entrance_z': 0., 'polarization_angle_rad': math.pi/2,
            'gaussian_fringe': True, 'fringe_relative_cutoff': 1e-9}]},
        'detectors': {'enabled': True, 'directory': 'detectors',
            'field_planes': [{'name': 'exit-field', 'z': FIELD, 'rhythm': DT/GAMMA,
                'buffer_samples': 8, 'compression': 0, 'particle_background': {'enabled': False}}],
            'particle_planes': [plane('release', RELEASE),
                plane('aux', RELEASE+16*GAMMA*DZ), plane('exit', FIELD)]},
        'stop': {'mode': 'reference-center-z', 'z': geo['stop_reference_lab_z_m']},
        'trajectory': {'enabled': True, 'directory': 'trajectory', 'basename': 'trajectory',
            'rhythm': DT, 'buffer_records': 4096, 'compression': 0},
        'energy_ledger': {'enabled': True, 'directory': 'ledger', 'filename': 'energy-ledger.h5',
            'sample_interval_steps': 1, 'buffer_records': 128, 'compression': 0,
            'warning_relative_tolerance': .05}}


def prepare(args):
    root, binary = args.root.resolve(), args.binary.resolve()
    require(not root.exists(), 'Use a new result directory')
    require(args.nodes <= args.ranks and args.ranks % args.nodes == 0, 'Ranks must divide over nodes')
    binding = args.binding_probe.resolve()
    for path in (binary, binding):
        require(path.is_file() and os.access(path, os.X_OK), 'Missing executable '+str(path))
    root.mkdir(parents=True)
    (root/'inputs').mkdir()
    source = fixture(root/'inputs/held.h5', args.particles)
    geo = geometry(source)
    geometry_replay = None
    if args.geometry_from:
        previous_path = args.geometry_from.resolve()
        previous = json.loads(previous_path.read_text())
        require(previous.get('format') in ('aprl-moving-window-integration-v1', 'aprl-moving-window-integration-v2'),
            'Unknown geometry replay plan format')
        for key in ('records', 'seed', 'represented_electrons', 'macro_electrons', 'sigma_arrival_time_s'):
            require(previous['source'][key] == source[key], 'Geometry replay fixture mismatch: '+key)
        old_input = previous_path.parent/'inputs/held.h5'
        with h5py.File(old_input) as old, h5py.File(root/'inputs/held.h5') as new:
            require(np.array_equal(old['particles/records'][:], new['particles/records'][:]),
                'Geometry replay requires identical complete fixture records')
        previous_card_path = previous_path.parent/previous['variants'][0]/'run.yaml'
        previous_card = yaml.safe_load(previous_card_path.read_text())
        replay_geometry = copy.deepcopy(previous['geometry'])
        proposed_card = card(root, replay_geometry, source)
        def normalized(value):
            value = copy.deepcopy(value)
            value['mesh'].pop('moving_window', None)
            value['beam']['input']['file'] = '<verified identical input>'
            value['output']['manifest'] = '<output path>'
            for name in ('detectors', 'trajectory', 'energy_ledger'):
                value[name]['directory'] = '<output path>'
            return value
        require(normalized(previous_card) == normalized(proposed_card),
            'Geometry replay changes units, physical constants, or another non-path card setting')
        geometry_replay = {'plan': str(previous_path), 'plan_sha256': sha(previous_path),
            'card': str(previous_card_path), 'card_sha256': sha(previous_card_path),
            'all_fixture_records_equal': True, 'nonpath_card_settings_equal': True,
            'purpose': 'Exact frozen duration/stop/source-grid replay for layout comparison; original geometry estimates remain labelled as estimates.'}
        shutil.copy2(previous_path, root/'replayed-plan.json')
        shutil.copy2(previous_card_path, root/'replayed-source-card.yaml')
        geo = replay_geometry
    base = card(root, geo, source)
    variants = ['moving', 'fixed']
    binaries = {name: str(binary) for name in variants}
    if args.baseline:
        old = args.baseline.resolve()
        require(old.is_file() and os.access(old, os.X_OK), 'Missing old executable')
        variants.append('old-fixed')
        binaries['old-fixed'] = str(old)
    if args.variants:
        selected = args.variants.split(',')
        require(selected and len(selected) == len(set(selected)) and set(selected) <= set(variants),
            'Unknown/duplicate variant or old-fixed without --baseline')
        variants = selected
        binaries = {name: binaries[name] for name in variants}
    for name in variants:
        directory = root/name
        directory.mkdir()
        current = copy.deepcopy(base)
        if name == 'moving':
            current['mesh']['moving_window'] = {'enabled': True, 'lower_guard_cells': 64,
                'upper_guard_cells': 128, 'block_layers': 32,
                'reference_directory': str(directory/'reference')}
        (directory/'run.yaml').write_text(yaml.safe_dump(current, sort_keys=False))
    for name in ('run.py', 'analyze.py', 'README.txt', 'run.sbatch', 'run-two-nodes.sbatch'):
        shutil.copy2(HERE/name, root/name)
    (root/'source.diff').write_bytes(subprocess.check_output(['git', '-C', str(SOURCE), 'diff', '--binary', 'HEAD']))
    source_files = [p for p in sorted((SOURCE/'src').glob('*')) if p.is_file()]
    source_files += [SOURCE/'CMakeLists.txt', SOURCE/'verification/CMakeLists.txt']
    save(root/'source-files.json', {str(p.relative_to(SOURCE)): sha(p) for p in source_files})
    for path in source_files:
        destination = root/'source'/path.relative_to(SOURCE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    plan = {'format': 'aprl-moving-window-integration-v2', 'variants': variants, 'binaries': binaries,
        'binding': str(binding), 'nodes': args.nodes, 'ranks': args.ranks, 'threads': args.threads,
        'case_timeout_s': 1200, 'total_timeout_s': 1500,
        'source_head': subprocess.check_output(['git', '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip(),
        'source': source, 'geometry': geo,
        'geometry_replay': geometry_replay,
        'scope': 'Functional integration, fixed-to-moving approximation reported without an invented physical PASS'}
    save(root/'plan.json', plan)
    frozen = [p for p in root.rglob('*') if p.is_file()] + [binding] + [Path(p) for p in binaries.values()]
    for filename in binaries.values():
        frozen += [p for p in (Path(filename).parent/'runtime-libs').glob('*') if p.is_file()]
    save(root/'frozen.json', {str(p): sha(p) for p in set(frozen)})
    print(root/'plan.json')


def launch(command, directory, env, limit, log_name):
    if command and Path(command[0]).name == 'srun':
        command = [command[0], '--chdir='+str(directory.resolve())]+command[1:]
    command = ['timeout', '--signal=TERM', '--kill-after=5', str(limit)] + command
    record = {'command': command, 'timeout_s': limit, 'started_unix': time.time()}
    started = time.monotonic()
    process = None
    try:
        with (directory/log_name).open('x') as stream:
            process = subprocess.Popen(command, cwd=directory, env=env, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            record['exit_code'] = process.wait(timeout=limit+10)
        require(record['exit_code'] == 0, 'Failed launch: '+str(directory/log_name))
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        record['seconds'] = time.monotonic()-started
        save(directory/(log_name+'.launch.json'), record)


def run(args):
    root = args.root.resolve()
    require(not (root/'STARTED.json').exists(), 'Campaign already started')
    for filename, digest in json.loads((root/'frozen.json').read_text()).items():
        require(sha(filename) == digest, 'Frozen input changed: '+filename)
    plan = json.loads((root/'plan.json').read_text())
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES', 0)) == plan['nodes'],
        'This full-program campaign requires its recorded Slurm allocation')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) == 128, 'Reserve exclusive 128 CPU nodes')
    save(root/'STARTED.json', {'host': socket.gethostname(), 'unix_time': time.time(),
        'slurm_job_id': os.environ['SLURM_JOB_ID'], 'allocated_cpus': 128*plan['nodes'],
        'used_cpus': plan['ranks']*plan['threads'], 'nodes': os.environ.get('SLURM_JOB_NODELIST')})
    launcher = ['srun', '--mpi='+os.environ.get('SRUN_MPI', 'pmix_v3'), '--exact',
        '--nodes='+str(plan['nodes']), '--ntasks='+str(plan['ranks']),
        '--ntasks-per-node='+str(plan['ranks']//plan['nodes']), '--cpus-per-task='+str(plan['threads']),
        '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(plan['threads']), OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
        OMP_PROC_BIND='close', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
        APRL_PROFILE_INTERVAL='0')
    deadline = time.monotonic()+plan['total_timeout_s']
    try:
        launch(launcher+[plan['binding'], str(plan['ranks']), str(plan['threads'])], root, env, 45, 'binding.log')
        require('[binding] result=PASS failures=0' in (root/'binding.log').read_text(), 'CPU binding failed')
        for name in plan['variants']:
            remaining = int(deadline-time.monotonic()-10)
            require(remaining >= 60, 'Campaign time budget exhausted')
            binary = Path(plan['binaries'][name])
            run_env = env.copy()
            if (binary.parent/'runtime-libs').is_dir():
                run_env['LD_LIBRARY_PATH'] = str(binary.parent/'runtime-libs')+':'+run_env.get('LD_LIBRARY_PATH', '')
            launch(launcher+[str(binary), str(root/name/'run.yaml')], root/name, run_env,
                min(plan['case_timeout_s'], remaining), 'aprl.log')
            require('Configured stop reached:' in (root/name/'aprl.log').read_text(), 'Missing configured stop')
            print('FINISHED '+name, flush=True)
        save(root/'COMPLETED.json', {'status': 'EXECUTION_COMPLETE_SCIENTIFIC_REVIEW_REQUIRED', 'unix_time': time.time()})
    except BaseException as error:
        save(root/'FAILED.json', {'error': str(error), 'unix_time': time.time()})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--binary', type=Path, required=True)
    prep.add_argument('--baseline', type=Path, help='Optional earlier binary for disabled-path replay; all variants share the 1500s budget')
    prep.add_argument('--binding-probe', type=Path, default=WORK/'build/cluster-tools/binding_probe')
    prep.add_argument('--nodes', type=int, choices=(1, 2), default=1)
    prep.add_argument('--ranks', type=int, choices=(1, 2), default=1)
    prep.add_argument('--threads', type=int, choices=range(1, 17), default=16)
    prep.add_argument('--particles', type=int, choices=(64, 31675), default=64)
    prep.add_argument('--variants', help='Explicit frozen comma-separated subset/order: moving,fixed,old-fixed')
    prep.add_argument('--geometry-from', type=Path,
        help='Replay an existing compatible frozen plan geometry exactly after verifying complete source records and non-path card settings')
    execute = sub.add_parser('run')
    execute.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    (prepare if args.action == 'prepare' else run)(args)
