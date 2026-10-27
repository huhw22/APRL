#!/usr/bin/env python3
"""Small full-program Cowan initializer integration checks; preparation is read-only computation."""
import argparse
import copy
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
HELPER = HERE / 'checkpoint_helpers.py'
if not HELPER.exists():
    HELPER = HERE.parents[2] / 'checkpoint_restart/run.py'
SPEC = importlib.util.spec_from_file_location('cowan_checkpoint_helpers', HELPER)
common = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(common)
require, save, digest = common.require, common.save, common.digest
OLD_MODEL = 'relativistic-poisson'
NEW_MODEL = 'cowan-relativistic-poisson'
VARIANTS = ('accepted', 'legacy', 'cowan')
E_CHARGE = 1.602176634e-19
REST_EV = 510998.9499961642


def cards(root, include_held):
    base = common.gaussian_card()
    base['mesh'].update(cells=[32, 32, 192], cell_size=[2, 2, .5], duration=600)
    base['boundary']['cells'] = [4, 4, 8]
    base['beam']['reference']['initial_center_z'] = 0
    base['beam']['input'].update(
        electrons=1e-12 / E_CHARGE, macroparticles=4096,
        gamma=2, sigma_position=[2, 2, 1], random_seed=271027)
    base['initial_self_field'].update(relative_tolerance=1e-10, maximum_iterations=2000)
    base['energy_ledger']['sample_interval_steps'] = 1
    base['trajectory'].update(rhythm=10, buffer_records=8192)
    base['detectors']['particle_planes'] = [
        {'name': 'downstream', 'z': 60, 'buffer_records': 4096, 'compression': 0}]
    base['detectors']['field_planes'][0].update(rhythm=.01, buffer_samples=4)
    static = copy.deepcopy(base)
    static['boundary'] = {'type': 'pec'}
    moving = copy.deepcopy(base)
    beta = .23
    moving['beam']['input']['gamma'] = 2 * (1 + math.sqrt(.75) * beta) / math.sqrt(1 - beta**2)
    result = {'gaussian-static-pec': static, 'gaussian-moving-cpml': moving}
    if include_held:
        # Runtime requires HDF5-v4 arrival events for actual lab-plane held.
        # This small established fixture supplements the generated-Gaussian base.
        held = common.gaussian_card()
        held['beam'] = {'reference': {
            'initial_center_z': 1, 'input_plane_z': 0,
            'frozen_lab_plane_injection': True, 'activation_plane_z': 0,
            'hold_transverse_until_activation': True}, 'input': {
            'type': 'hdf5', 'file': str(root / 'inputs/held.h5'),
            'electrons': 64, 'position_offset': [0, 0, 0]}}
        held['energy_ledger']['sample_interval_steps'] = 1
        held['detectors']['field_planes'][0]['rhythm'] = .01
        held['detectors']['particle_planes'] = [
            {'name': 'activation', 'z': 0, 'buffer_records': 64, 'compression': 0},
            {'name': 'downstream', 'z': 60, 'buffer_records': 64, 'compression': 0}]
        result['held-supplement'] = held
    return result


def physics(card):
    result = copy.deepcopy(card)
    result.pop('checkpoint', None)
    result.pop('output', None)
    for name in ('trajectory', 'detectors', 'energy_ledger'):
        if name in result:
            result[name].pop('directory', None)
    return result


def prepare(args):
    root = args.root.resolve()
    require(not root.exists(), 'Use a fresh result directory')
    require(args.ranks in (1, 2, 4), 'Use rank count 1, 2, or 4')
    binaries = {'accepted': args.baseline.resolve(), 'legacy': args.binary.resolve(),
                'cowan': args.binary.resolve()}
    for binary in binaries.values():
        require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable: ' + str(binary))
    launcher = json.loads(args.launcher_json)
    require(isinstance(launcher, list) and launcher and all(isinstance(v, str) for v in launcher)
            and any('{ranks}' in v for v in launcher), 'Launcher requires {ranks}')
    probe = args.binding_probe.resolve()
    require(probe.is_file(), 'Binding probe is required')
    root.mkdir(parents=True)
    (root / 'scripts').mkdir()
    shutil.copy2(HERE / 'run.py', root / 'scripts/run.py')
    shutil.copy2(HERE / 'README.txt', root / 'scripts/README.txt')
    if (HERE / 'analyze.py').is_file():
        shutil.copy2(HERE / 'analyze.py', root / 'scripts/analyze.py')
    shutil.copy2(HELPER, root / 'scripts/checkpoint_helpers.py')
    (root / 'inputs').mkdir()
    if args.include_held:
        common.held_input(root / 'inputs/held.h5')
    cases = cards(root, args.include_held)
    if args.cases:
        names = args.cases.split(',')
        require(set(names) <= cases.keys() and len(names) == len(set(names)), 'Unknown/duplicate case')
        cases = {name: cases[name] for name in names}
    plan = {'format': 'aprl-cowan-initial-integration-v1', 'ranks': args.ranks,
            'threads': 1, 'timeout_seconds': 180, 'timeout_total_seconds': 7200,
            'launcher': launcher, 'binding_probe': str(probe), 'variants': {}, 'cases': {},
            'scope': 'Small full-program integration, not a full-scale cost or radiation claim'}
    for name, binary in binaries.items():
        plan['variants'][name] = {'binary': str(binary), 'sha256': digest(binary)}
    for name, base in cases.items():
        steps = [0, 1, 7, 55] if name == 'held-supplement' else [0, 7, 55]
        for variant in VARIANTS:
            card = copy.deepcopy(base)
            card['initial_self_field']['model'] = NEW_MODEL if variant == 'cowan' else OLD_MODEL
            directory = root / name / variant
            written = common.write_card(directory, card, {
                'directory': str(directory / 'checkpoints'), 'steps': steps,
                'io_timeout_seconds': 60})
            normalized = physics(written)
            normalized['initial_self_field']['model'] = OLD_MODEL
            require(normalized == physics(base), 'A variant changed physics beyond initializer')
        plan['cases'][name] = {'base': base, 'steps': steps, 'restart_step': 7}
    save(root / 'plan.json', plan)
    frozen = list(root.rglob('*')) + list(set(binaries.values())) + [probe]
    for binary in set(binaries.values()):
        frozen.extend((binary.parent / 'runtime-libs').glob('*'))
    save(root / 'frozen.json', {str(p): digest(p) for p in sorted(set(frozen)) if p.is_file()})
    print(root / 'plan.json')


def verify(root):
    for name, expected in json.loads((root / 'frozen.json').read_text()).items():
        require(digest(Path(name)) == expected, 'Prepared artifact changed: ' + name)
    return json.loads((root / 'plan.json').read_text())


def idle_check():
    other = []
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            name = Path(os.fsdecode((entry / 'cmdline').read_bytes().split(b'\0')[0])).name
            if name in ('aprl', 'mpirun', 'mpiexec', 'orted'):
                other.append({'pid': int(entry.name), 'executable': name})
        except (OSError, ValueError):
            continue
    def counters():
        values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(values), values[3]
    before, idle_before = counters()
    time.sleep(2)
    after, idle_after = counters()
    fraction = (idle_after - idle_before) / (after - before)
    require(not other and fraction >= .95, 'Authorized CPU test host is busy: ' + str(other))
    return {'cpu_idle_fraction': fraction, 'compute_processes': other, 'loadavg': os.getloadavg()}


def binding(root, plan):
    directory = root / 'binding'
    directory.mkdir()
    command = [v.format(ranks=plan['ranks'], threads=1) for v in plan['launcher']]
    command += [plan['binding_probe'], str(plan['ranks']), '1']
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_PLACES='cores', OMP_PROC_BIND='close')
    completed = subprocess.run(['timeout', '--signal=TERM', '--kill-after=5', '45'] + command,
                               env=env, cwd=directory, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=55)
    (directory / 'aprl.log').write_text(completed.stdout)
    save(directory / 'launch.json', {'command': command, 'exit_code': completed.returncode})
    require(completed.returncode == 0 and '[binding] result=PASS failures=0' in completed.stdout,
            'CPU binding precheck failed')


def run(args):
    root = args.root.resolve()
    plan = verify(root)
    require(not (root / 'RUN_STARTED.json').exists(), 'Campaign already started')
    require(socket.gethostname().split('.')[0] == 'gnode1' and args.allow_gnode1 and
            os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1' and
            not os.environ.get('SLURM_JOB_ID'), 'Requires explicitly authorized gnode1 CPU session')
    lock = (WORKSPACE / 'results/node_survey/.cpu-calibration.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    save(root / 'initial-resources.json', idle_check())
    save(root / 'RUN_STARTED.json', {'host': socket.gethostname(), 'started_unix': time.time()})
    deadline = time.monotonic() + plan['timeout_total_seconds']
    binding(root, plan)
    for name, item in plan['cases'].items():
        for variant in VARIANTS:
            require(time.monotonic() + 180 < deadline, 'Campaign time budget exhausted')
            launch_plan = dict(plan, binary=plan['variants'][variant]['binary'])
            common.launch(root, launch_plan, root / name / variant)
        # A changed initializer must not be used to reinterpret a saved state.
        source = common.checkpoint_path(root / name / 'cowan', item['restart_step'])
        card = copy.deepcopy(item['base'])
        card['initial_self_field']['model'] = NEW_MODEL
        directory = root / name / 'restart'
        common.write_card(directory, card, {'restart_from': str(source),
                          'io_timeout_seconds': 60, 'reset_diagnostics': False})
        common.launch(root, dict(plan, binary=plan['variants']['cowan']['binary']), directory)
        if name == 'gaussian-moving-cpml':
            negative = root / name / 'restart-model-mismatch'
            common.write_card(negative, item['base'], {
                'restart_from': str(source), 'io_timeout_seconds': 60})
            common.launch(root, dict(plan, binary=plan['variants']['cowan']['binary']), negative,
                          rejection=r'Restart physics/diagnostic contract differs')
        # Stop before the next case if any completed case violates integration gates.
        review_case(root, plan, name)
    review(args)


def snapshot_compare(first, second, categories=('fields', 'particles', 'cpml', 'carriers')):
    checked = []
    with h5py.File(first) as a, h5py.File(second) as b:
        require(int(a['header'][2]) == int(b['header'][2]) == 1, 'Incomplete snapshot')
        for category in categories:
            prefix = 'blocks/' + category
            require((prefix in a) == (prefix in b), 'Checkpoint category mismatch')
            if prefix not in a:
                continue
            left, right = {}, {}
            for file, target in ((a, left), (b, right)):
                file[prefix].visititems(lambda n, obj: target.__setitem__(n, obj[()])
                                       if isinstance(obj, h5py.Dataset) else None)
            require(left.keys() == right.keys(), 'Checkpoint layout differs')
            for key in left:
                common.bitwise_equal(left[key], right[key], prefix + '/' + key)
                checked.append(prefix + '/' + key)
    return checked


def snapshot_particles(path):
    with h5py.File(path) as file:
        particles = file['blocks/particles/data'][:]
    return np.sort(particles, order='id')


def moments(particles):
    require(len(particles) > 0, 'Empty Gaussian cohort')
    weights = np.abs(particles['charge'])
    xyz = np.column_stack([particles[name] for name in ('x', 'y', 'z')])
    mean = np.average(xyz, weights=weights, axis=0)
    rms = np.sqrt(np.average((xyz - mean)**2, weights=weights, axis=0))
    gamma = np.sqrt(1 + sum(particles[name]**2 for name in ('ux', 'uy', 'uz')))
    mean_gamma = np.average(gamma, weights=weights)
    sigma_gamma = np.sqrt(np.average((gamma - mean_gamma)**2, weights=weights))
    return {'particles': len(particles), 'rms_m': rms.tolist(),
            'mean_kinetic_box_eV': float((mean_gamma - 1) * REST_EV),
            'sigma_energy_box_eV': float(sigma_gamma * REST_EV)}


def model_metrics(directory, steps):
    _, arrays, _ = common.inspect_file(directory / 'ledger/energy-ledger.h5')
    data = arrays['records']
    require(len(data) > 10, 'Too few ledger records')
    require(np.all(data['prescribed_work_J'] == 0), 'Free-drift case gained prescribed work')
    result = {'steps': data['step'].astype(int).tolist(),
              'time_box_fs': (data['time_box_s'] * 1e15).tolist(),
              'delta_particle_kinetic_J': (data['particle_kinetic_active_J'] - data['particle_kinetic_active_J'][0]).tolist(),
              'delta_mean_lab_energy_eV': ((data['mean_gamma_lab'] - data['mean_gamma_lab'][0]) * REST_EV).tolist(),
              'sigma_lab_energy_eV': (data['sigma_gamma_lab'] * REST_EV).tolist(),
              'sigma_box_energy_eV': (data['sigma_gamma_box'] * REST_EV).tolist(),
              'field_energy_interior_J': data['field_energy_interior_J'].tolist(),
              'balance_residual_J': data['balance_residual_J'].tolist(), 'snapshot_moments': {}}
    initial = snapshot_particles(common.checkpoint_path(directory, 0))
    config = yaml.safe_load((directory / 'run.yaml').read_text())
    boost = np.longdouble(config['mesh']['boost_gamma'])
    beta = np.sqrt(1 - 1 / boost**2)
    u0 = np.column_stack([initial[key] for key in ('ux', 'uy', 'uz')]).astype(np.longdouble)
    g0 = np.sqrt(1 + np.sum(u0**2, axis=1))
    for step in steps:
        particles = snapshot_particles(common.checkpoint_path(directory, step))
        value = moments(particles)
        common.bitwise_equal(initial['id'], particles['id'], 'Snapshot particle cohort')
        u = np.column_stack([particles[key] for key in ('ux', 'uy', 'uz')]).astype(np.longdouble)
        gamma = np.sqrt(1 + np.sum(u**2, axis=1))
        # Stable per-particle change; do not subtract two rounded beam totals.
        delta_gamma = np.sum((u-u0)*(u+u0), axis=1) / (gamma+g0)
        weights = np.abs(initial['charge']).astype(np.longdouble)
        value['mean_delta_box_energy_eV'] = float(np.average(delta_gamma, weights=weights) * REST_EV)
        value['mean_delta_lab_energy_eV'] = float(np.average(
            boost * (delta_gamma + beta*(u[:, 2]-u0[:, 2])), weights=weights) * REST_EV)
        result['snapshot_moments'][str(step)] = value
    return result


def review_case(root, plan, name):
    item = plan['cases'][name]
    directory = root / name
    for variant in VARIANTS:
        launch_record = json.loads((directory / variant / 'launch.json').read_text())
        require(launch_record['exit_code'] == 0, 'Simulation did not exit successfully')
        log = (directory / variant / 'aprl.log').read_text()
        require('Configured stop reached:' in log, 'Missing physical completion')
        manifest = yaml.safe_load((directory / variant / 'run-manifest.yaml').read_text())
        require(manifest['mpi_size'] == plan['ranks'], 'Wrong MPI count')
        require(yaml.safe_load(manifest['configuration_yaml']) ==
                yaml.safe_load((directory / variant / 'run.yaml').read_text()), 'Executed card differs')
        marker = ('Cowan relativistic-Poisson quasi-steady rigid-beam construction converged'
                  if variant == 'cowan' else 'relativistic-Poisson rigid-beam construction converged')
        require(marker in log, 'Requested initializer did not converge: ' + variant)
        for path in common.output_files(directory / variant).values():
            common.inspect_file(path)
    exact = common.compare_outputs(directory / 'accepted', directory / 'legacy')
    checked = {}
    for step in item['steps']:
        checked[str(step)] = snapshot_compare(common.checkpoint_path(directory / 'accepted', step),
                                              common.checkpoint_path(directory / 'legacy', step))
    # Initializer is allowed to change fields, never initial particle state.
    initial_equal = snapshot_compare(common.checkpoint_path(directory / 'legacy', 0),
                                     common.checkpoint_path(directory / 'cowan', 0), ('particles',))
    source = common.checkpoint_path(directory / 'cowan', item['restart_step'])
    require(json.loads((directory / 'restart/launch.json').read_text())['exit_code'] == 0,
            'Restart did not exit successfully')
    restarted = common.compare_outputs(directory / 'cowan', directory / 'restart',
                                       item['restart_step'], source)
    if name == 'gaussian-moving-cpml':
        rejected = directory / 'restart-model-mismatch'
        rejected_launch = json.loads((rejected / 'launch.json').read_text())
        require(rejected_launch['exit_code'] not in (0, 124, 137, 143) and
                'Restart physics/diagnostic contract differs' in (rejected / 'aprl.log').read_text(),
                'Changed-model restart was not explicitly rejected')
    coverage = {str(step): common.inspect_snapshot(common.checkpoint_path(directory / 'cowan', step))
                for step in item['steps']}
    if 'cpml' in name or name == 'held-supplement':
        require(coverage['55']['cpml_nonzero'] > 0, 'No nonzero CPML history covered')
    if name == 'held-supplement':
        log = (directory / 'cowan/aprl.log').read_text()
        match = re.search(r'Frozen lab-plane injection: activations=(\d+), '
                          r'cross-rank activation pushes=(\d+), still frozen=(\d+)', log)
        require(match and int(match[1]) == 64 and int(match[3]) == 0,
                'Supplemental held cohort did not fully activate')
        require(coverage['1']['held'] > 0, 'Early snapshot did not preserve any held particle')
        coverage['held_activations'] = int(match[1])
        coverage['cross_rank_activation_pushes'] = int(match[2])
    metrics = {variant: model_metrics(directory / variant, item['steps'])
               for variant in ('legacy', 'cowan')}
    orbit = {}
    for step in item['steps']:
        a = snapshot_particles(common.checkpoint_path(directory / 'legacy', step))
        b = snapshot_particles(common.checkpoint_path(directory / 'cowan', step))
        common.bitwise_equal(a['id'], b['id'], 'Model-comparison particle cohort')
        orbit[str(step)] = {component: relative_difference(a[component], b[component])
                            for component in ('x', 'y', 'z', 'ux', 'uy', 'uz')}
    _, field_data, _ = common.inspect_file(directory / 'cowan/detectors/field.h5')
    require(len(field_data['time_s']) > 10 and any(np.any(field_data['fields'][key] != 0)
            for key in field_data['fields'].dtype.names), 'Field-plane sampling is empty or identically zero')
    coverage['field_samples'] = len(field_data['time_s'])
    report = {'integration_status': 'PASS', 'exact_legacy_datasets': exact,
              'exact_legacy_snapshot_datasets': checked, 'initial_particle_state_exact': initial_equal,
              'exact_restart_suffix_datasets': restarted, 'coverage': coverage, 'metrics': metrics,
              'new_vs_legacy_particle_differences': orbit,
              'physics_status': 'MEASURED; approximation differences require scientific interpretation'}
    save(directory / 'REVIEW.json', report)
    return report


def review(args):
    root = args.root.resolve()
    plan = verify(root)
    reports = {name: review_case(root, plan, name) for name in plan['cases']}
    result = {'integration_status': 'PASS', 'ranks': plan['ranks'], 'cases': reports,
              'limitations': ['Same-rank legacy and restart checks are bitwise.',
                             'New versus legacy initializer physics is measured, not required bitwise.',
                             'No full-scale speed, pure-radiation or entire-production-window claim.']}
    save(root / 'RESULT.json', result)
    (root / 'PASS').write_text('Full-program integration gates passed; see measured physical differences.\n')
    print(root / 'RESULT.json')


def relative_difference(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    require(a.shape == b.shape, 'Cross-rank output shape differs')
    scale = np.linalg.norm(a.ravel())
    error = np.linalg.norm((b-a).ravel())
    return {'rms_relative': float(error / scale) if scale else (0.0 if error == 0 else None),
            'max_abs': float(np.max(np.abs(b-a))) if a.size else 0.0,
            'reference_l2': float(scale), 'difference_l2': float(error)}


def aggregate(args):
    roots = [p.resolve() for p in args.roots]
    plans = [verify(p) for p in roots]
    require(plans[0]['ranks'] == 1, 'First reference must be the one-rank campaign')
    for root, plan in zip(roots, plans):
        require(json.loads((root / 'RESULT.json').read_text())['integration_status'] == 'PASS',
                'Cross-rank comparison requires each integration campaign to pass')
        require(plan['variants'] == plans[0]['variants'], 'Cross-rank executables differ')
    result = {'reference': str(roots[0]), 'comparisons': [],
              'interpretation': 'Report new-model cross-rank differences alongside unchanged legacy differences; no relaxed tolerance or speed claim.'}
    for root, plan in zip(roots[1:], plans[1:]):
        require(plan['cases'].keys() == plans[0]['cases'].keys(), 'Campaign case sets differ')
        for name in plan['cases']:
            a_card = physics(plans[0]['cases'][name]['base'])
            b_card = physics(plan['cases'][name]['base'])
            if name == 'held-supplement':
                a_card['beam']['input']['file'] = b_card['beam']['input']['file'] = '<same-fixed-fixture>'
                with h5py.File(roots[0] / 'inputs/held.h5') as first, h5py.File(root / 'inputs/held.h5') as second:
                    common.bitwise_equal(first['particles/records'][:], second['particles/records'][:],
                                         'Cross-campaign fixed held fixture')
            require(a_card == b_card, 'Cross-rank physics differs')
            report = {'case': name, 'ranks': plan['ranks'], 'variants': {},
                      'initial_particles_by_id_bitwise': {}}
            for variant in VARIANTS:
                first = snapshot_particles(common.checkpoint_path(roots[0] / name / variant, 0))
                second = snapshot_particles(common.checkpoint_path(root / name / variant, 0))
                common.bitwise_equal(first, second, 'Cross-rank initial particle state: '+name+'/'+variant)
                report['initial_particles_by_id_bitwise'][variant] = {
                    'particles': len(first), 'fields': list(first.dtype.names), 'status': 'PASS'}
            for variant in ('legacy', 'cowan'):
                _, a, _ = common.inspect_file(roots[0] / name / variant / 'detectors/field.h5')
                _, b, _ = common.inspect_file(root / name / variant / 'detectors/field.h5')
                common.bitwise_equal(a['time_s'], b['time_s'], 'Cross-rank sampling times')
                groups = {component: relative_difference(a['fields'][component], b['fields'][component])
                          for component in a['fields'].dtype.names}
                components = {}
                for group, prefix in (('electric_V_per_m', 'E'), ('magnetic_T', 'B')):
                    require(group in a['fields'].dtype.names and a['fields'][group].shape[-1] == 3,
                            'Unexpected field-vector layout')
                    for axis, suffix in enumerate(('x', 'y', 'z')):
                        components[prefix+suffix] = relative_difference(
                            a['fields'][group][..., axis], b['fields'][group][..., axis])
                report['variants'][variant] = {'field_groups': groups, 'field_components': components}
            result['comparisons'].append(report)
    save(args.output.resolve(), result)
    print(args.output.resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--root', type=Path, required=True)
    prep.add_argument('--binary', type=Path, default=WORKSPACE / 'build/cowan-integration-20261027/aprl')
    prep.add_argument('--baseline', type=Path, default=WORKSPACE / 'build/cowan-integration-control-20261027/aprl')
    prep.add_argument('--binding-probe', type=Path, default=WORKSPACE / 'build/node-survey/binding_probe')
    prep.add_argument('--ranks', type=int, required=True)
    prep.add_argument('--cases', default='')
    prep.add_argument('--include-held', action='store_true')
    prep.add_argument('--launcher-json', default='["mpirun","--host","gnode1:16","-np","{ranks}","--map-by","slot:PE=1","--bind-to","core","--nooversubscribe","--report-bindings"]')
    for action in ('run', 'review'):
        child = sub.add_parser(action)
        child.add_argument('--root', type=Path, required=True)
        child.add_argument('--allow-gnode1', action='store_true')
    agg = sub.add_parser('aggregate')
    agg.add_argument('--roots', type=Path, nargs='+', required=True)
    agg.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        globals()[args.action](args)
    except Exception as error:
        if getattr(args, 'root', None) and args.root.exists() and args.action in ('run', 'review'):
            save(args.root / 'FAIL.json', {'error': str(error), 'time_unix': time.time()})
        raise


if __name__ == '__main__':
    main()
