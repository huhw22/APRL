#!/usr/bin/env python3
"""Freeze/check four small 256-rank held/retirement shutdown acceptance cases."""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil

import h5py
import numpy as np
import yaml


HERE = Path(__file__).resolve().parent
HELPER = HERE / 'science_helpers.py'
if not HELPER.is_file():
    HELPER = HERE.parent / 'full_scale_cowan/run.py'
spec = importlib.util.spec_from_file_location('duration_acceptance_science', HELPER)
science = importlib.util.module_from_spec(spec)
spec.loader.exec_module(science)
require, digest, save = science.require, science.digest, science.save
WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
VARIANTS = ('normal-control', 'normal-fixed', 'bounded-1', 'bounded-2')
STEPS = 200
CONTROL_SHA = '56cc1e0931103687dcc440c3ad390a923b2e74f5da0d68ef1806e3992ccbfeac'
CONTROL_REVISION = 'ac58e67e33d6b4ecfd12c3d00f4a9947d0e8500d'


def card(path):
    return yaml.safe_load(path.read_text())


def physical_configuration(configuration):
    result = copy.deepcopy(configuration)
    for section, key in (('output', 'manifest'), ('detectors', 'directory'),
                         ('energy_ledger', 'directory'), ('mesh', 'duration')):
        del result[section][key]
    return result


def prepare(root, binary, revision):
    require(re.fullmatch(r'[0-9a-f]{40}', revision), 'Use the exact compiled source revision')
    reference = WORKSPACE / 'results/elegant_retirement_20261024/compact-observed-window/preflight-r2'
    require((reference / 'PASS').is_file(), 'Accepted small held/retirement reference is missing')
    original = card(reference / 'functional/run.yaml')
    require(original['units'] == {'length': 'micrometer', 'time': 'femtosecond'}, 'Unexpected card units')
    require(original['mesh']['field_solver'] == 'cowan-z', 'Cowan dt=dz/c required')
    control = WORKSPACE / 'build/candidate/aprl'
    require(digest(control) == CONTROL_SHA, 'Accepted control binary changed')
    tools = {'control': science.build_metadata(control), 'fixed': science.build_metadata(binary)}
    require(tools['control']['selected'] == tools['fixed']['selected'], 'Control/fixed toolchains differ')
    require(not root.exists(), 'Use a new acceptance directory')
    root.mkdir(parents=True)
    scripts = root / 'scripts'
    scripts.mkdir()
    shutil.copy2(Path(__file__), scripts / 'cases.py')
    shutil.copy2(HELPER, scripts / 'science_helpers.py')
    shutil.copy2(reference / 'functional/particles.h5', root / 'particles.h5')
    shutil.copy2(reference / 'functional/run.yaml', root / 'reference-run.yaml')
    shutil.copy2(reference / 'functional/run-manifest.yaml', root / 'reference-manifest.yaml')
    original['beam']['input']['file'] = str(root / 'particles.h5')
    original['detectors']['field_planes'][0]['name'] = 'duration-field'
    # An empty writer must also flush its HDF5 metadata before expected abort.
    original['detectors']['particle_planes'].append({
        'name': 'after-stop-empty', 'z': 2000, 'buffer_records': 512, 'compression': 0})
    dt = original['mesh']['cell_size'][2] * 1e-6 / 299792458.0
    variants = {}
    for name in VARIANTS:
        directory = root / name
        directory.mkdir()
        configuration = science.redirect(original, directory)
        bounded = name.startswith('bounded-')
        if bounded:
            configuration['mesh']['duration'] = (STEPS - 0.5) * dt / 1e-15
        require(physical_configuration(configuration) == physical_configuration(original),
                'Unexpected physics/diagnostics difference between cases')
        (directory / 'run.yaml').write_text(yaml.safe_dump(configuration, sort_keys=False))
        chosen = control if name == 'normal-control' else binary
        variants[name] = {'binary': str(chosen), 'binary_sha256': digest(chosen),
                          'source_revision': CONTROL_REVISION if name == 'normal-control' else revision,
                          'bounded': bounded, 'expected_exit_code': 1 if bounded else 0,
                          'delay_h5_basename': 'duration-field.h5' if bounded else None}
    plan = {'format': 'aprl-duration-acceptance-v1', 'variants': variants,
            'nodes': 2, 'ranks': 256, 'threads_per_rank': 1, 'particles': 512,
            'expected_mesh': [16, 16, 4096], 'bounded_steps': STEPS, 'dt_box_s': dt,
            'input_sha256': digest(root / 'particles.h5'), 'toolchains': tools,
            'comparison': 'Same-rank scientific datasets and attributes bitwise; no relaxed tolerances',
            'scope': 'Small held/retirement and deterministic delayed-close correctness, not performance or radiation convergence'}
    save(root / 'plan.json', plan)
    frozen = sorted(p for p in root.rglob('*') if p.is_file()) + [binary, control]
    frozen += [Path(v['cache']) for v in tools.values()]
    (root / 'frozen.sha256').write_text(''.join(f'{digest(p)}  {p}\n' for p in frozen))
    print(root)


def verify(root):
    for line in (root / 'frozen.sha256').read_text().splitlines():
        sha, filename = line.split('  ', 1)
        require(digest(filename) == sha, 'Frozen acceptance evidence changed: ' + filename)


def check_variant(root, name):
    plan = json.loads((root / 'plan.json').read_text())
    entry = plan['variants'][name]
    bounded = entry['bounded']
    directory = root / name
    require(int((directory / 'EXIT_CODE').read_text()) == entry['expected_exit_code'],
            'Unexpected simulator exit: ' + name)
    log = (directory / 'aprl.log').read_text(errors='replace')
    if bounded:
        require(science.EXHAUSTED in log, 'Missing expected bounded-duration exception')
        require('phase=duration-exhaustion-readiness-sync-enter' in log,
                'Duration close-readiness barrier was not exercised')
        for marker in ('DURATION_CLOSE_DELAY_BEGIN file=duration-field.h5 seconds=3',
                       'DURATION_CLOSE_DELAY_COMPLETE file=duration-field.h5'):
            require(marker in log, 'Deterministic delayed HDF5 close was not exercised: ' + marker)
    else:
        require('Configured stop reached:' in log and 'Simulation failure' not in log,
                'Normal physical stop did not complete')
    # Post-barrier exception lines can be interrupted by the first MPI_Abort.
    # Root close/final reports and the complete HDF5 checks are mandatory.
    for phase in ('detector-close-complete', 'ledger-close-complete', 'final-resource-report-complete'):
        require(re.search(r'\[teardown\] rank=0\b[^\n]* phase=' + phase + r'\b', log),
                'Missing root scientific close marker: ' + phase)
    final = re.search(r'\[resource\] phase=final step=(\d+)\b', log)
    require(final is not None, 'Missing root final resource report')
    final_step = int(final.group(1))
    require(final_step == STEPS if bounded else final_step > STEPS, 'Unexpected final step')
    if bounded:
        progress = [int(x) for x in re.findall(r'\[resource\] phase=progress step=(\d+)\b', log)]
        require(progress == list(range(50, STEPS + 1, 50)), 'Missing bounded progress steps')
    manifest = card(directory / 'run-manifest.yaml')
    configuration = card(directory / 'run.yaml')
    require(manifest['mpi_size'] == 256 and manifest['source_revision'] == entry['source_revision'],
            'Incorrect runtime MPI size or compiled revision')
    require(manifest['resolved_mesh']['cells'] == plan['expected_mesh'], 'Resolved mesh changed')
    require(manifest['resolved_beam_reference']['startup_mode'] == 'held', 'Held startup changed')
    require(yaml.safe_load(manifest['configuration_yaml']) == configuration, 'Manifest/card mismatch')
    require(digest(entry['binary']) == entry['binary_sha256'], 'Executable changed')
    outputs = science.inspect_outputs(directory, not bounded)
    with h5py.File(directory / 'ledger/energy-ledger.h5') as file:
        records = file['energy_ledger/records'][:]
        require(np.array_equal(records['step'], np.arange(final_step + 1)), 'Ledger step sequence incomplete')
        expected_time = np.longdouble(final_step) * np.longdouble(plan['dt_box_s'])
        unit = np.longdouble(np.finfo(float).eps) / 2
        bound = (final_step * unit / (1 - final_step * unit)) * abs(expected_time)
        require(abs(np.longdouble(records[-1]['time_box_s']) - expected_time) <= bound,
                'Ledger final time differs beyond accumulated-step rounding')
        require(int(records[-1]['removed_macroparticles']) == 512 and
                int(records[-1]['active_macroparticles']) == 0, 'Retirement cohort incomplete')
        empty = np.flatnonzero(records['active_macroparticles'] == 0)
        require(len(empty) and final_step + 1 - int(empty[0]) >= 20,
                'Insufficient continued field propagation after retirement')
    with h5py.File(directory / 'detectors/after-stop-empty.h5') as file:
        require(len(file['particle_plane/records']) == 0, 'Expected empty particle plane is nonempty')
    with h5py.File(root / 'particles.h5') as file:
        source = file['particles/records'][:]
    for detector in ('activation', 'retirement-entry'):
        with h5py.File(directory / 'detectors' / (detector + '.h5')) as file:
            observed = file['particle_plane/records'][:]
        observed = observed[np.argsort(observed['source_id'])]
        require(len(observed) == 512 and np.array_equal(observed['source_id'], source['source_id']),
                'Particle cohort differs at ' + detector)
        if detector == 'activation':
            require(np.max(np.abs(observed['position_m'][:, :2] - source['plane_position_m'])) <= 1e-12,
                    'Held activation transverse positions differ from reference check')
    with h5py.File(directory / 'detectors/duration-field.h5') as file:
        require(len(file['field_plane/time_s']) > 20, 'Missing field propagation samples')
    report = {'pass': True, 'variant': name, 'bounded': bounded, 'steps': final_step,
              'outputs': outputs, 'exit_code': entry['expected_exit_code'],
              'deterministic_close_delay_exercised': bounded}
    save(directory / 'acceptance-result.json', report)
    return report


def compare_prefix(short, full):
    """Compare bounded physics to the exact matching prefix of a normal run."""
    compared = []
    ignored_attributes = {'run_id', 'run_manifest', 'configuration_digest_fnv1a64', 'source_revision'}
    for relative in science.output_files(short):
        with h5py.File(short / relative) as a, h5py.File(full / relative) as b:
            require(science.dataset_names(a) == science.dataset_names(b), 'Prefix dataset schemas differ')
            objects = ['']
            a.visit(objects.append)
            for name in objects:
                aa, bb = (a, b) if not name else (a[name], b[name])
                keys = set(aa.attrs) - ignored_attributes
                require(keys == set(bb.attrs) - ignored_attributes, 'Prefix scientific attributes differ')
                for key in keys:
                    science.bitwise(np.asarray(aa.attrs[key]), np.asarray(bb.attrs[key]), str(relative) + '@' + key)
            for name in science.dataset_names(a):
                if name.rsplit('/', 1)[-1] in ('complete', 'committed_samples', 'committed_records'):
                    continue
                aa, bb = a[name], b[name]
                require(aa.dtype == bb.dtype and aa.shape[1:] == bb.shape[1:] and
                        aa.shape[0] <= bb.shape[0], 'Scientific prefix layout differs')
                for selection in science.slices(aa):
                    science.bitwise(np.asarray(aa[selection]), np.asarray(bb[selection]), str(relative) + '/' + name)
                compared.append(str(relative) + '/' + name)
    return {'pass': True, 'datasets': compared, 'comparison': 'bounded datasets exactly equal normal scientific prefixes'}


def check(root, variant=None):
    verify(root)
    if variant:
        print(json.dumps(check_variant(root, variant), indent=2))
        return
    reports = {name: check_variant(root, name) for name in VARIANTS}
    comparisons = {
        'normal_control_vs_fixed': science.compare_outputs(root / 'normal-control', root / 'normal-fixed'),
        'bounded_repeat': science.compare_outputs(root / 'bounded-1', root / 'bounded-2'),
        'bounded_vs_normal_prefix': compare_prefix(root / 'bounded-1', root / 'normal-fixed')}
    save(root / 'acceptance-result.json', {'pass': True, 'cases': reports, 'comparisons': comparisons,
                                         'ranks': 256, 'threads_per_rank': 1,
                                         'scope': 'Small held/retirement functionality and complete scientific shutdown; not a performance run'})
    (root / 'PASS').write_text('Normal same-rank equivalence, bounded repeats/prefixes, and all scientific output closes passed.\n')
    print(root / 'PASS')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'check'))
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--revision')
    parser.add_argument('--variant', choices=VARIANTS)
    args = parser.parse_args()
    if args.mode == 'prepare':
        if args.binary is None or args.revision is None:
            parser.error('prepare requires --binary and --revision')
        prepare(args.root.resolve(), args.binary.resolve(), args.revision)
    else:
        check(args.root.resolve(), args.variant)
