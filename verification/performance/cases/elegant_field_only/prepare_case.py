#!/usr/bin/env python3
"""Freeze the authorized held, field-only ELEGANT case; never submit jobs."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

import h5py
import yaml


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def prepare_free_drift(workspace, root, smoke, reference):
    """Freeze a no-magnet background case from an identified held production."""
    if reference is None:
        raise ValueError('--free-drift requires --reference')
    reference = reference.resolve()
    parent = json.loads((reference / 'run-plan.json').read_text())
    original = yaml.safe_load((reference / 'run.yaml').read_text())
    manifest = yaml.safe_load((reference / 'run-manifest.yaml').read_text())
    if yaml.safe_load(manifest['configuration_yaml']) != original:
        raise ValueError('Reference manifest/configuration mismatch')
    if not json.loads((reference / 'geometry-preflight.json').read_text())['pass']:
        raise ValueError('Reference geometry has not passed')
    if manifest['resolved_beam_reference']['startup_mode'] != 'held':
        raise ValueError('Reference must use held startup')
    fields = original['detectors']['field_planes']
    magnets = original['sources']['magnetic_elements']
    if (len(fields) != 1 or fields[0]['z'] != 4.5 or len(magnets) != 1
            or magnets[0]['type'] != 'planar-undulator'
            or original['sources']['incident_waves']
            or original['detectors']['particle_planes']
            or original['particle_retirement']['enabled']
            or fields[0].get('particle_background', {}).get('enabled', False)):
        raise ValueError('Expected the authorized 4.5m held, raw field-only reference')
    if original['stop']['mode'] != 'after-last-element':
        raise ValueError('Expected particle-state physical endpoint')
    accepted = json.loads((workspace / 'results/detector_sampling_20261027/PROVENANCE.json').read_text())
    binary = Path(parent['binary'])
    binary_sha = digest(binary)
    if binary_sha != parent['binary_sha256'] or binary_sha != accepted['binaries_sha256'][str(binary)]:
        raise ValueError('Reference executable differs from accepted MPI sampling build')
    smoke_plan = json.loads((smoke / 'plan.json').read_text())
    if smoke_plan['binary_sha256'] != binary_sha or smoke_plan.get('case_variant') != 'free-drift':
        raise ValueError('Require free-drift smoke acceptance for this executable')
    input_path = Path(parent['input_file'])
    stat = input_path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (parent['input_size_bytes'], parent['input_mtime_ns']):
        raise ValueError('Full particle input changed')

    card = copy.deepcopy(original)
    card['output'] = {'overwrite': False, 'manifest': str(root / 'run-manifest.yaml')}
    card['sources'] = {'incident_waves': [], 'magnetic_elements': []}
    card['detectors']['directory'] = str(root / 'detectors')
    card['detectors']['field_planes'][0]['z'] = 0.5
    card['energy_ledger']['directory'] = str(root / 'ledger')
    card.pop('checkpoint', None)
    # Keep every input, solver, boost, boundary, held, resolution and guard
    # parameter unchanged. The core naturally selects one particle substep
    # with no magnetic period; it also resolves a new compact longitudinal box.
    comparable = copy.deepcopy(card)
    for key in ('output', 'sources', 'detectors', 'energy_ledger', 'checkpoint'):
        comparable.pop(key, None)
    parent_comparable = copy.deepcopy(original)
    for key in ('output', 'sources', 'detectors', 'energy_ledger', 'checkpoint'):
        parent_comparable.pop(key, None)
    if comparable != parent_comparable:
        raise ValueError('Unintended reference physics change')
    nominal_exit = magnets[0]['entrance_z'] + magnets[0]['period'] * magnets[0]['periods']
    comparison = {
        'reference_root': str(reference), 'reference_config_sha256': digest(reference / 'run.yaml'),
        'reference_run_id': manifest['run_id'], 'reference_field_z_m': 4.5,
        'reference_nominal_magnet_exit_z_m': nominal_exit,
        'reference_distance_from_nominal_exit_m': 4.5 - nominal_exit,
        'control_field_z_m': 0.5, 'control_magnetic_elements': [],
        'same_input_binary_boost_resolution_held_self_field_and_guard_parameters': True,
        'reference_resolved_grid': manifest['resolved_mesh']['cells'],
        'reference_particle_substeps': 5, 'control_particle_substeps': 1,
        'interpretation': 'Free-drift near-field/startup background comparison, not a certified radiation-free solution.',
        'limitations': [
            'The requested0.5m is the approximate downstream distance of the4m device; the reference nominal exit is3.978m.',
            'Automatic longitudinal bounds change with the detector position, including the finite-domain initial Poisson boundary.',
            'With no magnetic period the existing core uses one Boris substep; the field/current timestep and detector cadence are unchanged.',
            'Held activation and self-consistent particle motion can produce field transients; original input modulation remains.',
            'These runs do not supply an exact common-trajectory near-field subtraction; their raw fields cannot simply be identified as radiation after subtraction.',
        ],
    }
    root.mkdir(parents=True)
    scripts = root / 'scripts'; scripts.mkdir()
    for file in sorted(Path(__file__).resolve().parent.iterdir()):
        if file.is_file():
            shutil.copy2(file, scripts / file.name)
    (root / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
    save(root / 'reference-comparison.json', comparison)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=workspace / 'source/aprl', text=True).strip()
    plan = copy.deepcopy(parent)
    plan.update({
        'case_variant': 'free-drift', 'preparation_revision': revision,
        'reference_root': str(reference), 'reference_comparison': str(root / 'reference-comparison.json'),
        'preflight_result': str(smoke), 'geometry_report': str(root / 'geometry-preflight.json'),
        'field_plane_z_m': 0.5, 'checkpoint_steps': [], 'checkpoint_directory': str(root / 'checkpoints'),
        'expected_particle_substeps': 1,
        'physics': 'Full original beam, held, no magnetic elements or incident waves, self field/feedback, no retirement; raw lab field at0.5m.',
        'changes': ['Remove every magnetic element and incident source.',
                    'Move the single raw lab field plane from4.5m to user-specified0.5m.',
                    'Disable the reference front-exit checkpoint, which has no magnet-exit meaning in this control.',
                    'Retain original compact longitudinal policy; actual bounds are resolved for0.5m.',
                    'Existing automatic particle subcycling resolves to1 with no magnetic period.'],
        'radiation_free_claim': False,
    })
    save(root / 'run-plan.json', plan)
    frozen = [root / 'run.yaml', root / 'run-plan.json', root / 'reference-comparison.json'] + sorted(p for p in scripts.iterdir() if p.is_file())
    (root / 'frozen.sha256').write_text(''.join(f'{digest(p)}  {p.relative_to(root)}\n' for p in frozen))
    print(root / 'run.yaml')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--smoke', type=Path, required=True)
    parser.add_argument('--free-drift', action='store_true', help='No magnetic elements; raw field at0.5m, no checkpoint')
    parser.add_argument('--reference', type=Path, help='Frozen4.5m held field-only production to compare')
    args = parser.parse_args()
    w, root, smoke = args.workspace.resolve(), args.root.resolve(), args.smoke.resolve()
    if root.exists() or not (smoke / 'PASS').is_file():
        raise ValueError('Use a fresh production directory after small held field-only acceptance')
    if args.free_drift:
        prepare_free_drift(w, root, smoke, args.reference)
        return
    if args.reference is not None:
        raise ValueError('--reference is only used with --free-drift')
    old = w / 'results/elegant_retirement_20261024/compact-observed-window/production'
    historical = w / 'results/elegant_4m_analysis_20261024/original/runs/20261021T070618Z-861230'
    old_plan = json.loads((old / 'run-plan.json').read_text())
    accepted = json.loads((w / 'results/detector_sampling_20261027/PROVENANCE.json').read_text())
    binary = w / 'build/detector-sampling-20261027/aprl'
    binary_sha = digest(binary)
    if binary_sha != accepted['binaries_sha256'][str(binary)]:
        raise ValueError('Latest executable no longer matches sampling acceptance')
    smoke_plan = json.loads((smoke / 'plan.json').read_text())
    if smoke_plan['binary_sha256'] != binary_sha:
        raise ValueError('Smoke tested a different executable')
    input_path = Path(old_plan['input_file'])
    stat = input_path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (old_plan['input_size_bytes'], old_plan['input_mtime_ns']):
        raise ValueError('Full particle input changed; fresh audit required')

    card = copy.deepcopy(yaml.safe_load((old / 'run.yaml').read_text()))
    card['output'] = {'overwrite': False, 'manifest': str(root / 'run-manifest.yaml')}
    card['beam']['input']['file'] = str(input_path)
    card['detectors'] = {'enabled': True, 'directory': str(root / 'detectors'),
        'field_planes': [{'name': 'exit-field', 'z': 4.5,
                          'rhythm': card['detectors']['field_planes'][0]['rhythm'],
                          'buffer_samples': 8, 'compression': 0,
                          'particle_background': {'enabled': False}}],
        'particle_planes': []}
    card['particle_retirement'] = {'enabled': False}
    card['stop'] = {'mode': 'after-last-element'}
    card['energy_ledger']['directory'] = str(root / 'ledger')
    card['trajectory'] = {'enabled': False}

    # A measured leading-particle event from the same held beam, prior to any
    # retirement action. The user permits an approximately equivalent time;
    # this avoids a new full-particle reduction at every production step.
    source = old / 'detectors/pre-retirement-particles.h5'
    geometry = json.loads((old / 'geometry-preflight.json').read_text())
    with h5py.File(source, 'r') as file:
        group = file['particle_plane']
        if int(group['complete'][()]) != 1:
            raise ValueError('Historical crossing output incomplete')
        prefix = group['records'][:65536]
        first_t = float(prefix['time_s'].min())
        crossing_z = float(group.attrs['plane_z_m'])
    gamma = card['mesh']['boost_gamma']
    beta = math.sqrt(1 - 1 / gamma ** 2)
    t0 = geometry['held_geometry']['frame_origin_lab_time_s']
    z0 = geometry['held_geometry']['frame_origin_lab_z_m']
    dt = card['mesh']['cell_size'][2] / 299792458.
    event_box = gamma * (first_t - t0 - beta * (crossing_z - z0) / 299792458.)
    checkpoint_step = math.ceil(event_box / dt)
    if checkpoint_step != 3466 or not 4.0 < crossing_z < 4.5:
        raise ValueError('Unexpected front-exit checkpoint estimate')
    card['checkpoint'] = {'directory': str(root / 'checkpoints'),
                          'steps': [checkpoint_step], 'io_timeout_seconds': 3600}
    root.mkdir(parents=True)
    scripts = root / 'scripts'
    scripts.mkdir()
    here = Path(__file__).resolve().parent
    for file in sorted(here.iterdir()):
        if file.is_file():
            shutil.copy2(file, scripts / file.name)
    (root / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
    snapshot = {
        'kind': 'one complete synchronized state, scheduled by accepted step API',
        'step': checkpoint_step, 'time_box_s': checkpoint_step * dt,
        'path': str(root / 'checkpoints' / f'step-{checkpoint_step:012d}.h5'),
        'basis_source': str(source), 'basis_prefix_records': len(prefix),
        'basis_plane_z_m': crossing_z, 'basis_first_event_lab_s': first_t,
        'basis_first_event_box_s': event_box,
        'expected_front_z_m_approximately': 4.03,
        'field_plane_z_m': 4.5, 'approximate_remaining_distance_m': .47,
        'interpretation': 'Leading particles have exited; trailing particles may remain in the undulator or held. This is a common boost-time slice.',
        'limitation': 'Step selected from prior same-held trajectories, not a new exact crossing trigger. Verify actual front position from the saved state.',
        'future_particle_diagnostics': 'Same256ranks; new output directory; checkpoint.reset_diagnostics=true permits new downstream particle planes. Crossings before the snapshot cannot be recreated.'}
    save(root / 'snapshot-request.json', snapshot)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=w / 'source/aprl', text=True).strip()
    plan = {
        'format': 'aprl-held-field-only-production-v1',
        'binary': str(binary), 'binary_sha256': binary_sha,
        'accepted_source_revision': accepted['commit'],
        'compiled_revision': accepted['compiled_revision'], 'preparation_revision': revision,
        'input_file': str(input_path), 'input_size_bytes': stat.st_size,
        'input_mtime_ns': stat.st_mtime_ns, 'input_sha256': old_plan['input_sha256'],
        'records': 41732775, 'electrons': 41732775,
        'preflight_result': str(smoke), 'historical_plan': str(historical / 'input-plan.json'),
        'geometry_report': str(root / 'geometry-preflight.json'),
        'nodes': 2, 'ranks': 256, 'threads_per_rank': 1, 'exclusive': True,
        'time_limit': '144:00:00', 'field_plane_name': 'exit-field', 'field_plane_z_m': 4.5,
        'checkpoint_steps': [checkpoint_step], 'checkpoint_directory': str(root / 'checkpoints'),
        'physics': 'Full original beam, held startup, Cowan, self field/feedback, no retirement; raw field plane only.',
        'stop_policy': 'after-last-element: surviving physical particles pass the field location; independent of particle-plane export.',
        'quiet_tail_claim': False,
        'changes': ['Field plane fixed at user-specified4.5m.', 'All particle-plane/reference outputs disabled.',
                    'Retirement and retirement frequency protection disabled.',
                    'Original compact64/128guards and256-cell alignment retained.',
                    'One complete snapshot near leading-particle undulator exit.']}
    save(root / 'run-plan.json', plan)
    frozen = [root / 'run.yaml', root / 'run-plan.json', root / 'snapshot-request.json'] + sorted(p for p in scripts.iterdir() if p.is_file())
    (root / 'frozen.sha256').write_text(''.join(f'{digest(p)}  {p.relative_to(root)}\n' for p in frozen))
    print(root / 'run.yaml')


if __name__ == '__main__':
    main()
