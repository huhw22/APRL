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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--smoke', type=Path, required=True)
    args = parser.parse_args()
    w, root, smoke = args.workspace.resolve(), args.root.resolve(), args.smoke.resolve()
    if root.exists() or not (smoke / 'PASS').is_file():
        raise ValueError('Use a fresh production directory after small held field-only acceptance')
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
