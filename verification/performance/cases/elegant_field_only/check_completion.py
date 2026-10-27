#!/usr/bin/env python3
"""Check a completed held field-only run without requiring particle-plane output."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import h5py
import numpy as np
import yaml

C = 299792458.0


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    sha = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def finite(value, name):
    value = np.asarray(value)
    if value.dtype.names:
        for key in value.dtype.names:
            finite(value[key], name + '/' + key)
    elif value.dtype.kind in 'fc':
        require(np.isfinite(value).all(), 'Nonfinite data in ' + name)


def string(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def check_identity(group, manifest):
    for key in ('run_id', 'source_revision', 'configuration_digest_fnv1a64'):
        require(string(group.attrs.get(key, '')) == str(manifest[key]),
                'Diagnostic identity mismatch: ' + key)


def check_time_axis(times, cadence):
    require(len(times) > 1 and np.isfinite(times).all() and np.all(np.diff(times) > 0),
            'Invalid field time axis')
    tolerance = 8 * abs(np.spacing(float(np.max(np.abs(times))))) + 64 * np.finfo(float).eps * abs(cadence)
    require(np.all(np.abs(np.diff(times) - cadence) <= tolerance),
            'Field cadence differs beyond timestamp rounding')
    return tolerance


def checkpoint_summary(path, ranks, step, dt):
    with h5py.File(path) as file:
        header = file['header'][:]
        require(len(header) == 12 and int(header[0]) == 4706352098102169649 and int(header[1]) == 1 and int(header[2]) == 1,
                'Snapshot is not a completed schema-v1 checkpoint')
        require(int(header[3]) == ranks, 'Checkpoint MPI layout differs')
        blocks = []
        def visit(name, obj):
            if isinstance(obj, h5py.Group) and 'rank_counts' in obj and 'data' in obj:
                counts = obj['rank_counts'][:]
                require(counts.shape == (ranks,), 'Checkpoint rank-count shape: ' + name)
                require(sum(map(int, counts)) == len(obj['data']), 'Checkpoint payload length: ' + name)
                blocks.append(name)
        file['blocks'].visititems(visit)
        require(all(name in blocks for name in ('fields/ex', 'fields/ey', 'fields/ez',
                                               'fields/bx', 'fields/by', 'fields/bz',
                                               'particles', 'state/scalars', 'metadata/diagnostic_prefixes')),
                'Snapshot missing full physical-state blocks')
        prefixes = file['blocks/metadata/diagnostic_prefixes']
        raw, counts, cursor, output_names = prefixes['data'][:], prefixes['rank_counts'][:], 0, []
        for count in counts:
            count = int(count)
            item = yaml.safe_load(raw[cursor:cursor+count].tobytes().decode())
            cursor += count
            require(int(item['step']) == step, 'Checkpoint prefix step differs')
            require(np.isclose(float(item['time_box_s']), step * dt, rtol=2e-12, atol=0),
                    'Checkpoint prefix time differs from requested step')
            output_names.extend(str(output['file']) for output in item['outputs'])
        require(cursor == len(raw) and len(set(output_names)) == len(output_names),
                'Invalid snapshot diagnostic prefixes')
    return dict(file=str(path), size_bytes=path.stat().st_size, complete=True,
                step=step, ranks=ranks, physical_blocks=blocks,
                diagnostic_prefixes=output_names,
                validation_scope='Committed header and all block extents; no full state reread or restart in this production audit')


def check(root):
    report = dict(format='aprl-held-field-only-completion-v1', pass_=False,
                  output_transactions_complete=False, pure_radiation_validated=False,
                  quiet_tail_validated=False, particle_plane_output_required=False)
    try:
        card = yaml.safe_load((root / 'run.yaml').read_text())
        plan = json.loads((root / 'run-plan.json').read_text())
        variant = plan.get('case_variant', 'undulator')
        require(variant in ('undulator', 'free-drift'), 'Unknown field-only case variant')
        report['case_variant'] = variant
        if variant == 'free-drift':
            require(card['sources'].get('magnetic_elements') == [] and
                    card['sources'].get('incident_waves') == [], 'Free drift must have no prescribed sources')
            require(plan['field_plane_z_m'] == .5, 'Free-drift plane must be at 0.5 m')
            require(not card.get('checkpoint') and not plan.get('checkpoint_steps', []),
                    'Unexpected free-drift checkpoint configuration')
        geometry = json.loads(Path(plan['geometry_report']).read_text())
        require(geometry['pass'], 'Full-input geometry did not pass')
        manifest = yaml.safe_load((root / 'run-manifest.yaml').read_text())
        require(card['units'] == {'length': 'm', 'time': 's'}, 'Expected production SI card')
        require(manifest['source_revision'] == plan['compiled_revision'], 'Manifest binary revision differs')
        require(digest(plan['binary']) == plan['binary_sha256'], 'Frozen binary changed')
        require(manifest['mpi_size'] == plan['ranks'], 'Manifest MPI layout differs')
        require(yaml.safe_load(manifest['configuration_yaml']) == card, 'Manifest card differs')
        require(not card.get('particle_retirement', {}).get('enabled', False), 'Unexpected particle retirement')
        require(not card['detectors'].get('particle_planes', []), 'Unexpected particle-plane configuration')
        require(not card.get('trajectory', {}).get('enabled', False), 'Unexpected trajectory output')
        require(card['stop']['mode'] == 'after-last-element', 'Expected dynamic independent physical stop')
        held = manifest['resolved_beam_reference']
        require(held['startup_mode'] == 'held' and held['frozen_lab_plane_injection'] and
                held['hold_transverse_until_activation'], 'Held startup not preserved')
        watchdog = manifest['resolved_runtime']['shutdown_watchdog']
        require(watchdog['enabled'] and watchdog['timeout_seconds'] == 600 and
                watchdog['finalize_timeout_seconds'] == 60, 'MPI shutdown protection changed')
        mesh = manifest['resolved_mesh']
        cells = [*card['mesh']['cells'][:2], geometry['estimated_cells_z']]
        require(mesh['cells'] == cells and mesh['cell_size_m'] == card['mesh']['cell_size'],
                'Resolved grid differs from full-input geometry')
        require(mesh['center_m'] == card['mesh']['center'] and not mesh['ballistic_outside_transverse'],
                'Unexpected mesh origin or transverse ballistic policy')
        fields = card['detectors']['field_planes']
        require(len(fields) == 1 and fields[0]['name'] == plan['field_plane_name'] and
                fields[0]['z'] == plan['field_plane_z_m'], 'Field plane differs from frozen plan')
        require(not fields[0].get('particle_background', {}).get('enabled', False), 'Unexpected field subtraction')
        files = sorted((root / 'detectors').glob('*.h5'))
        path = root / 'detectors' / (fields[0]['name'] + '.h5')
        require(files == [path], 'Unexpected particle/reference/other detector file')
        log = (root / 'aprl.log').read_text(errors='replace')
        require('Simulation failure' not in log and 'mesh.duration was exhausted' not in log,
                'Run failed or exhausted the hard duration')
        stop = re.findall(r'Configured stop reached: all (\d+) valid particles passed the last interaction boundary at z=([^ ]+) m', log)
        require(len(stop) == 1 and int(stop[0][0]) > 0, 'Missing normal after-last-element stop')
        # The log prints six significant digits; exact target is frozen in the card.
        require(np.isclose(float(stop[0][1]), fields[0]['z'], rtol=2e-6, atol=0), 'Logged physical stop target differs')
        final = re.findall(r'\[resource\] phase=final step=(\d+)\b', log)
        require(len(final) == 1, 'Missing unique complete-step report')
        final_step = int(final[0]); dt = card['mesh']['cell_size'][2] / C
        require('phase=solve-complete' in log and 'phase=mpi-finalize-complete' in log,
                'Normal solve/finalization not confirmed')
        with h5py.File(path) as file:
            group = file['field_plane']; check_identity(group, manifest)
            require('adjacent MPI half-cell planes included' in string(group.attrs.get('spatial_sampling', '')),
                    'Corrected MPI spatial sampling is not identified')
            require(string(group.attrs.get('frame', '')) == 'laboratory', 'Field frame is not lab')
            require(float(group.attrs['plane_z_m']) == fields[0]['z'], 'Field HDF plane differs')
            n = int(group['committed_samples'][()]); data = group['fields']; times = group['time_s'][:]
            require(int(group['complete'][()]) == 1 and n == final_step + 1 and
                    data.shape == (n, cells[1], cells[0]) and times.shape == (n,),
                    'Field sample coverage/committed extent incomplete')
            cadence = dt / card['mesh']['boost_gamma']
            require(np.isclose(fields[0]['rhythm'], cadence, rtol=1e-12), 'Expected per-step field cadence')
            tolerance = check_time_axis(times, cadence)
            require(abs(float(times[0]) - geometry['detector_window']['first_lab_time_s']) <= tolerance,
                    'Initial detector timestamp differs from held-frame geometry')
            require(abs((times[-1] - times[0]) - final_step * cadence) <=
                    tolerance + 128 * np.finfo(float).eps * final_step * cadence,
                    'Field time window does not cover actual completed steps')
            first_box_z = geometry['detector_window']['first_box_z_m']
            last_box_z = first_box_z - geometry['beta_boost'] * C * final_step * dt
            lower, upper = geometry['physical_bounds_z_m']
            spatial_tolerance = 128 * np.finfo(float).eps * max(1., abs(lower), abs(upper))
            require(lower - spatial_tolerance <= min(first_box_z, last_box_z) and
                    max(first_box_z, last_box_z) <= upper + spatial_tolerance,
                    'Detector left the physical mesh during actual run')
            # Read bounded chunks; check every value without holding the full field history.
            chunk = max(1, (32 * 1024 * 1024) // (cells[0] * cells[1] * data.dtype.itemsize))
            for begin in range(0, n, chunk):
                finite(data[begin:min(n, begin+chunk)], 'fields')
            report['field'] = dict(file=str(path), complete=True, samples=n,
                first_lab_time_s=float(times[0]), last_lab_time_s=float(times[-1]),
                cadence_s=cadence, finite_check_scope='every scalar in every frame',
                spatial_sampling=string(group.attrs['spatial_sampling']))
        ledger_path = root / 'ledger' / card['energy_ledger']['filename']
        with h5py.File(ledger_path) as file:
            group = file['energy_ledger']; check_identity(group, manifest)
            records = group['records'][:]; count = int(group['committed_records'][()])
            require(int(group['complete'][()]) == 1 and count == len(records) and count > 0,
                    'Energy ledger incomplete')
            finite(records, 'ledger')
            if variant == 'free-drift':
                require(np.all(records['prescribed_work_J'] == 0), 'Nonzero prescribed-field work in free drift')
            require(int(records[-1]['step']) == final_step and
                    np.all(np.diff(records['step'].astype(np.int64)) > 0), 'Ledger step coverage differs')
            require(np.isclose(float(records[-1]['time_box_s']), final_step * dt, rtol=2e-12, atol=0),
                    'Final field/ledger time mismatch')
            report['ledger'] = dict(file=str(ledger_path), complete=True, records=count, final_step=final_step)
        checkpoint_dir = Path(plan['checkpoint_directory'])
        expected_paths = [checkpoint_dir / f'step-{step:012d}.h5' for step in plan['checkpoint_steps']]
        require(sorted(checkpoint_dir.glob('*.h5')) == sorted(expected_paths), 'Unexpected or missing checkpoint file')
        require(card.get('checkpoint', {}).get('steps', []) == plan['checkpoint_steps'], 'Checkpoint schedule changed')
        require(not list(checkpoint_dir.glob('*.partial*')), 'Uncommitted checkpoint artifact remains')
        report['checkpoints'] = []
        for step, checkpoint_path in zip(plan['checkpoint_steps'], expected_paths):
            require(step < final_step, 'Checkpoint was not taken before final state')
            report['checkpoints'].append(checkpoint_summary(checkpoint_path, plan['ranks'], step, dt))
        report.update(pass_=True, output_transactions_complete=True,
                      final_step=final_step, physical_stop=stop[0], cells=cells,
                      particle_boundary_reports=[line for line in log.splitlines() if line.startswith(
                          ('Frozen lab-plane injection:', 'CPML physical-domain exits:', 'Direct outer-domain exits:'))],
                      limitations=['Dynamic stop follows surviving particle state; no particle diagnostic file is required.',
                                   'This validates configured stop and recorded field coverage, not quiet tail, pure radiation or grid convergence.'])
        if plan['checkpoint_steps']:
            report['limitations'].append('Scheduled snapshot time was estimated from a previous held run; actual leading-particle position needs separate state analysis.')
        if variant == 'free-drift':
            report['limitations'].append('No undulator or incident wave is prescribed; held startup and self-consistent motion mean this is a background comparison, not a validated zero-radiation solution.')
    except (ValueError, KeyError, OSError, TypeError, IndexError) as error:
        report['error'] = str(error) or type(error).__name__
    (root / 'completion-check.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    result = check(args.root.resolve())
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['pass_'] else 1)


if __name__ == '__main__':
    main()
