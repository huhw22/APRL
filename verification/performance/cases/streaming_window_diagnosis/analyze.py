#!/usr/bin/env python3
"""Independently audit native snapshots and quantify moving-window differences.

No physical-error tolerance is imposed. Spatial differences locate observable
onsets, but do not by themselves identify a unique mechanism or prove stability.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct

import numpy as np

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
PRIOR = WORKSPACE / 'results/streaming_cpml_residual_20261027/main-r1'
NAMES = ('Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz')
E_KEYS = ('Ex_V_m', 'Ey_V_m', 'Ez_V_m')
B_KEYS = ('Bx_T', 'By_T', 'Bz_T')
CSV_NAMES = ('initial-particles.csv', 'release-particles.csv',
             'detector-crossings.csv', 'history.csv', 'field-plane.csv')
Z_STAGGER = (0., 0., .5, .5, .5, 0.)
C = np.longdouble('299792458')
ME = np.longdouble('9.1093837139e-31')
QE = np.longdouble('1.602176634e-19')
BOOST = np.longdouble(5) / 3
BETA = np.longdouble('.8')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    out = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            out.update(block)
    return out.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def table(path):
    values = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
    require(values.size > 0 and values.dtype.names, 'Empty CSV: ' + str(path))
    require(all(np.isfinite(values[key]).all() for key in values.dtype.names),
            'Nonfinite CSV: ' + str(path))
    return values


def expected_offset(step, guard):
    return math.floor(218. - .8 * step) + 84 - guard - 20


def read_exact(stream, size):
    result = stream.read(size)
    require(len(result) == size, 'Truncated field-section binary')
    return result


def sections(directory, metrics):
    manifest = json.loads((directory / 'field-sections.json').read_text())
    require(manifest['format'] == 'aprl-native-field-xz-sections-v1' and
            manifest['completed'] is True and manifest['byte_order'] == 'little',
            'Invalid native-section manifest')
    maximum = int(manifest['maximum_state_step'])
    require(0 < maximum <= metrics['steps'] and
            manifest['records'] == 1 + 3 * maximum and manifest['section_j'] == 14,
            'Invalid section record count or scope')
    expected = [(0, 0)] + [(step, stage) for step in range(1, maximum + 1)
                           for stage in (1, 2, 3)]
    nx, ny, nz = metrics['nx'], metrics['ny'], metrics['window_nz']
    dimensions = ((nx, ny+1, nz+1), (nx+1, ny, nz+1), (nx+1, ny+1, nz),
                  (nx+1, ny, nz), (nx, ny+1, nz), (nx, ny, nz+1))
    result = {}
    with (directory / 'field-sections.bin').open('rb') as stream:
        require(read_exact(stream, 8) == b'APRLXZ01', 'Incorrect section magic')
        for expected_step, expected_stage in expected:
            step, stage, offset, count = struct.unpack('<4Q', read_exact(stream, 32))
            require((step, stage) == (expected_step, expected_stage) and count == 6,
                    'Incorrect native stage/order or component count')
            lattice_step = step-1 if stage in (1, 2) else step
            wanted_offset = (0 if metrics['mode'] == 'full' else
                             expected_offset(lattice_step, metrics['guard_cells']))
            require(offset == wanted_offset, 'Wrong world offset at native stage')
            fields = []
            for component, shape in enumerate(dimensions):
                a, cx, cy, cz, j = struct.unpack('<5Q', read_exact(stream, 40))
                require(a == component and (cx, cy, cz) == shape and j == 14,
                        'Section dimensions/component order changed')
                field = np.frombuffer(read_exact(stream, 8*cx*cz), dtype='<f8').reshape(cz, cx)
                require(np.isfinite(field).all(), 'Nonfinite native field section')
                fields.append(field)
            e_step = step-1 if stage == 1 else step
            b_step = step-.5
            result[(step, stage)] = {'offset': offset, 'fields': fields,
                                     'E_step': e_step, 'B_step': b_step}
        require(not stream.read(1), 'Trailing bytes after expected section records')
    return result, manifest


def difference(reference, candidate):
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    require(reference.shape == candidate.shape and reference.size > 0, 'Invalid difference inputs')
    delta = candidate-reference
    source = float(np.linalg.norm(reference.ravel()))
    error = float(np.linalg.norm(delta.ravel()))
    return {'relative_l2': error/source if source else None,
            'absolute_rms': float(np.sqrt(np.mean(delta*delta))),
            'maximum_absolute': float(np.max(np.abs(delta))),
            'reference_l2': source, 'difference_l2': error,
            'values': int(delta.size)}


def particles(path):
    values = np.sort(table(path), order='id')
    require(len(np.unique(values['id'])) == len(values), 'Repeated particle ID')
    return values


def particle_moments(events, initial):
    require(np.array_equal(events['id'], initial['id']), 'Different event and initial cohorts')
    weights = np.abs(events['charge_C']).astype(np.longdouble)
    require(weights.sum() > 0, 'Zero charge weight')
    weights /= weights.sum()
    u = np.column_stack([events[k] for k in ('ux', 'uy', 'uz')]).astype(np.longdouble)
    u0 = np.column_stack([initial[k] for k in ('ux', 'uy', 'uz')]).astype(np.longdouble)
    gamma = np.sqrt(1 + (u*u).sum(axis=1))
    gamma0 = np.sqrt(1 + (u0*u0).sum(axis=1))
    delta_gamma = ((u-u0)*(u+u0)).sum(axis=1)/(gamma+gamma0)
    rest = ME*C*C/QE
    delta_box = delta_gamma*rest
    delta_lab = BOOST*(delta_gamma+BETA*(u[:, 2]-u0[:, 2]))*rest
    require(np.max(np.abs(delta_box-events['delta_energy_box_eV'])) < 1e-9 and
            np.max(np.abs(delta_lab-events['delta_energy_lab_eV'])) < 1e-9,
            'Event energy columns disagree with independently computed momentum change')
    lab_energy = BOOST*(gamma+BETA*u[:, 2])*rest
    mean_energy = np.sum(weights*lab_energy)
    spread = np.sqrt(np.sum(weights*(lab_energy-mean_energy)**2))
    mean_delta = np.sum(weights*delta_lab)
    return {'mean_delta_energy_lab_eV': float(mean_delta),
            'mean_energy_loss_lab_eV': -float(mean_delta),
            'sigma_energy_lab_eV': float(spread),
            'count': len(events)}


def plane_review(full, candidate, full_events, events):
    keys = ('step', 'pixel_i', 'pixel_j', 'time_E_box_s', 'time_B_box_s',
            'time_E_lab_s', 'x_m', 'y_m', 'z_box_m')
    for key in keys:
        require(np.array_equal(full[key], candidate[key]), 'Plane sample differs: ' + key)
    vectors = lambda a, names: np.column_stack([a[key] for key in names])
    first_detector = min(float(np.min(full_events['step']+full_events['fraction'])),
                         float(np.min(events['step']+events['fraction'])))
    last_detector = min(float(np.max(full_events['step']+full_events['fraction'])),
                        float(np.max(events['step']+events['fraction'])))
    intervals = {'whole': np.ones(len(full), dtype=bool),
                 'before_first_detector': full['step'] < first_detector,
                 'until_common_last_detector': full['step'] <= last_detector}
    result = {'intervals': {}, 'by_step': [], 'first_detector_continuous_step': first_detector,
              'common_last_detector_continuous_step': last_detector}
    for name, mask in intervals.items():
        result['intervals'][name] = {
            label: difference(vectors(full[mask], names), vectors(candidate[mask], names))
            for label, names in (('E', E_KEYS), ('B', B_KEYS))}
    for step in np.unique(full['step']):
        mask = full['step'] == step
        row = {'step': int(step)}
        for label, names in (('E', E_KEYS), ('B', B_KEYS)):
            row[label] = difference(vectors(full[mask], names), vectors(candidate[mask], names))
        delta = vectors(candidate[mask], B_KEYS)-vectors(full[mask], B_KEYS)
        row['cDeltaB_vector_rms_V_m'] = float(C*np.sqrt(np.mean(np.sum(delta*delta, axis=1))))
        result['by_step'].append(row)
    return result


def inspect_case(directory):
    metrics = json.loads((directory / 'metrics.json').read_text())
    require(metrics['completed'] is True and metrics['mode'] in
            ('full', 'residual', 'residual-refresh'), 'Incomplete or unexpected case')
    require(metrics['beta_box'] == .23 and metrics['guard_cells'] == 24 and
            metrics['steps'] == 176 and metrics['cpml_cells'] == [0, 0, 20],
            'This diagnostic campaign must use the fixed moving g24 case')
    native, manifest = sections(directory, metrics)
    initial = particles(directory / 'initial-particles.csv')
    events = particles(directory / 'detector-crossings.csv')
    release = particles(directory / 'release-particles.csv')
    require(len(initial) == len(events) == len(release) == 1024, 'Incomplete particle cohort')
    plane = table(directory / 'field-plane.csv')
    require(len(plane) == metrics['field_plane_samples'] == 45*400, 'Incomplete field-plane CSV')
    require(np.allclose(plane['time_E_box_s'], plane['step']*metrics['dt_s'], rtol=1e-14, atol=1e-27)
            and np.allclose(plane['time_B_box_s'], (plane['step']-.5)*metrics['dt_s'], rtol=1e-14, atol=1e-27),
            'Plane native time layers changed')
    return {'directory': directory, 'metrics': metrics, 'native': native, 'manifest': manifest,
            'initial': initial, 'events': events, 'release': release, 'plane': plane,
            'moments': particle_moments(events, initial)}


def compare_sections(full, candidate, metrics_writer, rows_writer):
    require(set(full['native']) == set(candidate['native']), 'Different diagnostic stage coverage')
    summary, onset = [], {}
    width = candidate['metrics']['window_nz']
    name = candidate['directory'].name
    heatmaps = {a: [] for a in (3, 4, 5)}
    for key, record in candidate['native'].items():
        step, stage = key
        reference = full['native'][key]
        require(record['E_step'] == reference['E_step'] and record['B_step'] == reference['B_step'],
                'Native field time mismatch')
        start = record['offset']-reference['offset']
        for a, field in enumerate(record['fields']):
            source = reference['fields'][a][start:start+field.shape[0]]
            require(source.shape == field.shape and start >= 0, 'Reference does not cover candidate world cells')
            delta = field-source
            q = np.arange(len(field), dtype=float)+Z_STAGGER[a]
            masks = {'all': np.ones(len(field), dtype=bool), 'low_cpml': q < 20,
                     'interior': (q >= 20) & (q <= width-20), 'high_cpml': q > width-20}
            scale = 1. if a < 3 else float(C)
            for region, mask in masks.items():
                diff = difference(source[mask], field[mask])
                row = {'case': name, 'step': step, 'stage': stage, 'component': NAMES[a],
                       'region': region, 'offset': record['offset'],
                       'native_time_step': record['E_step'] if a < 3 else record['B_step'], **diff,
                       'E_or_cB_absolute_rms_V_m': diff['absolute_rms']*scale,
                       'E_or_cB_maximum_absolute_V_m': diff['maximum_absolute']*scale}
                metrics_writer.writerow(row)
                if region != 'all':
                    for threshold in (1e-6, 1e-3, 1.):
                        identity = f'{NAMES[a]}:{region}:{threshold:g}_V_m'
                        if identity not in onset and diff['maximum_absolute']*scale > threshold:
                            onset[identity] = row
                if region == 'interior' and a >= 3:
                    summary.append(row)
            row_rms = np.sqrt(np.mean(delta*delta, axis=1))
            for k, value in enumerate(row_rms):
                region = 'low_cpml' if q[k] < 20 else ('high_cpml' if q[k] > width-20 else 'interior')
                rows_writer.writerow({'case': name, 'step': step, 'stage': stage,
                    'component': NAMES[a], 'local_k': k, 'world_grid_z': record['offset']+q[k],
                    'world_z_m': candidate['metrics']['full_origin_z_m']+(record['offset']+q[k])*1e-6,
                    'region': region, 'absolute_row_rms': float(value),
                    'E_or_cB_row_rms_V_m': float(value*scale)})
            if stage in (0, 3) and a >= 3:
                heatmaps[a].append((step, row_rms*scale))
    return {'first_threshold_crossings': onset,
            'maximum_interior_cB_component_error_V_m': max(
                row['E_or_cB_maximum_absolute_V_m'] for row in summary)}, heatmaps


def plots(output, comparisons, cases, heatmaps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True)
    for comparison in comparisons:
        records = comparison['plane']['by_step']
        label = comparison['mode']
        axes[0].plot([v['step'] for v in records], [100*v['B']['relative_l2'] for v in records], label=label)
        axes[1].plot([v['step'] for v in records], [v['cDeltaB_vector_rms_V_m'] for v in records], label=label)
    axes[0].set_ylabel('Plane B relative L2 difference (%)')
    axes[1].set_ylabel('Plane cΔB vector RMS (V/m)')
    axes[1].set_xlabel('Completed step; native B at step - 1/2')
    for axis in axes:
        axis.grid(alpha=.25)
        axis.legend()
    fig.savefig(output / 'native_B_error.png', dpi=180)
    fig.savefig(output / 'native_B_error.pdf')
    plt.close(fig)
    fig, axes = plt.subplots(len(cases), 3, figsize=(12, 3.6*len(cases)), squeeze=False,
                             constrained_layout=True)
    for row, case in enumerate(cases):
        width = case['metrics']['window_nz']
        name = case['directory'].name
        for col, a in enumerate((3, 4, 5)):
            records = heatmaps[name][a]
            steps = np.array([v[0] for v in records])
            data = np.stack([v[1] for v in records])
            coordinates = np.arange(data.shape[1])+Z_STAGGER[a]
            maximum = max(float(data.max()), 1e-8)
            picture = axes[row, col].pcolormesh(coordinates, steps,
                np.maximum(data, maximum*1e-10), shading='nearest',
                norm=LogNorm(vmin=maximum*1e-10, vmax=maximum), cmap='magma')
            axes[row, col].axvline(20, color='cyan', linewidth=1)
            axes[row, col].axvline(width-20, color='cyan', linewidth=1)
            axes[row, col].set_title(case['metrics']['mode']+' '+NAMES[a])
            axes[row, col].set_xlabel('Native z coordinate relative to window (cells)')
            axes[row, col].set_ylabel('Completed step; after remap')
            fig.colorbar(picture, ax=axes[row, col], label='cΔB component row RMS (V/m)')
    fig.suptitle('Same-world-cell difference versus full; j=14; cyan: CPML interfaces')
    fig.savefig(output / 'native_B_time_z.png', dpi=180)
    fig.savefig(output / 'native_B_time_z.pdf')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    require(not (root / 'FAILED.json').exists(), 'Campaign has failed')
    require((root / 'COMPLETED.json').is_file(), 'Campaign is not complete')
    cases = [inspect_case(path.parent) for path in sorted(root.glob('*/field-sections.json'))]
    require(len(cases) == 3 and {v['metrics']['mode'] for v in cases} ==
            {'full', 'residual', 'residual-refresh'}, 'Expected exactly the three moving g24 cases')
    full = next(case for case in cases if case['metrics']['mode'] == 'full')
    candidates = [case for case in cases if case is not full]
    output = root / 'diagnosis'
    output.mkdir(exist_ok=True)
    report = {'format': 'aprl-streaming-window-diagnosis-analysis-v1',
              'script_sha256': digest(Path(__file__)), 'root': str(root),
              'legacy_csv_replay_status': 'PASS', 'legacy_replay': [], 'comparisons': [],
              'input_sha256': {}, 'full_particle_moments': full['moments']}
    for case in cases:
        directory = case['directory']
        for filename in CSV_NAMES+('metrics.json', 'field-sections.bin', 'field-sections.json'):
            report['input_sha256'][str(directory/filename)] = digest(directory/filename)
        if case['metrics']['mode'] in ('full', 'residual'):
            previous = PRIOR / ('b023-'+case['metrics']['mode']+'-g24')
            for filename in CSV_NAMES:
                require(digest(directory/filename) == digest(previous/filename),
                        'Instrumented legacy CSV differs: '+str(directory/filename))
            report['legacy_replay'].append({'current': str(directory), 'previous': str(previous),
                                            'five_csv_byte_identical': True})
    metric_columns = ('case', 'step', 'stage', 'component', 'region', 'offset', 'native_time_step',
                      'relative_l2', 'absolute_rms', 'maximum_absolute', 'reference_l2',
                      'difference_l2', 'values', 'E_or_cB_absolute_rms_V_m', 'E_or_cB_maximum_absolute_V_m')
    row_columns = ('case', 'step', 'stage', 'component', 'local_k', 'world_grid_z', 'world_z_m',
                   'region', 'absolute_row_rms', 'E_or_cB_row_rms_V_m')
    all_heatmaps = {}
    with (output/'section-metrics.csv').open('w', newline='') as stream, \
            (output/'section-z-errors.csv').open('w', newline='') as rows:
        writer = csv.DictWriter(stream, fieldnames=metric_columns)
        row_writer = csv.DictWriter(rows, fieldnames=row_columns)
        writer.writeheader()
        row_writer.writeheader()
        for case in candidates:
            for cohort in ('initial', 'release'):
                for key in full[cohort].dtype.names:
                    require(np.array_equal(full[cohort][key], case[cohort][key]),
                            'Prescribed input/release differs: '+cohort+':'+key)
            native, heatmap = compare_sections(full, case, writer, row_writer)
            all_heatmaps[case['directory'].name] = heatmap
            plane = plane_review(full['plane'], case['plane'], full['events'], case['events'])
            report['comparisons'].append({'path': str(case['directory']), 'mode': case['metrics']['mode'],
                'native_sections': native, 'plane': plane, 'particle_moments': case['moments'],
                'particle_moment_differences': {key: case['moments'][key]-full['moments'][key]
                                               for key in full['moments'] if key != 'count'},
                'diagnostic_maximum_step': case['manifest']['maximum_state_step']})
    report['interpretation'] = [
        'Field differences are measured, not assigned physical acceptance by this script.',
        'Each section component is compared at identical native lattice position and time; components are not collocated.',
        'Regions use native z: low <20, interior 20 through width-20 inclusive, high >width-20.',
        'Onsets use absolute E or cB thresholds 1e-6, 1e-3, and 1 V/m; these are diagnostic markers, not accuracy limits.',
        'Stage 1 includes the magnetic boundary hook; stage 2 includes the electric boundary hook; stage 3 includes remap and optional endpoint refresh.',
        'Spatial norms alone cannot identify the unique source of an error; no pure-radiation or energy-flux inference is made.',
        'Partial section duration cannot establish long-time stability; field-plane comparisons cover the complete 176-step fixture.']
    if not args.no_plots:
        plots(output, report['comparisons'], candidates, all_heatmaps)
    save(output/'REVIEW.json', report)
    print(output/'REVIEW.json')


if __name__ == '__main__':
    main()
