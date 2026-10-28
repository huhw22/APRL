#!/usr/bin/env python3
"""Inspect complete native diagnostics; report approximation separately from replay."""
import argparse
import json
from pathlib import Path
import re

import h5py
import numpy as np

REST_EV = np.longdouble('510998.9499961642')


def require(test, message):
    if not test:
        raise RuntimeError(message)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def finite(data):
    if data.dtype.names:
        return all(finite(data[name]) for name in data.dtype.names)
    return bool(np.isfinite(data).all()) if data.dtype.kind in 'fciu' else True


def read(path):
    with h5py.File(path) as stream:
        require(len(stream) == 1, 'Unexpected groups: '+str(path))
        group = stream[next(iter(stream))]
        require(int(group['complete'][()]) == 1, 'Incomplete HDF5: '+str(path))
        data = {key: group[key][()] for key in group if isinstance(group[key], h5py.Dataset)}
        require(all(finite(value) for value in data.values()), 'Non-finite HDF5: '+str(path))
        for commit, name in [('committed_records', 'records'), ('committed_samples', 'time_s')]:
            if commit in data:
                require(int(data[commit]) == len(data[name]), 'Truncated dataset: '+str(path))
        if 'fields' in data:
            require(len(data['fields']) == len(data['time_s']), 'Field/time count mismatch')
        return data, dict(group.attrs)


def ordered(data, count, label):
    require(len(data) == count and len(np.unique(data['source_id'])) == count,
            'Missing/duplicate particle cohort: '+label)
    data = data[np.argsort(data['source_id'])]
    require(np.array_equal(data['source_id'], np.arange(1, count+1)), 'Wrong source IDs: '+label)
    return data


def ld(data):
    return np.asarray(data, dtype=np.longdouble)


def mean(x, weights):
    return np.sum(ld(x)*weights)/np.sum(weights)


def sigma(x, weights):
    return np.sqrt(mean((ld(x)-mean(x, weights))**2, weights))


def covariance(x, y, weights):
    return mean((ld(x)-mean(x, weights))*(ld(y)-mean(y, weights)), weights)


def energy_change(u, u0):
    u, u0 = ld(u), ld(u0)
    g, g0 = np.sqrt(1+np.sum(u*u, axis=1)), np.sqrt(1+np.sum(u0*u0, axis=1))
    return np.sum((u-u0)*(u+u0), axis=1)/(g+g0)*REST_EV


def particles(release, records):
    weights = np.abs(ld(release['charge_C']))
    delta = energy_change(records['proper_velocity'], release['proper_velocity'])
    gamma = np.sqrt(1+np.sum(ld(records['proper_velocity'])**2, axis=1))
    gamma0 = np.sqrt(1+np.sum(ld(release['proper_velocity'])**2, axis=1))
    result = {'count': len(records), 'mean_delta_K_eV': float(mean(delta, weights)),
        'rms_delta_K_eV': float(np.sqrt(mean(delta**2, weights))),
        'std_delta_K_eV': float(sigma(delta, weights)),
        'sigma_K_eV': float(sigma(gamma, weights)*REST_EV),
        'delta_sigma_K_eV': float((sigma(gamma, weights)-sigma(gamma0, weights))*REST_EV),
        'relative_sigma_gamma': float(sigma(gamma, weights)/mean(gamma, weights)), 'transverse': {}}
    ballistic = ld(release['position_m']) + (
        ld(release['proper_velocity'])/ld(release['proper_velocity'])[:, 2:3] *
        (ld(records['position_m'])[:, 2:3]-ld(release['position_m'])[:, 2:3]))
    for axis, name in enumerate(('x', 'y')):
        position, momentum = records['position_m'][:, axis], records['proper_velocity'][:, axis]
        result['transverse'][name] = {
            'sigma_m': float(sigma(position, weights)),
            'ballistic_sigma_m': float(sigma(ballistic[:, axis], weights)),
            'sigma_minus_ballistic_m': float(sigma(position, weights)-sigma(ballistic[:, axis], weights)),
            'cov_position_u_m': float(covariance(position, momentum, weights)),
            'cov_minus_ballistic_m': float(covariance(position, momentum, weights)-
                covariance(ballistic[:, axis], release['proper_velocity'][:, axis], weights))}
    return result


def array_error(reference, candidate):
    difference = ld(candidate)-ld(reference)
    norm = np.sqrt(np.sum(ld(reference)**2))
    return {'bitwise_equal': reference.dtype == candidate.dtype and reference.tobytes() == candidate.tobytes(),
        'max_abs': float(np.max(np.abs(difference))) if difference.size else 0.,
        'rms_abs': float(np.sqrt(np.mean(difference**2))) if difference.size else 0.,
        'relative_l2_reference_norm': float(np.sqrt(np.sum(difference**2))/norm) if norm else None}


def inspect_case(directory, plan, source):
    count = plan['source']['records']
    log = (directory/'aprl.log').read_text()
    require('Configured stop reached:' in log, 'Configured stop missing: '+str(directory))
    launch = json.loads((directory/'aprl.log.launch.json').read_text())
    require(launch['exit_code'] == 0, 'Nonzero launch: '+str(directory))
    files = {str(path.relative_to(directory)): read(path) for folder in ('detectors', 'trajectory', 'ledger')
             for path in sorted((directory/folder).glob('*.h5'))}
    require(len(files) == 5+plan['ranks'], 'Unexpected HDF5 set: '+str(directory))
    planes = {name: ordered(files['detectors/'+name+'.h5'][0]['records'], count, name)
              for name in ('release', 'aux', 'exit')}
    for name, data in planes.items():
        require(np.all(data['crossing_direction'] == 1), 'Unexpected backward crossing: '+name)
        for key in ('particle_id', 'source_id', 'charge_C', 'mass_kg', 'weight'):
            require(np.array_equal(planes['release'][key], data[key]), 'Particle identity/weight changed: '+key)
    trajectory = np.concatenate([data[0]['records'] for name, data in files.items() if name.startswith('trajectory/')])
    samples = trajectory[trajectory['event_type'] == 0]
    samples = samples[np.lexsort((samples['time_s'], samples['source_id']))]
    _, first = np.unique(samples['source_id'], return_index=True)
    first_recorded = ordered(samples[first], count, 'first recorded released trajectory')
    result = {'hdf5_files_complete': len(files), 'elapsed_s_including_initialization_and_IO': launch['seconds'],
        'source_release_u_error': array_error(source['proper_velocity'], planes['release']['proper_velocity']),
        'source_release_xy_error_m': array_error(source['plane_position_m'], planes['release']['position_m'][:, :2]),
        'particle_planes': {name: particles(planes['release'], data) for name, data in planes.items()}}
    # Every ID must cross the auxiliary plane. Its entire preceding full step
    # stays before the fringe even at the maximum relativistic lab displacement.
    dz, boost = 1.2543262415359034e-6, 929.1305492858543
    beta = np.sqrt(1-boost**-2)
    maximum_step_lab_dz = boost*(1+beta)*dz
    aux_max_z = float(np.max(planes['aux']['position_m'][:, 2]))
    clearance = plan['geometry']['fringe_start_lab_z_m']-aux_max_z-maximum_step_lab_dz
    require(clearance > 0, 'Auxiliary-plane fringe clearance proof failed')
    result['auxiliary_ballistic_scope'] = {'all_ids_crossed': True,
        'worst_full_step_lab_dz_m': float(maximum_step_lab_dz), 'remaining_fringe_clearance_m': clearance,
        'interpretation': 'Aux response excludes direct magnetic fringe. Exit response also includes the undulator.'}
    energy = files['ledger/energy-ledger.h5'][0]['records']
    result['ledger_final'] = {key: float(energy[key][-1]) for key in (
        'closure_valid', 'field_energy_interior_J', 'outward_field_energy_J',
        'particle_kinetic_active_J', 'particle_kinetic_removed_J', 'balance_residual_J')}
    if directory.name == 'moving':
        match = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)', log)
        require(match is not None, 'Missing moving-window completion report')
        result['moving_exit'] = dict(zip(('steps', 'particle_exits', 'cropped_carriers'), map(int, match.groups())))
        require(result['moving_exit']['particle_exits'] > 0 and result['moving_exit']['cropped_carriers'] > 0,
            'Moving right particle exit/carrier crop was not exercised')
        ledgers = sorted((directory/'reference').glob('window-rank-*.csv'))
        require(len(ledgers) == plan['ranks'], 'Missing moving-volume rank ledger')
        last = [np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))[-1] for path in ledgers]
        result['moving_volume_final_totals'] = {name: float(sum(row[name] for row in last)) for name in (
            'inserted_field_J', 'cropped_field_J', 'particle_exits', 'exit_charge_C', 'cropped_carriers', 'cropped_charge_C')}
        geometry_match = re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+).*fixed_lab_particle_exit_m=([^ ]+)', log)
        require(geometry_match is not None, 'Missing actual moving geometry')
        result['actual_moving_geometry'] = dict(zip(('reference_cells', 'active_cells', 'initial_offset'),
            map(int, geometry_match.groups()[:3])))
        result['actual_moving_geometry']['safe_exit_lab_z_m'] = float(geometry_match[4])
        geo = plan['geometry']
        if geo.get('geometry_contract_version', 1) >= 2:
            actual = result['actual_moving_geometry']
            require(actual['reference_cells'] == geo['reference_cells'] and
                actual['active_cells'] == geo['active_cells_estimate'] and
                actual['initial_offset'] == geo['initial_offset_estimate'], 'Runtime geometry differs from prepared contract')
            require(abs(actual['safe_exit_lab_z_m']-geo['safe_exit_lab_z_estimate_m']) <= 2e-10,
                'Runtime safe exit differs from held-centred geometry prediction')
            require(geo['witness_relative_to_initial_inlet_cells']['2'] < 0 and
                geo['witness_entry_estimates']['2']['entry_type'] == 'physical_drift',
                'Prepared source ID2 no longer exercises the left physical entry')
            result['actual_moving_geometry']['prepared_contract_verified'] = True
        else:
            result['actual_moving_geometry']['prepared_contract_verified'] = False
            result['actual_moving_geometry']['note'] = 'Version1 estimates omitted held boxTranslation; use actual log values.'
        result['ledger_scope'] = 'Stationary ledger residual alone is not closure for a moving volume; native insert/crop accounting is reported separately.'
    return {'summary': result, 'files': files, 'planes': planes, 'first_recorded': first_recorded}


def field_regions(reference, candidate, prefix):
    """Keep E and B in separate SI units; never concatenate unscaled E/B."""
    first, attrs = reference['files']['detectors/exit-field.h5']
    second = candidate['files']['detectors/exit-field.h5'][0]
    times = first['time_s']
    nt, ny, nx = first['fields'].shape
    x = float(attrs['x_first_m'])+np.arange(nx)*float(attrs['dx_m'])
    y = float(attrs['y_first_m'])+np.arange(ny)*float(attrs['dy_m'])
    xx, yy = np.meshgrid(x, y)
    area = float(attrs['dx_m'])*float(attrs['dy_m'])
    release = reference['planes']['release']
    weights = np.abs(ld(release['charge_C']))
    centre = np.array([float(mean(release['position_m'][:, axis], weights)) for axis in (0, 1)])
    width = np.array([float(sigma(release['position_m'][:, axis], weights)) for axis in (0, 1)])
    masks = {'all_aperture': np.ones((ny, nx), dtype=bool),
        'physical_xy': (np.abs(xx)<350e-6) & (np.abs(yy)<350e-6),
        'beamcore_2sigma': ((xx-centre[0])/width[0])**2+((yy-centre[1])/width[1])**2<=4}
    periods = {'all_output': np.ones(nt, dtype=bool), 'through_complete_detector_cohort': prefix,
               'after_complete_detector_cohort': ~prefix}
    result = {'normalization': 'E and B separately. Local/plane L2 uses fixed vector norm; no mixed SI norm. Area sums of amplitude squared are not Poynting power.',
        'core_definition': {'centroid_m': centre.tolist(), 'sigma_m': width.tolist(), 'ellipse_radius_sigma': 2},
        'physical_xy_definition': '|x|<350um and |y|<350um excludes the original 10-cell transverse CPML', 'fields': {}}
    for field, labels in [('electric_V_per_m', ('Ex', 'Ey', 'Ez')), ('magnetic_T', ('Bx', 'By', 'Bz'))]:
        a, b = first['fields'][field], second['fields'][field]
        output = {}
        for name, mask in masks.items():
            points = np.argwhere(mask)
            aa, bb = a[:, mask, :], b[:, mask, :]
            difference = bb-aa
            a2, b2, d2 = (np.sum(value*value, axis=(1, 2)) for value in (aa, bb, difference))
            count = len(points)
            regions = {'pixel_count': count, 'area_m2': count*area, 'windows': {}}
            for period, selected in periods.items():
                ids = np.flatnonzero(selected)
                if not len(ids):
                    continue
                av, bv, delta = aa[selected], bb[selected], difference[selected]
                index = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
                ti, pi, component = int(ids[index[0]]), index[1], index[2]
                j, i = points[pi]
                local_norm = np.linalg.norm(aa[ti, pi])
                norm = np.sum(a2[selected])
                nonzero = a2[selected]>max(a2.max()*1e-12, np.finfo(float).tiny)
                relative = np.sqrt(d2[selected][nonzero]/a2[selected][nonzero])
                peak_a, peak_b = int(ids[np.argmax(a2[selected])]), int(ids[np.argmax(b2[selected])])
                regions['windows'][period] = {'sample_count': len(ids), 'sample_index_range': [int(ids[0]), int(ids[-1])],
                    'relative_L2_difference_over_fixed': float(np.sqrt(np.sum(d2[selected])/norm)),
                    'fixed_vector_rms': float(np.sqrt(norm/(len(ids)*count))),
                    'moving_vector_rms': float(np.sqrt(np.sum(b2[selected])/(len(ids)*count))),
                    'relative_trace_sum_area_amplitude_squared_change': float((np.sum(b2[selected])-norm)/norm),
                    'fixed_trace_sum_area_amplitude_squared': float(norm*area),
                    'moving_trace_sum_area_amplitude_squared': float(np.sum(b2[selected])*area),
                    'fixed_peak_sample': peak_a, 'moving_peak_sample': peak_b,
                    'peak_time_shift_s': float(times[peak_b]-times[peak_a]),
                    'max_instantaneous_relative_L2': float(relative.max()) if len(relative) else None,
                    'max_instantaneous_relative_L2_sample': int(ids[nonzero][np.argmax(relative)]) if len(relative) else None,
                    'instantaneous_denominator_floor_relative_to_peak_squared_norm': 1e-12,
                    'peak_plane_error_sample': int(ids[np.argmax(d2[selected])]),
                    'component_relative_L2': {labels[k]: float(np.sqrt(np.sum(delta[:, :, k]**2)/np.sum(av[:, :, k]**2)))
                        if np.any(av[:, :, k]) else None for k in range(3)},
                    'max_absolute_component_difference': {'component': labels[component], 'value': float(abs(delta[index])),
                        'fixed_value': float(aa[ti, pi, component]), 'moving_value': float(bb[ti, pi, component]),
                        'sample_index': ti, 'lab_time_s': float(times[ti]), 'time_relative_first_fs': float((times[ti]-times[0])*1e15),
                        'x_m': float(xx[j, i]), 'y_m': float(yy[j, i]),
                        'local_vector_relative_difference': float(np.linalg.norm(difference[ti, pi])/local_norm) if local_norm else None}}
            regions['time_traces'] = {'relative_lab_time_fs': ((times-times[0])*1e15).tolist(),
                'fixed_vector_rms': np.sqrt(a2/count).tolist(), 'moving_vector_rms': np.sqrt(b2/count).tolist(),
                'difference_vector_rms': np.sqrt(d2/count).tolist()}
            output[name] = regions
        # Compare the same previously identified site AND the new global maxima.
        # This prevents a dense run from appearing improved merely because the
        # worst pixel or time moved to a different location.
        i, j = int(np.argmin(abs(x+37.5e-6))), int(np.argmin(abs(y+22.5e-6)))
        if nt>89:
            av, bv = a[89, j, i], b[89, j, i]
            output['historical_sparse_max_site'] = {'sample_index': 89, 'x_m': float(x[i]), 'y_m': float(y[j]),
                'fixed_vector': av.tolist(), 'moving_vector': bv.tolist(), 'difference_vector': (bv-av).tolist(),
                'relative_vector_difference': float(np.linalg.norm(bv-av)/np.linalg.norm(av)) if np.linalg.norm(av) else None,
                'scope': 'Same site as single-r1 sparse maximum; always also inspect the present global maximum.'}
        result['fields'][field] = output
    return result


def particle_response_differences(reference, candidate, comparison):
    result = {}
    for name in ('aux', 'exit'):
        a, b = reference['summary']['particle_planes'][name], candidate['summary']['particle_planes'][name]
        item = {'fixed_mean_delta_K_eV': a['mean_delta_K_eV'], 'moving_mean_delta_K_eV': b['mean_delta_K_eV'],
            'difference_mean_delta_K_eV': b['mean_delta_K_eV']-a['mean_delta_K_eV'],
            'per_id_K_difference_rms_eV': comparison['planes'][name]['rms_candidate_minus_reference_K_eV'],
            'per_id_difference_rms_over_fixed_self_delta_rms': comparison['planes'][name]['rms_candidate_minus_reference_K_eV']/a['rms_delta_K_eV'] if a['rms_delta_K_eV'] else None,
            'sigma_K_difference_eV': b['sigma_K_eV']-a['sigma_K_eV'],
            'fixed_delta_sigma_K_eV': a['delta_sigma_K_eV'], 'moving_delta_sigma_K_eV': b['delta_sigma_K_eV'],
            'transverse_ballistic_corrected': {}}
        for axis in ('x', 'y'):
            item['transverse_ballistic_corrected'][axis] = {}
            for key in ('sigma_minus_ballistic_m', 'cov_minus_ballistic_m'):
                av, bv = a['transverse'][axis][key], b['transverse'][axis][key]
                item['transverse_ballistic_corrected'][axis][key] = {'fixed': av, 'moving': bv, 'difference': bv-av,
                    'relative_difference_over_fixed_response': (bv-av)/av if av else None}
        result[name] = item
    result['scope'] = 'Compare approximation error with the measured self-field response, not total beam energy/width. Exit includes undulator dynamics; auxiliary plane excludes fringe. No automatic physical tolerance is invented.'
    return result


def compare(reference, candidate, label):
    result = {'comparison': label, 'first_recorded_state': {}, 'planes': {},
        'first_recorded_state_scope': 'Formal trajectory omits held records. These are the first released samples, not the t_box=0 state; release invariance is checked separately.'}
    for key in ('source_id', 'particle_id', 'position_m', 'proper_velocity', 'charge_C', 'weight', 'time_s'):
        result['first_recorded_state'][key] = array_error(reference['first_recorded'][key], candidate['first_recorded'][key])
    for name in ('release', 'aux', 'exit'):
        a, b = reference['planes'][name], candidate['planes'][name]
        weights = np.abs(ld(a['charge_C']))
        for key in ('particle_id', 'source_id', 'charge_C', 'mass_kg', 'weight'):
            require(np.array_equal(a[key], b[key]), label+': identity/weight changed')
        delta = energy_change(b['proper_velocity'], a['proper_velocity'])
        result['planes'][name] = {'mean_candidate_minus_reference_K_eV': float(mean(delta, weights)),
            'rms_candidate_minus_reference_K_eV': float(np.sqrt(mean(delta**2, weights))),
            'position_m': array_error(a['position_m'], b['position_m']),
            'proper_velocity': array_error(a['proper_velocity'], b['proper_velocity']),
            'crossing_time_s': array_error(a['time_s'], b['time_s'])}
    a, attrs_a = reference['files']['detectors/exit-field.h5']
    b, attrs_b = candidate['files']['detectors/exit-field.h5']
    require(a['fields'].shape == b['fields'].shape, label+': field shapes differ')
    require(np.allclose(a['time_s'], b['time_s'], rtol=64*np.finfo(float).eps, atol=0),
            label+': lab sample clocks differ; no interpolation is allowed')
    for attr in ('plane_z_m', 'x_first_m', 'y_first_m', 'dx_m', 'dy_m', 'magnetic_time_stagger'):
        require(attr in attrs_a and attr in attrs_b and np.array_equal(attrs_a[attr], attrs_b[attr]),
            label+': aperture or native field-time convention differs: '+attr)
    end = min(float(reference['planes']['exit']['time_s'].max()),
              float(candidate['planes']['exit']['time_s'].max()))
    prefix = a['time_s'] <= end
    require(np.count_nonzero(prefix) > 1, 'No scientific field samples before the complete cohort exit')
    result['fields'] = {'clock': array_error(a['time_s'], b['time_s']),
        'full_output_samples': len(prefix), 'through_cohort_exit_samples': int(np.count_nonzero(prefix)),
        'cohort_cutoff_lab_s': end, 'native_B_time': str(attrs_a['magnetic_time_stagger'])}
    for field in ('electric_V_per_m', 'magnetic_T'):
        result['fields'][field] = {'full_output': array_error(a['fields'][field], b['fields'][field]),
            'through_cohort_exit': array_error(a['fields'][field][prefix], b['fields'][field][prefix])}
    result['field_regions'] = field_regions(reference, candidate, prefix)
    result['particle_response_differences'] = particle_response_differences(reference, candidate, result)
    return result


def replay(a, b):
    require(a['files'].keys() == b['files'].keys(), 'Replay HDF5 file set differs')
    datasets = 0
    for filename, (first, _) in a['files'].items():
        second = b['files'][filename][0]
        require(first.keys() == second.keys(), 'Replay dataset set differs: '+filename)
        for key in first:
            require(first[key].dtype == second[key].dtype and first[key].shape == second[key].shape and
                first[key].tobytes() == second[key].tobytes(), 'Replay scientific mismatch: '+filename+'/'+key)
            datasets += 1
    return {'scientific_datasets_bitwise_equal': datasets}


def analyze(args):
    root, output = args.root.resolve(), args.output.resolve()
    require(not output.exists(), 'Use a fresh review directory')
    plan = json.loads((root/'plan.json').read_text())
    with h5py.File(root/'inputs/held.h5') as stream:
        source = stream['particles/records'][:]
    paths = {name: root/name for name in plan['variants']}
    for name in ('fixed', 'moving', 'old-fixed'):
        supplied = getattr(args, name.replace('-', '_')+'_case')
        if supplied:
            supplied = supplied.resolve()
            other_plan = json.loads((supplied.parent/'plan.json').read_text())
            require(other_plan['ranks'] == plan['ranks'], 'Resumed same-layout comparison has different ranks')
            with h5py.File(supplied.parent/'inputs/held.h5') as stream:
                require(np.array_equal(source, stream['particles/records'][:]), 'Resumed case input differs')
            paths[name] = supplied
    cases = {name: inspect_case(path, plan, source) for name, path in paths.items()}
    report = {'format': 'aprl-moving-window-integration-review-v1', 'root': str(root),
        'scope': 'Functional structure checked; approximation metrics reported, not automatic physical acceptance or production speedup.',
        'cases': {name: case['summary'] for name, case in cases.items()}}
    if 'fixed' in cases and 'moving' in cases:
        report['fixed_to_moving'] = compare(cases['fixed'], cases['moving'], 'same-rank fixed-to-moving full-path approximation')
    if 'old-fixed' in cases and 'fixed' in cases:
        report['disabled_path_replay'] = replay(cases['old-fixed'], cases['fixed'])
    if args.other_root:
        other = args.other_root.resolve()
        other_plan = json.loads((other/'plan.json').read_text())
        with h5py.File(other/'inputs/held.h5') as stream:
            require(np.array_equal(source, stream['particles/records'][:]), 'Cross-layout inputs differ')
        report['other_layout'] = {'root': str(other), 'comparisons': {}}
        for name in set(cases).intersection(('fixed', 'moving')):
            other_case = inspect_case(other/name, other_plan, source)
            report['other_layout']['comparisons'][name] = compare(cases[name], other_case, name+' across layouts')
    output.mkdir(parents=True)
    save(output/'REVIEW.json', report)
    if not args.no_plot and 'fixed' in cases and 'moving' in cases:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
        first = cases['fixed']['files']['detectors/exit-field.h5'][0]
        second = cases['moving']['files']['detectors/exit-field.h5'][0]
        times = (first['time_s']-first['time_s'][0])*1e15
        for row, (field, unit) in enumerate((('electric_V_per_m', 'V/m'), ('magnetic_T', 'T'))):
            a, b = first['fields'][field], second['fields'][field]
            for label, data in [('fixed', a), ('moving', b)]:
                axes[row, 0].plot(times, np.sqrt(np.mean(data**2, axis=(1, 2, 3))), label=label)
            scale = np.sqrt(np.mean(a**2))
            axes[row, 1].plot(times, np.sqrt(np.mean((b-a)**2, axis=(1, 2, 3)))/scale)
            axes[row, 0].set_ylabel('Plane RMS '+unit)
            axes[row, 1].set_ylabel('Difference RMS / fixed whole-trace RMS')
            for axis in axes[row]:
                axis.set_xlabel('Lab time relative to first field sample (fs)')
                axis.grid(alpha=.25)
        axes[0, 0].legend()
        figure.suptitle('Full integration: fixed vs moving, native E and staggered B; not radiation-only')
        figure.savefig(output/'FIELD_COMPARISON.png', dpi=160)
        figure.savefig(output/'FIELD_COMPARISON.pdf')
        plt.close(figure)
    print(output/'REVIEW.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--other-root', type=Path)
    parser.add_argument('--fixed-case', type=Path, help='Completed fixed case from a separate frozen continuation')
    parser.add_argument('--moving-case', type=Path)
    parser.add_argument('--old-fixed-case', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--no-plot', action='store_true')
    analyze(parser.parse_args())
