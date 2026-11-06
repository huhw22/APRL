"""Analyze only closed physical-electron FEL outputs in bounded memory.

Source IDs identify immutable physical electrons. Local bins are the original
generation-time bins, not propagated time slices. Native substep gamma_lab is
used as the affine four-momentum scalar, never reconstructed from interpolated u.
The frozen report_discrete.py/review_discrete.py pair must accompany this file.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import traceback

import h5py
import numpy as np
import yaml


C = np.longdouble('299792458')
MEC2_EV = np.longdouble('510998.95069')
BLOCK = 131072


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(item) for item in value]
    return value


def save(path, value):
    Path(path).write_text(json.dumps(native(value), indent=2, allow_nan=False) + '\n')


def finite(array, where):
    if array.dtype.names:
        require(all(np.isfinite(array[name]).all() for name in array.dtype.names
                    if array[name].dtype.kind in 'fc'), 'Nonfinite records: ' + str(where))
    elif array.dtype.kind in 'fc':
        require(np.isfinite(array).all(), 'Nonfinite values: ' + str(where))


def phase(time, anchor, wavelength):
    cycle = (np.asarray(time, dtype=np.longdouble) - anchor) * C / np.longdouble(wavelength)
    angle = np.asarray(np.remainder(cycle, np.longdouble(1)) *
                       np.longdouble(2 * np.pi), dtype=float)
    return np.exp(1j * angle)


class Moment:
    """Long-double moments of b, native gamma, and per-ID gamma changes."""

    def __init__(self):
        self.s = np.zeros(12, dtype=np.longdouble)
        self.count = 0

    def add(self, weight, harmonic, gamma, delta):
        if not len(weight):
            return
        w = np.asarray(weight, dtype=np.longdouble)
        g = np.asarray(gamma, dtype=np.longdouble)
        d = np.asarray(delta, dtype=np.longdouble)
        qr = np.asarray(harmonic.real, dtype=np.longdouble)
        qi = np.asarray(harmonic.imag, dtype=np.longdouble)
        self.count += len(w)
        self.s += [np.sum(w), np.sum(w * w), np.sum(w * qr), np.sum(w * qi),
                   np.sum(w * g), np.sum(w * g * g), np.sum(w * d), np.sum(w * d * d),
                   np.sum(w * d * qr), np.sum(w * d * qi),
                   np.sum(w * g * qr), np.sum(w * g * qi)]

    def result(self):
        s = self.s
        require(s[0] > 0, 'Empty BF scope')
        b = complex(float(s[2] / s[0]), float(s[3] / s[0]))
        mean = s[6] / s[0]
        variance = max(np.longdouble(0), s[7] / s[0] - mean * mean)
        h = complex(float(s[8] / s[0] - mean * s[2] / s[0]),
                    float(s[9] / s[0] - mean * s[3] / s[0]))
        noise = s[1] / (s[0] * s[0])
        return dict(records=self.count, represented_electrons=float(s[0]),
                    b_real=b.real, b_imag=b.imag, BF=abs(b),
                    random_phase_BF_rms=float(np.sqrt(noise)),
                    random_phase_BF_power=float(noise),
                    BF_power_over_noise=float(abs(b)**2 / noise),
                    effective_electrons_per_macro=float(s[1] / s[0]),
                    mean_gamma=float(s[4] / s[0]), mean_delta_gamma=float(mean),
                    std_delta_gamma=float(np.sqrt(variance)),
                    mean_delta_energy_eV=float(mean * MEC2_EV),
                    std_delta_energy_eV=float(np.sqrt(variance) * MEC2_EV),
                    centered_delta_energy_h1_eV=float(abs(h) * MEC2_EV))


def input_data(plan, cohort):
    """Keep only gamma0 (8N bytes); no 64N-byte source-record cache."""
    n = int(plan['expected_particles'])
    first_id = int(cohort.get('parent_crop_original_row_start', 0)) + 1
    gamma0 = np.empty(n, dtype=np.float64)
    tmin, tmax = math.inf, -math.inf
    with h5py.File(plan['input'], 'r') as handle:
        data = handle['particles/records']
        if 'generation' in handle:
            require(int(handle['generation/complete'][()]) == 1, 'Physical source is not closed')
            counts = handle['generation/bin_counts'][:].astype(np.int64)
        else:
            require('selection' in handle and int(handle['selection'].attrs['complete']) == 1,
                    'Physical cropped source is not closed')
            counts = handle['selection/selected_records_per_original_generation_bin'][:].astype(np.int64)
        require(len(data) == n and int(counts.sum()) == n, 'Source record/bin count differs')
        require(np.all(counts > 0), 'Empty physical generation bin')
        for start in range(0, n, BLOCK):
            records = data[start:start + BLOCK]
            finite(records, plan['input'])
            require(np.array_equal(records['source_id'], np.arange(start + first_id,
                    start + len(records) + first_id, dtype=np.uint64)), 'Source IDs are not contiguous physical IDs')
            require(np.all(records['macro_weight'] == 1.), 'Source is not one electron per particle')
            proper = records['proper_velocity']
            gamma0[start:start + len(records)] = np.sqrt(1. + np.sum(proper * proper, axis=1))
            tmin = min(tmin, float(records['arrival_time_offset_s'].min()))
            tmax = max(tmax, float(records['arrival_time_offset_s'].max()))
    require(np.array_equal(counts, np.asarray(plan['source_metadata']['original_generation_bin_counts'])),
            'Source generation counts differ from prepared plan')
    edges = np.r_[0, np.cumsum(counts)] + first_id
    first = int(plan['first_original_slice'])
    core_start = int(plan['core_first_slice_inclusive']) - first
    core_stop = int(plan['core_end_slice_exclusive']) - first
    require(0 <= core_start < core_stop <= len(counts), 'Invalid fixed core bin range')
    require(int(edges[core_start]) == cohort['core_original_source_id_first'] and
            int(edges[core_stop]) - 1 == cohort['core_original_source_id_last'],
            'Fixed core IDs do not match physical generation bins')
    require(int(edges[core_stop] - edges[core_start]) == plan['expected_core_records'],
            'Fixed core particle count differs')
    return gamma0, edges, core_start, core_stop, np.longdouble(tmin), dict(
        records=n, bin_count=len(counts), core_bin_count=core_stop - core_start,
        arrival_time_min_s=tmin, arrival_time_max_s=tmax,
        cache_bytes=int(gamma0.nbytes), first_source_id=first_id,
        source_id_contract='Contiguous original physical IDs; weight1.')


def closed_group(handle, groupname, count=None, ranks=None):
    group = handle[groupname]
    require(int(group['complete'][()]) == 1, 'Unclosed detector group: ' + groupname)
    for counter, values in [('committed_records', 'records'), ('committed_samples', 'time_s')]:
        if counter in group:
            require(int(group[counter][()]) == len(group[values]), 'Committed detector length differs')
    if count is not None:
        require(len(group['records']) == count, 'Incomplete full source through detector')
    if ranks is not None:
        if groupname in ['particle_plane', 'particle_substep_plane']:
            require('mpi_size' in group.attrs and int(group.attrs['mpi_size']) == ranks,
                    'Particle detector MPI size differs')
        elif 'mpi_size' in group.attrs:
            require(int(group.attrs['mpi_size']) == ranks, 'Detector MPI size differs')
    return group


def validate_ids(records, seen, count, path, first_id=1):
    ids = records['source_id']
    require(np.all((ids >= first_id) & (ids < first_id + count)), 'Unknown source ID: ' + str(path))
    index = ids.astype(np.int64) - first_id
    require(len(np.unique(index)) == len(index) and not np.any(seen[index]),
            'Repeated source ID through detector: ' + str(path))
    seen[index] = True
    require(np.all(records['weight'] == 1.), 'Detector particle weight differs from1')
    require(np.all(records['crossing_direction'] == 1), 'Unexpected backward detector crossing')
    return index


def plane_statistics(path, plan, gamma0, edges, core_start, core_stop, tmin, z):
    count = len(gamma0)
    anchor = tmin + np.longdouble(z) / C
    moments = {key: Moment() for key in ['full', 'core', 'core_gaussian_100um']}
    bins = [Moment() for _ in range(core_stop - core_start)]
    seen = np.zeros(count, dtype=np.bool_)
    with h5py.File(path, 'r') as handle:
        group = closed_group(handle, 'particle_substep_plane', count, plan['ranks'])
        require(abs(float(group.attrs['plane_z_m']) - z) < 1e-12, 'Wrong lab detector location')
        data = group['records']
        for start in range(0, count, BLOCK):
            records = data[start:start + BLOCK]
            finite(records, path)
            index = validate_ids(records, seen, count, path, int(edges[0]))
            harmonic = phase(records['time_s'], anchor, plan['wavelength_m'])
            gamma = records['gamma_lab']
            delta = gamma.astype(np.longdouble) - gamma0[index].astype(np.longdouble)
            weight = records['weight']
            original_bin = np.searchsorted(edges, records['source_id'], side='right') - 1
            core = (original_bin >= core_start) & (original_bin < core_stop)
            moments['full'].add(weight, harmonic, gamma, delta)
            moments['core'].add(weight[core], harmonic[core], gamma[core], delta[core])
            radius2 = np.sum(records['position_m'][:, :2]**2, axis=1)
            transverse_weight = weight * np.exp(-radius2 / (2 * (100e-6)**2))
            moments['core_gaussian_100um'].add(transverse_weight[core], harmonic[core],
                                               gamma[core], delta[core])
            # Stable grouping avoids37 complete boolean scans of the whole block.
            core_indices = np.flatnonzero(core)
            order = np.argsort(original_bin[core_indices], kind='stable')
            ordered = core_indices[order]
            groups = original_bin[ordered]
            boundaries = np.searchsorted(groups, np.arange(core_start, core_stop + 1))
            for k, moment in enumerate(bins):
                selected = ordered[boundaries[k]:boundaries[k + 1]]
                moment.add(weight[selected], harmonic[selected], gamma[selected], delta[selected])
    require(np.all(seen), 'Missing source ID through particle detector: ' + str(path))
    result = {key: value.result() for key, value in moments.items()}
    local = [value.result() for value in bins]
    electrons = np.asarray([value['represented_electrons'] for value in local])
    result['core']['local_RMS_BF_1lambda'] = float(np.sqrt(np.sum(electrons *
                                                np.array([value['BF']**2 for value in local])) / electrons.sum()))
    result['core']['local_random_phase_BF_rms_1lambda'] = float(np.sqrt(np.sum(electrons *
                                        np.array([value['random_phase_BF_power'] for value in local])) / electrons.sum()))
    result['core']['local_RMS_definition'] = (
        f'Charge-weighted RMS across{len(local)} fixed parent generation bins, '
        'each one optical wavelength; original IDs, not propagated time bins.')
    for k, value in enumerate(local):
        value['original_optical_bin'] = int(plan['first_original_slice']) + core_start + k
        value['source_id_first'] = int(edges[core_start + k])
        value['source_id_last'] = int(edges[core_start + k + 1]) - 1
    result['core_bins'] = local
    return result


def verify_legacy_and_fields(directory, plan):
    verified = []
    n = int(plan['expected_particles'])
    cohort = json.loads(Path(plan['cohort']).read_text())
    first_id = int(cohort.get('parent_crop_original_row_start', 0)) + 1
    expected_legacy = {plane['name'] + '-particles.h5': plane['z'] for plane in plan['particle_planes']}
    for path in sorted((directory / 'detectors').glob('*.h5')):
        if path.name.endswith('-particles-substeps.h5'):
            continue
        with h5py.File(path, 'r') as handle:
            roots = list(handle)
            require(len(roots) == 1, 'Unexpected detector group count: ' + str(path))
            group = closed_group(handle, roots[0], n if path.name in expected_legacy else None,
                                 plan['ranks'])
            if path.name in expected_legacy:
                require(roots[0] == 'particle_plane', 'Wrong legacy particle group')
                require(abs(float(group.attrs['plane_z_m']) - expected_legacy[path.name]) < 1e-12,
                        'Legacy particle plane location differs')
                seen = np.zeros(n, dtype=np.bool_)
                for start in range(0, n, BLOCK):
                    records = group['records'][start:start + BLOCK]
                    finite(records, path)
                    validate_ids(records, seen, n, path, first_id)
                require(np.all(seen), 'Missing legacy plane source IDs')
                del seen
            else:
                def visit(name, dataset):
                    if not isinstance(dataset, h5py.Dataset):
                        return
                    if not dataset.shape:
                        finite(np.asarray(dataset[()]), path)
                        return
                    block = 4 if name == 'fields' else BLOCK
                    for start in range(0, len(dataset), block):
                        finite(dataset[start:start + block], path)
                group.visititems(visit)
        verified.append(str(path.relative_to(directory)))
    require(set(expected_legacy).issubset({Path(value).name for value in verified}),
            'Missing configured legacy particle file')
    return verified


def verify_detector_paths(directory):
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    detectors = card['detectors']
    expected = {plane['name'] + '.h5' for key in ['particle_planes', 'field_planes']
                for plane in detectors.get(key, [])}
    expected.update(plane['name'] + '-energy.h5' for plane in detectors.get('energy_planes', []))
    if detectors.get('substep_particle_planes'):
        expected.update(plane['name'] + '-substeps.h5' for plane in detectors['particle_planes'])
    actual = {path.name for path in (directory / 'detectors').glob('*.h5')}
    require(actual == expected, 'Configured detector file set differs: missing=' +
            str(sorted(expected - actual)) + '; extra=' + str(sorted(actual - expected)))


def ledger_and_coupling(directory, plan, output):
    import report_discrete
    rows, sums, ranks, steps = report_discrete.load(directory / 'ledger',
                                                   directory / 'ledger/energy-ledger.h5')
    require(ranks == plan['ranks'] and abs(steps - plan['expected_steps']) <= 2,
            'MPI layout or final step differs')
    index, quantities, summary = report_discrete.calculate(rows, sums, steps)
    denominator = float(np.max(np.abs(np.asarray(rows['delta_particle_kinetic_box_J'], dtype=np.longdouble))))
    denominator_source = 'max absolute stable same-cohort delta_particle_kinetic_box_J over full record'
    if denominator == 0:
        denominator = float(max(np.max(np.abs(sums['native_begin_J'])),
                                np.max(np.abs(sums['native_end_J']))))
        denominator_source = 'zero particle change: peak native field energy over complete audit record'
    require(denominator > 0, 'All-zero ledger does not define relative normalization')
    names = list(quantities)
    for suffix, relative in [('absolute', False), ('relative', True)]:
        with (output / ('energy_ledger_' + suffix + '.csv')).open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['step'] + names)
            writer.writerows([[int(step)] + [str(value[k] * (100 / denominator if relative else 1))
                                             for value in quantities.values()]
                              for k, step in enumerate(index)])
    summary['normalization'] = dict(denominator_J=denominator, definition=denominator_source,
                                     unit='percent', sign='Signed; same fixed denominator for every curve and time.')
    summary['final_percent_of_D'] = {name: float(value[-1] * 100 / denominator)
                                     for name, value in quantities.items()}
    summary['physical_conservation_accepted'] = False
    summary['interpretation'] = ('Discrete numerical accounting includes time/Cowan/curl-exchange channels; '
                                  'these are not emitted radiation or physical external work.')
    paths = sorted((directory / 'reference').glob('source-coupling-rank-*.csv'))
    require([path.name for path in paths] ==
            [f'source-coupling-rank-{k:05d}.csv' for k in range(ranks)], 'Missing/extra coupling histories')
    trials = committed = 0
    maximum_force = maximum_work = maximum_iterations = 0.
    for path in paths:
        require(path.with_suffix('.complete').read_text().split() == ['1', str(steps)],
                'Unclosed conservative-source history: ' + str(path))
        records = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        require(np.array_equal(records['step'], np.arange(1, steps + 1)), 'Missing/repeated coupling step')
        finite(records, path)
        require(np.all(records['self_force_defect'] < 1e-13) and
                np.all(records['relative_fixed_point_work_defect'] < 1e-8),
                'Conservative-source iteration did not converge')
        trials += int(records['trial_pushes'].sum())
        committed += int(records['committed_pushes'].sum())
        maximum_force = max(maximum_force, float(np.max(records['self_force_defect'])))
        maximum_work = max(maximum_work, float(np.max(records['relative_fixed_point_work_defect'])))
        if 'iterations' in records.dtype.names:
            maximum_iterations = max(maximum_iterations, float(np.max(records['iterations'])))
    coupling = dict(trial_pushes=trials, committed_pushes=committed,
                    maximum_self_force_defect=maximum_force,
                    maximum_relative_fixed_point_work_defect=maximum_work,
                    maximum_iterations=maximum_iterations, ranks=ranks, steps=steps)
    return steps, summary, coupling


def figures_and_csv(output, plan, planes):
    rows = []
    for plane in planes:
        for scope in ['full', 'core', 'core_gaussian_100um']:
            rows.append(dict(plane=plane['name'], z_m=plane['z_m'], scope=scope,
                        **{key: value for key, value in plane['scopes'][scope].items()
                           if not isinstance(value, str)}))
    names = list(dict.fromkeys(key for row in rows for key in row))
    with (output / 'bf_curves.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)
    local_rows = [dict(plane=plane['name'], z_m=plane['z_m'], **value)
                  for plane in planes for value in plane['scopes']['core_bins']]
    with (output / 'local_bin_curves.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(local_rows[0]))
        writer.writeheader()
        writer.writerows(local_rows)
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-binomial3-fel-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), layout='constrained')
    z = np.asarray([plane['z_m'] for plane in planes])
    core = [plane['scopes']['core'] for plane in planes]
    bf = np.asarray([value['BF'] for value in core])
    noise = np.asarray([value['random_phase_BF_rms'] for value in core])
    local = np.asarray([value['local_RMS_BF_1lambda'] for value in core])
    local_noise = np.asarray([value['local_random_phase_BF_rms_1lambda'] for value in core])
    energy_harmonic = np.asarray([value['centered_delta_energy_h1_eV'] for value in core])
    for ax, values in [(axes[0, 0], bf), (axes[0, 1], bf / noise),
                       (axes[1, 0], local), (axes[1, 1], energy_harmonic)]:
        ax.plot(z, values, '-o', color='#005A8D', markerfacecolor='white',
                linewidth=1.8, markersize=4, label='Fixed central37-wavelength cohort')
        ax.grid(alpha=.25)
        ax.set_xlabel('Fixed laboratory plane z (m)')
    axes[0, 0].plot(z, noise, '--', color='#555555', linewidth=1.2, label='Random-phase noise RMS')
    axes[0, 1].axhline(1., linestyle='--', color='#555555', linewidth=1.2)
    axes[1, 0].plot(z, local_noise, '--', color='#555555', linewidth=1.2,
                    label='Local random-phase noise RMS')
    titles = ['Global bunching of fixed core', 'Global bunching / noise RMS',
              'Local RMS bunching:37 fixed one-wavelength bins',
              'Centered energy-change harmonic at resonance']
    labels = ['|b|', '|b| / noise RMS', 'Charge-weighted RMS |b|', 'Energy modulation (eV)']
    for ax, title, label in zip(axes.flat, titles, labels):
        ax.set_title(title, fontsize=10.5)
        ax.set_ylabel(label)
    axes[0, 0].set_yscale('log')
    axes[1, 0].set_yscale('log')
    axes[0, 0].legend(fontsize=8, loc='best')
    axes[1, 0].legend(fontsize=8, loc='best')
    fig.suptitle('37-period FEL: independent one-electron loading, matched three-point F3 kernel\n'
                 'Energy-conservative particle coupling; noise lines are expectations, not subtracted backgrounds',
                 fontsize=10.5)
    fig.savefig(output / 'bf_growth.png', dpi=240)
    fig.savefig(output / 'bf_growth.pdf')
    plt.close(fig)


def analyze_case(directory):
    directory = Path(directory).resolve()
    require((directory / 'PIC_COMPLETE.json').exists(), 'PIC must finish and close HDF5 before analysis')
    closed = json.loads((directory / 'PIC_COMPLETE.json').read_text())
    require(closed.get('complete') is True, 'PIC completion is not successful')
    frozen = json.loads((directory.parent / 'FROZEN_INPUTS.json').read_text())
    for path, expected in frozen.items():
        digest = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
                digest.update(block)
        require(digest.hexdigest() == expected, 'Frozen identity differs before analysis: ' + path)
    plan = json.loads((directory / 'run-plan.json').read_text())
    cohort = json.loads(Path(plan['cohort']).read_text())
    verify_detector_paths(directory)
    output = directory.parent / 'analysis'
    output.mkdir(exist_ok=True)
    require(not (output / 'ANALYSIS_COMPLETE.json').exists(), 'Refusing to overwrite completed analysis')
    gamma0, edges, core_start, core_stop, tmin, input_summary = input_data(plan, cohort)
    planes = []
    for plane in plan['particle_planes']:
        path = directory / 'detectors' / (plane['name'] + '-particles-substeps.h5')
        values = plane_statistics(path, plan, gamma0, edges, core_start, core_stop, tmin, plane['z'])
        planes.append(dict(name=plane['name'], z_m=plane['z'], scopes=values))
        print(f'Closed particle plane analyzed: {plane["name"]}', flush=True)
    del gamma0
    verified = verify_legacy_and_fields(directory, plan)
    steps, ledger, coupling = ledger_and_coupling(directory, plan, output)
    figures_and_csv(output, plan, planes)
    result = dict(complete=True, program_integrity_pass=True, physical_acceptance=False,
        name=plan['name'], target_weight=1, actual_weight=1., selection_seed=plan['selection_seed'],
        steps=steps, ranks=plan['ranks'], planes=planes, input=input_summary,
        other_closed_scientific_files=verified, strict_ledger=ledger, coupling=coupling,
        coupling_trial_pushes=coupling['trial_pushes'],
        coupling_committed_pushes=coupling['committed_pushes'],
        shape='linear-binomial3', shape_fourier_response='F3(theta)=cos(theta/2)^2',
        definition='b=sum(weight*exp(i*2pi*c*(time_s-t_anchor)/lambda))/sum(weight); fixed original source IDs.',
        gamma_definition='Native affine endpoint gamma_lab, minus initial source gamma by unchanged source_id.',
        energy_unit_conversion_eV_per_gamma=str(MEC2_EV),
        interpretation='One physical random seed; complete current37-period benchmark. This is not statistical/theoretical FEL acceptance or a pure-radiation energy measurement.')
    save(directory / 'RESULT.json', result)
    save(output / 'SUMMARY.json', result)
    save(output / 'ANALYSIS_COMPLETE.json', dict(complete=True, program_integrity_pass=True,
         physical_acceptance=False, planes=len(planes), source_records=input_summary['records']))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    if (args.directory / 'COMPLETE.json').exists() and (args.directory.parent / 'analysis/ANALYSIS_COMPLETE.json').exists():
        require(json.loads((args.directory / 'COMPLETE.json').read_text()).get('analysis_complete') is True,
                'Unexpected existing completion record')
        print('Completed analysis already exists; retaining original outputs.')
        return
    try:
        result = analyze_case(args.directory)
        save(args.directory / 'COMPLETE.json', dict(complete=True, program_integrity_pass=True,
             analysis_complete=True, physical_acceptance=False))
        print(json.dumps(dict(complete=True, planes=len(result['planes']), steps=result['steps'],
                              ranks=result['ranks'], physical_acceptance=False), indent=2))
    except Exception as error:
        save(args.directory / 'ANALYSIS_FAILED.json', dict(complete=False, error=str(error),
             traceback=traceback.format_exc(), physical_acceptance=False))
        raise


if __name__ == '__main__':
    main()
