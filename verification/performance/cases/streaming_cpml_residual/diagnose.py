#!/usr/bin/env python3
"""Offline time-window comparison for the completed residual-CPML campaign.

Read the final ten-case REVIEW.json and its original CSV files. No PIC is run,
no physical-error threshold is introduced, and numerical first inequality is
not treated as a physical onset. Frozen runner/reviewer files are not changed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


C = 299792458.0
E = ('Ex_V_m', 'Ey_V_m', 'Ez_V_m')
B = ('Bx_T', 'By_T', 'Bz_T')
CENTER = ('field_Ex_center_V_m', 'field_Ey_center_V_m', 'field_Ez_center_V_m')
CSV = {'history': 'history.csv', 'field': 'field-plane.csv',
       'initial': 'initial-particles.csv', 'release': 'release-particles.csv',
       'detector': 'detector-crossings.csv'}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


class Inputs:
    def __init__(self):
        self.hashes = {}

    def record(self, path):
        path = path.resolve()
        self.hashes[str(path)] = digest(path)
        return path

    def json(self, path):
        return json.loads(self.record(path).read_text())

    def table(self, path):
        a = np.atleast_1d(np.genfromtxt(self.record(path), names=True, delimiter=','))
        require(a.size > 0, 'Empty CSV: ' + str(path))
        require(all(np.isfinite(a[k]).all() for k in a.dtype.names),
                'Nonfinite CSV: ' + str(path))
        return a


def load_case(inputs, directory, completed):
    m = inputs.json(directory / 'metrics.json')
    require(m['completed'], 'Incomplete case: ' + str(directory))
    require(inputs.hashes[str((directory/'metrics.json').resolve())] == completed['metrics_sha256'],
            'Completed metrics changed: ' + str(directory))
    case = {'path': str(directory.resolve()), 'metrics': m}
    for key, name in CSV.items():
        a = inputs.table(directory/name)
        require(inputs.hashes[str((directory/name).resolve())] == completed['csv_sha256'][name],
                'Completed CSV changed: ' + str(directory/name))
        if key in ('initial', 'release', 'detector'):
            require(np.unique(a['id']).size == a.size, 'Duplicate particle ID')
            a = np.sort(a, order='id')
        elif key == 'field':
            a = np.sort(a, order=('step', 'pixel_j', 'pixel_i'))
        else:
            a = np.sort(a, order='step')
        case[key] = a
    return case


def vectors(table, columns):
    return np.column_stack([table[k] for k in columns])


def difference(reference, candidate):
    require(reference.shape == candidate.shape and reference.size, 'Unmatched or empty vector data')
    delta = candidate-reference
    source_norm = np.linalg.norm(reference.ravel())
    error_norm = np.linalg.norm(delta.ravel())
    return {
        'relative_l2': float(error_norm/source_norm) if source_norm else None,
        'vector_rms_difference': float(np.sqrt(np.mean(np.sum(delta*delta, axis=1)))),
        'maximum_absolute_component_difference': float(np.max(np.abs(delta))),
        'reference_vector_rms': float(np.sqrt(np.mean(np.sum(reference*reference, axis=1)))),
        'candidate_vector_rms': float(np.sqrt(np.mean(np.sum(candidate*candidate, axis=1)))),
        'difference_l2': float(error_norm), 'reference_l2': float(source_norm),
    }


def field_interval(full, candidate, mask):
    steps = candidate['field']['step'][mask]
    require(steps.size, 'Empty field interval')
    return {
        'included_field_steps': [int(v) for v in np.unique(steps)],
        'pixel_records_per_run': int(steps.size),
        'E_V_m': difference(vectors(full['field'][mask], E), vectors(candidate['field'][mask], E)),
        'cB_V_m': difference(C*vectors(full['field'][mask], B), C*vectors(candidate['field'][mask], B)),
    }


def first_event(rows, dt):
    row = rows[np.argmin(rows['time_box_s'])]
    return {'id': int(row['id']), 'loop_step': int(row['step']), 'fraction': float(row['fraction']),
            'continuous_step_from_time': float(row['time_box_s']/dt),
            'time_box_s': float(row['time_box_s'])}


def moments(case):
    p, initial = case['detector'], case['initial']
    weights = np.abs(p['charge_C']).astype(np.longdouble)
    weights /= weights.sum()
    def mean(a):
        return np.sum(weights*np.asarray(a, dtype=np.longdouble))
    def std(a):
        a = np.asarray(a, dtype=np.longdouble)
        return np.sqrt(mean((a-mean(a))**2))
    radius = np.hypot(std(p['x_m']), std(p['y_m']))
    initial_radius = np.hypot(std(initial['x_m']), std(initial['y_m']))
    # These fixtures prescribe one common initial proper velocity to all
    # particles: sigma(delta E_lab) is therefore also sigma(E_lab), without
    # subtracting a large common rest/drift energy numerically.
    for k in ('ux', 'uy', 'uz'):
        require(np.all(initial[k] == initial[k][0]), 'Nonuniform initial momentum needs absolute-energy variance')
    return {'mean_delta_energy_lab_eV': float(mean(p['delta_energy_lab_eV'])),
            'sigma_energy_lab_eV': float(std(p['delta_energy_lab_eV'])),
            'mean_delta_energy_box_eV': float(mean(p['delta_energy_box_eV'])),
            'rms_radius_m': float(radius),
            'radius_growth_percent': float((radius/initial_radius-1)*100)}


def compare(full, candidate):
    m, fm = candidate['metrics'], full['metrics']
    require(m['beta_box'] == fm['beta_box'] and m['dt_s'] == fm['dt_s'], 'Different control beta/time step')
    for key in ('step', 'pixel_i', 'pixel_j', 'time_E_box_s', 'time_B_box_s', 'z_box_m', 'x_m', 'y_m'):
        require(np.array_equal(full['field'][key], candidate['field'][key]), 'Mismatched plane sampling: '+key)
    require(np.array_equal(full['history']['step'], candidate['history']['step']), 'Mismatched central sample steps')
    dt = m['dt_s']
    release = first_event(candidate['release'], dt)
    full_release = first_event(full['release'], dt)
    detector = first_event(candidate['detector'], dt)
    pre_release_limit = min(release['continuous_step_from_time'], full_release['continuous_step_from_time'])
    last_candidate = float(np.max(candidate['detector']['time_box_s'])/dt)
    last_full = float(np.max(full['detector']['time_box_s'])/dt)
    last_limit = min(last_candidate, last_full)
    steps = candidate['field']['step']
    by_time = []
    for step in np.unique(steps):
        row = field_interval(full, candidate, steps == step)
        by_time.append({'step': int(step), 'E_V_m': row['E_V_m'], 'cB_V_m': row['cB_V_m']})
    central_delta = vectors(candidate['history'], CENTER)-vectors(full['history'], CENTER)
    central_steps = candidate['history']['step']
    before = central_steps < pre_release_limit
    anchors = []
    for step in (20, 23, 24, 26, 28, 32, 40, 48):
        mask = central_steps == step
        if np.any(mask):
            anchors.append({'step': step, 'candidate_minus_full_Ex_Ey_Ez_V_m': central_delta[mask][0].tolist()})
    removed = candidate['history'][candidate['history']['evicted'] > 0]
    candidate_moments, full_moments = moments(candidate), moments(full)
    return {
        'candidate_path': candidate['path'], 'full_path': full['path'],
        'beta_box': m['beta_box'], 'mode': m['mode'], 'guard_cells': m['guard_cells'],
        'initial_and_release_by_ID_numeric_bits_identical': all(
            candidate[k].dtype == full[k].dtype and candidate[k].tobytes() == full[k].tobytes()
            for k in ('initial', 'release')),
        'first_release': release, 'first_detector': detector,
        'first_removed_frame': None if removed.size == 0 else int(removed[0]['step']),
        'first_ballistic_observed_step': None,
        'first_ballistic_earliest_possible_loop_step': detector['loop_step']+1,
        'before_first_release': {
            'exclusive_upper_continuous_step': pre_release_limit,
            **field_interval(full, candidate, steps < pre_release_limit)},
        'until_common_last_particle_detection': {
            'last_candidate_continuous_step': last_candidate,
            'last_full_continuous_step': last_full,
            'inclusive_upper_continuous_step': last_limit,
            **field_interval(full, candidate, steps <= last_limit)},
        'whole_recorded_interval': field_interval(full, candidate, np.ones(steps.size, dtype=bool)),
        'before_release_central_E_maximum_absolute_component_difference_V_m': float(np.max(np.abs(central_delta[before]))),
        'central_E_at_fixed_anchor_steps': anchors,
        'field_by_time': by_time,
        'candidate_detector_moments': candidate_moments, 'full_detector_moments': full_moments,
        'candidate_minus_full_detector_moments': {key: candidate_moments[key]-full_moments[key] for key in full_moments},
    }


def plot_onset(comparisons, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, layout='constrained')
    for column, beta in enumerate((0.0, 0.23)):
        for mode, color, marker in [('cpml', '#b2472c', 's'), ('residual', '#20639b', 'o')]:
            item = next(r for r in comparisons if r['beta_box'] == beta and r['mode'] == mode and r['guard_cells'] == 24)
            points = [r for r in item['field_by_time'] if r['step'] <= 48]
            for row, component in enumerate(('E_V_m', 'cB_V_m')):
                values = [r[component]['vector_rms_difference'] for r in points]
                axes[row, column].plot([r['step'] for r in points],
                    [v if v > 0 else math.nan for v in values],
                    marker=marker, markersize=3.5, color=color, linewidth=1.4,
                    label=('Total-field CPML' if mode == 'cpml' else 'Residual-input CPML'))
        representative = next(r for r in comparisons if r['beta_box'] == beta and r['mode'] == 'residual' and r['guard_cells'] == 24)
        release = representative['first_release']['continuous_step_from_time']
        for row in range(2):
            ax = axes[row, column]
            ax.axvline(release, color='#617527', linestyle='--', linewidth=1.1,
                       label='First release: %.3f' % release)
            if beta == 0.0:
                ax.axvline(23, color='#777777', linestyle=':', linewidth=1,
                           label='Prior total-field central E difference: n=23')
            ax.set_yscale('log')
            ax.set_xlim(0, 49)
            ax.grid(True, which='both', alpha=.17)
            ax.legend(loc='best', fontsize=7)
            ax.set_ylabel(('E' if row == 0 else 'cB')+' plane RMS difference (V/m)')
        axes[0, column].set_title('Guard 24, particle box beta = %s' % beta)
        axes[1, column].set_xlabel('Complete box-frame step n')
    fig.suptitle('Same-binary full-domain control: early plane-field differences')
    fig.text(.01, -.035,
        'Plane data only, sampled every 4 steps; exact zeros omitted on log axes. No onset threshold is applied.\n'
        'The n=23 line marks the prior central sample, not a plane point. Native E(n), B(n-1/2); cB is a unit conversion.', fontsize=8)
    directory.mkdir(parents=True, exist_ok=True)
    outputs = []
    for extension in ('png', 'pdf'):
        path = directory/('causal_onset.'+extension)
        fig.savefig(path, dpi=180, bbox_inches='tight')
        outputs.append(str(path.resolve()))
    plt.close(fig)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', type=Path, required=True)
    args = parser.parse_args()
    inputs = Inputs()
    official = inputs.json(args.review)
    require(official['cohort_charge_and_input_consistency_status'] == 'PASS' and
            official['legacy_csv_replay_status'] == 'PASS', 'Primary campaign review has not passed')
    require(official['coverage']['complete_planned_ten_case_set'], 'Final ten-case review is required')
    campaigns = sorted({Path(entry['path']).parent for entry in official['cases']})
    completed, binary_hashes = {}, set()
    for campaign in campaigns:
        plan = inputs.json(campaign/'plan.json')
        binary_hashes.add(plan['binary_sha256'])
        done = inputs.json(campaign/'COMPLETED.json')
        require(done['cohort_and_charge_gates'] == 'COMPLETED', 'Campaign not complete')
        for row in done['cases']:
            completed[str((campaign/row['case']).resolve())] = row
        inputs.record(campaign/'source-files.json')
        for name in ('streaming_cpml_boundary.h', 'streaming_cpml_boundary.cpp', 'streaming_cpml_bunch_probe.cpp'):
            inputs.record(campaign/'source/verification'/name)
    require(len(binary_hashes) == 1, 'Candidate and controls did not use one binary')
    cases = [load_case(inputs, Path(entry['path']), completed[entry['path']]) for entry in official['cases']]
    full = {case['metrics']['beta_box']: case for case in cases if case['metrics']['mode'] == 'full'}
    require(set(full) == {0.0, 0.23}, 'Missing full-domain controls')
    comparisons = [compare(full[case['metrics']['beta_box']], case) for case in cases if case['metrics']['mode'] != 'full']
    require(all(c['initial_and_release_by_ID_numeric_bits_identical'] for c in comparisons),
            'Prescribed initial/release cohort changed')
    report = {
        'format': 'aprl-streaming-cpml-residual-time-diagnostic-v1',
        'script_path': str(Path(__file__).resolve()), 'script_sha256': digest(Path(__file__)),
        'input_sha256': inputs.hashes, 'binary_sha256': next(iter(binary_hashes)),
        'primary_review': str(args.review.resolve()), 'comparisons': comparisons,
        'plots': plot_onset(comparisons, args.plots),
        'limitations': [
            'Offline analysis only. No new PIC, numerical-method change, or physical-error acceptance threshold.',
            'All primary comparisons use the independent same-beta full control from this binary.',
            'Release/detection interval limits come from both actual event files; plane samples are common native samples every four steps.',
            'No first-nonzero numerical difference is called a physical onset; floating-point differences are retained and shown without an invented noise floor.',
            'Before-first-release excludes particle release and post-detector outflow, but moving held sources and moving boundaries still operate.',
            'Until the last particle detection can include earlier particles that are already ballistic or removed. It is not an outflow-free control.',
            'No exact first ballistic transition was logged. Only its detector-gate lower bound is available.',
            'E/B retain staggered native times. cB supplies comparable units, not simultaneous lab energy flux or pure radiation.',
            'Relative B errors can be misleading for a nearly zero reference B; use the reported absolute cB differences.',
            'Differences from separate time intervals must not be subtracted as a reflection coefficient or isolated mechanism contribution.',
            'One bunch length and 176 steps do not establish long-time stability, weak-radiation accuracy, asymptotic scaling, or production acceptance.',
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(args.output.resolve())
    for row in comparisons:
        print(row['beta_box'], row['mode'], row['guard_cells'],
              'pre-release E L2', row['before_first_release']['E_V_m']['relative_l2'],
              'through-last-detection E L2', row['until_common_last_particle_detection']['E_V_m']['relative_l2'],
              'whole E L2', row['whole_recorded_interval']['E_V_m']['relative_l2'])


if __name__ == '__main__':
    main()
