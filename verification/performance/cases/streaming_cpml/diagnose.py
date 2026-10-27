#!/usr/bin/env python3
"""Read existing CSV only: boundary-onset chronology and historical full controls.

This does not run PIC, impose a new physics threshold, or infer unrecorded
ballistic event times. Run after review.py has audited the complete campaign.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct


C = 299792458.0
E = ('Ex_V_m', 'Ey_V_m', 'Ez_V_m')
B = ('Bx_T', 'By_T', 'Bz_T')
CENTER = tuple('field_' + key.replace('_V_m', '_center_V_m') for key in E)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(value):
    return struct.pack('>d', float(value))


class Inputs:
    def __init__(self):
        self.hashes = {}

    def record(self, path):
        path = path.resolve()
        self.hashes[str(path)] = sha256(path)
        return path

    def table(self, path):
        path = self.record(path)
        with path.open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        require(rows, 'Empty CSV: ' + str(path))
        for row in rows:
            require(all(math.isfinite(float(value)) for value in row.values()),
                    'Nonfinite CSV: ' + str(path))
        return rows

    def json(self, path):
        return json.loads(self.record(path).read_text())


def case(inputs, directory):
    result = {'path': str(directory.resolve()),
              'metrics': inputs.json(directory / 'metrics.json')}
    for key, name in [('history', 'history.csv'), ('field', 'field-plane.csv'),
                      ('initial', 'initial-particles.csv'),
                      ('release', 'release-particles.csv'),
                      ('detector', 'detector-crossings.csv')]:
        result[key] = inputs.table(directory / name)
    require(result['metrics']['completed'], 'Incomplete case: ' + str(directory))
    return result


def by_id(rows):
    result = {int(row['id']): row for row in rows}
    require(len(result) == len(rows), 'Duplicate particle ID')
    return result


def cohort_match(a, b):
    a, b = by_id(a), by_id(b)
    require(a.keys() == b.keys(), 'Different particle ID sets')
    require(all(a[i].keys() == b[i].keys() for i in a), 'Different event columns')
    return {
        'same_ids': True,
        'all_numeric_double_bits_identical': all(
            bits(a[i][key]) == bits(b[i][key]) for i in a for key in a[i]),
        'all_canonical_csv_values_identical': all(a[i] == b[i] for i in a),
    }


def difference(a, b, columns, scale=1.0):
    require(len(a) == len(b), 'Different comparison row counts')
    pairs = [(float(x[k]) * scale, float(y[k]) * scale)
             for x, y in zip(a, b) for k in columns]
    source2 = math.fsum(x*x for x, _ in pairs)
    error2 = math.fsum((y-x)**2 for x, y in pairs)
    return {
        'relative_l2': math.sqrt(error2/source2) if source2 else None,
        'scalar_rms_absolute_difference': math.sqrt(error2/len(pairs)),
        'vector_rms_absolute_difference': math.sqrt(error2/len(a)),
        'maximum_absolute_difference': max(abs(y-x) for x, y in pairs),
        'reference_l2': math.sqrt(source2), 'difference_l2': math.sqrt(error2),
        'numeric_double_bits_identical': all(
            bits(x[k]) == bits(y[k]) for x, y in zip(a, b) for k in columns),
    }


def field_groups(a, b):
    key = lambda r: (int(r['step']), int(r['pixel_j']), int(r['pixel_i']))
    a, b = sorted(a, key=key), sorted(b, key=key)
    require([key(r) for r in a] == [key(r) for r in b], 'Different detector sampling')
    groups = {}
    for x, y in zip(a, b):
        group = groups.setdefault(int(x['step']), ([], []))
        group[0].append(x)
        group[1].append(y)
    return a, b, groups


def first_event(rows):
    r = min(rows, key=lambda x: int(x['step']) + float(x['fraction']))
    return {'particle_id': int(r['id']), 'loop_step': int(r['step']),
            'fraction': float(r['fraction']),
            'continuous_step': int(r['step']) + float(r['fraction']),
            'time_box_s': float(r['time_box_s']),
            'time_lab_s': float(r['time_lab_s'])}


def chronology(full, candidate):
    a, b = full['history'], candidate['history']
    require([int(r['step']) for r in a] == [int(r['step']) for r in b],
            'History steps differ')
    first_center = None
    for x, y in zip(a, b):
        delta = [float(y[k])-float(x[k]) for k in CENTER]
        if any(v != 0.0 for v in delta):
            first_center = {'step': int(y['step']), 'time_box_s': float(y['time_box_s']),
                            'candidate_minus_full_Ex_Ey_Ez_V_m': delta}
            break
    _, _, groups = field_groups(full['field'], candidate['field'])
    first_plane = None
    last_equal_plane = None
    by_time = []
    for step, (x, y) in groups.items():
        item = {'step': step, 'E_V_m': difference(x, y, E),
                'cB_V_m': difference(x, y, B, C)}
        by_time.append(item)
        if (item['E_V_m']['maximum_absolute_difference'] or
                item['cB_V_m']['maximum_absolute_difference']):
            if first_plane is None:
                first_plane = item
        elif first_plane is None:
            last_equal_plane = step
    release, detector = first_event(candidate['release']), first_event(candidate['detector'])
    eviction = next(({'step': int(r['step']), 'time_box_s': float(r['time_box_s']),
                      'cumulative_evicted': int(r['evicted'])}
                     for r in b if int(r['evicted'])), None)
    prior = [r for r in b if int(r['step']) < release['continuous_step']]
    gauss_rows = [{k: float(r[k]) for k in ('step', 'origin_z_m',
                  'total_gauss_before_shift', 'total_gauss_after_shift')}
                 for r in b[:9]]
    initial_origin = float(b[0]['origin_z_m'])
    dz = candidate['metrics']['dt_s'] * C
    pml = candidate['metrics']['cpml_z_cells']
    nz = candidate['metrics']['window_nz']
    z = [float(r['z_box_m']) for r in candidate['initial']]
    initial_counts = {
        'lower_cpml_particle_centres': sum(initial_origin <= v < initial_origin+pml*dz for v in z),
        'nonabsorbing_particle_centres': sum(initial_origin+pml*dz <= v <= initial_origin+(nz-pml)*dz for v in z),
        'upper_cpml_particle_centres': sum(initial_origin+(nz-pml)*dz < v <= initial_origin+nz*dz for v in z),
        'active_including_admission_lookahead_and_collar': int(b[0]['active']),
        'scope': 'Centre positions, not CIC support or macrocharge fraction; every particle is initially held.',
    }
    return {
        'path': candidate['path'], 'baseline': full['path'],
        'guard_cells': candidate['metrics']['guard_cells'],
        'first_center_numeric_difference': first_center,
        'last_equal_sampled_plane_step': last_equal_plane,
        'first_sampled_plane_numeric_difference': first_plane,
        'first_release': release, 'first_detector': detector,
        'first_particle_removal_frame': eviction,
        'first_ballistic_transition': {
            'observed_step': None,
            'earliest_possible_loop_step_from_detector_gate': detector['loop_step']+1,
            'reason': 'Only aggregate ballistic_outgoing_steps was logged. The branch requires detected=true before that step\'s detector event check; no exact first ballistic step can be recovered.',
        },
        'center_difference_precedes_any_release': first_center is not None and first_center['step'] < release['continuous_step'],
        'plane_difference_precedes_any_release': first_plane is not None and first_plane['step'] < release['continuous_step'],
        'center_difference_precedes_any_detector': first_center is not None and first_center['step'] < detector['continuous_step'],
        'plane_difference_precedes_any_detector': first_plane is not None and first_plane['step'] < detector['continuous_step'],
        'initial_cohort': cohort_match(full['initial'], candidate['initial']),
        'release_cohort': cohort_match(full['release'], candidate['release']),
        'all_admitted_particles_still_held_in_pre_release_frames': all(r['active'] == r['held'] for r in prior),
        'pre_release_max_continuity_relative': max(float(r['continuity_relative']) for r in prior),
        'pre_release_max_open_charge_balance_C': max(abs(float(r['charge_balance_C'])) for r in prior),
        'initial_region_particle_counts': initial_counts,
        'early_gauss_before_and_after_remap': gauss_rows,
        'field_by_time': by_time,
    }


def historical_comparison(old, new):
    common = ('beta_box', 'steps', 'dt_s', 'nx', 'ny', 'particles', 'seed', 'charge_C',
              'release_plane_lab_m', 'detector_plane_lab_m', 'window_beta')
    same = {k: old['metrics'][k] == new['metrics'][k] for k in common}
    require(all(same.values()), 'Historical full-control inputs differ')
    dz = new['metrics']['dt_s'] * C
    old_lo = float(old['history'][0]['origin_z_m'])
    new_lo = float(new['history'][0]['origin_z_m'])
    old_extent = [old_lo, old_lo + old['metrics']['window_nz']*dz]
    pml = new['metrics']['cpml_z_cells']
    new_extent = [new_lo+pml*dz, new_lo+(new['metrics']['window_nz']-pml)*dz]
    a, b, groups = field_groups(old['field'], new['field'])
    old_particles, new_particles = by_id(old['detector']), by_id(new['detector'])
    require(old_particles.keys() == new_particles.keys(), 'Historical detector cohort differs')
    ids = sorted(old_particles)
    x, y = [old_particles[i] for i in ids], [new_particles[i] for i in ids]
    total_charge = math.fsum(abs(float(r['charge_C'])) for r in x)
    weights = [abs(float(r['charge_C']))/total_charge for r in x]
    def mean_delta(rows):
        return math.fsum(w*float(r['delta_energy_lab_eV']) for w, r in zip(weights, rows))
    old_mean, new_mean = mean_delta(x), mean_delta(y)
    first_step = min(groups)
    return {
        'old_path': old['path'], 'new_path': new['path'],
        'identical_model_input_metrics': same,
        'old_PEC_domain_m': old_extent, 'new_nonabsorbing_domain_m': new_extent,
        'nonabsorbing_extent_max_absolute_difference_m': max(abs(x-y) for x, y in zip(old_extent, new_extent)),
        'old_boundary': '384-cell fixed full-domain PEC',
        'new_boundary': '424-cell fixed full-domain total-field CPML; 20 cells outside each original end',
        'initial_cohort': cohort_match(old['initial'], new['initial']),
        'release_cohort': cohort_match(old['release'], new['release']),
        'initial_solve_iterations': [old['metrics']['initial_iterations'], new['metrics']['initial_iterations']],
        'initial_gauss_relative': [old['metrics']['initial_gauss_relative'], new['metrics']['initial_gauss_relative']],
        'initial_detector_plane_E_V_m': difference(*groups[first_step], E),
        'initial_detector_plane_B_T': difference(*groups[first_step], B),
        'all_time_detector_plane_E_V_m': difference(a, b, E),
        'all_time_detector_plane_B_T': difference(a, b, B),
        'sampling_coordinates_max_abs_difference': {
            k: max(abs(float(v[k])-float(u[k])) for u, v in zip(a, b))
            for k in ('time_E_box_s', 'time_B_box_s', 'time_E_lab_s', 'x_m', 'y_m', 'z_box_m')},
        'detector_particle_differences_by_ID': {
            k: difference(x, y, (k,)) for k in ('x_m', 'y_m', 'z_box_m', 'ux', 'uy', 'uz',
                'time_box_s', 'time_lab_s', 'delta_energy_box_eV', 'delta_energy_lab_eV')},
        'mean_delta_energy_lab_eV_old': old_mean,
        'mean_delta_energy_lab_eV_new': new_mean,
        'mean_delta_energy_lab_eV_new_minus_old': new_mean-old_mean,
        'mean_paired_delta_energy_lab_difference_eV': math.fsum(
            w*(float(b['delta_energy_lab_eV'])-float(a['delta_energy_lab_eV']))
            for w, a, b in zip(weights, x, y)),
        'scope': 'Historical consistency supplement, not identical-boundary validation. Entire initial mesh was not exported; only initial detector-plane fields and the prescribed particle cohort can be compared directly.',
    }


def plot_onset(item, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    directory.mkdir(parents=True, exist_ok=True)
    observed = [row for row in item['field_by_time'] if row['step'] <= 48]
    fig, ax = plt.subplots(figsize=(9.0, 4.6), layout='constrained')
    for key, label, color, marker in [('E_V_m', 'E difference', '#20639b', 'o'),
                                       ('cB_V_m', 'cB difference', '#b2472c', 's')]:
        # Exact zeros are omitted from a logarithmic ordinate, never replaced
        # with a noise floor. Both series use only the every-four-step plane.
        values = [r[key]['vector_rms_absolute_difference'] for r in observed]
        ax.plot([r['step'] for r in observed], [v if v > 0.0 else math.nan for v in values],
                color=color, marker=marker, label=label, linewidth=1.6, markersize=4)
    center_step = item['first_center_numeric_difference']['step']
    release = item['first_release']['continuous_step']
    detector = item['first_detector']['continuous_step']
    ax.axvline(center_step, color='#555555', linestyle=':', linewidth=1.3,
               label='First central E difference: step %s' % center_step)
    ax.axvline(release, color='#6b8e23', linestyle='--', linewidth=1.3,
               label='First release: %.3f' % release)
    ax.axvline(detector, color='#8064a2', linestyle='-.', linewidth=1.1,
               label='First particle detection: %.3f' % detector)
    ax.set_yscale('log')
    ax.set_xlim(0, 49)
    ax.set_ylim(100, 1e6)
    ax.set_xlabel('Complete box-frame step n')
    ax.set_ylabel('Plane vector RMS difference (V/m)')
    ax.set_title('Static bunch, guard 24: field disagreement precedes release')
    ax.text(.025, .95, 'Plane samples at n = 0, 4, ..., 20:\nE and B differences are exactly zero.\nZeros omitted on the log scale.',
            transform=ax.transAxes, va='top', fontsize=9)
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, which='both', alpha=.18)
    fig.text(.01, -.04,
             'Native E(n), B(n-1/2); cB only converts units. Curves use plane data; the central sample is a timing marker.',
             fontsize=8)
    output = []
    for extension in ('png', 'pdf'):
        path = directory / ('causal_onset.' + extension)
        fig.savefig(path, dpi=180, bbox_inches='tight')
        output.append(str(path.resolve()))
    plt.close(fig)
    return output


def field_before_last_detection(full, candidate):
    dt = full['metrics']['dt_s']
    require(dt == candidate['metrics']['dt_s'], 'Different time steps')
    last_full = max(float(r['time_box_s'])/dt for r in full['detector'])
    last_candidate = max(float(r['time_box_s'])/dt for r in candidate['detector'])
    limit = min(last_full, last_candidate)
    a, b, _ = field_groups(full['field'], candidate['field'])
    pairs = [(x, y) for x, y in zip(a, b) if int(x['step']) <= limit]
    require(pairs, 'No fields before the last detection')
    x, y = [p[0] for p in pairs], [p[1] for p in pairs]
    return {
        'candidate_path': candidate['path'], 'full_path': full['path'],
        'mode': candidate['metrics']['mode'], 'beta_box': candidate['metrics']['beta_box'],
        'guard_cells': candidate['metrics']['guard_cells'],
        'last_full_detector_continuous_step_from_time': last_full,
        'last_candidate_detector_continuous_step_from_time': last_candidate,
        'upper_cutoff_step_minimum': limit,
        'included_common_field_steps': sorted({int(r['step']) for r in x}),
        'plane_pixel_records_per_run': len(x),
        'E_V_m': difference(x, y, E), 'cB_V_m': difference(x, y, B, C),
        'scope': 'Common recorded fields no later than either run\'s last particle detection. Earlier detected particles may already be ballistic or removed. This interval is not an outflow-free control, and its L2 must not be subtracted from another interval as a reflection coefficient.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', type=Path)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    root = workspace / 'results/streaming_cpml_20261027'
    inputs = Inputs()
    official = inputs.json(root / 'REVIEW.json')
    require(official['cohort_charge_and_input_consistency_status'] == 'PASS',
            'Official cohort/charge/input review not passed')
    require(official['coverage']['complete_planned_ten_case_set'], 'Ten-case set incomplete')
    for campaign in ('preflight-r1', 'main-r1'):
        inputs.record(root / campaign / 'COMPLETED.json')
        for name in ('streaming_cpml_bunch_probe.cpp', 'streaming_cpml_boundary.cpp',
                     'streaming_reference_boundary.cpp', 'moving_window_remap.cpp'):
            inputs.record(root / campaign / 'source/verification' / name)
    for name in ('eb_cpml.cpp', 'eb_initial_field.cpp'):
        inputs.record(workspace / 'source/aprl/src' / name)
    full = case(inputs, root / 'preflight-r1/b0-full-g24')
    paths = {16: 'main-r1/b0-cpml-g16', 24: 'preflight-r1/b0-cpml-g24',
             40: 'main-r1/b0-cpml-g40'}
    static = [chronology(full, case(inputs, root / paths[g])) for g in (16, 24, 40)]
    old_root = workspace / 'results/streaming_bunch_20261027'
    dynamic_full = case(inputs, root / 'main-r1/b023-full-g24')
    history = [historical_comparison(case(inputs, old_root / 'static-preflight-r1/b0-full-g16'), full),
               historical_comparison(case(inputs, old_root / 'guard-beta-r1/b023-full-g16'),
                                     dynamic_full)]
    before_last = []
    for entry in official['cases']:
        if entry['metrics']['mode'] != 'full':
            candidate = case(inputs, Path(entry['path']))
            baseline = full if candidate['metrics']['beta_box'] == 0.0 else dynamic_full
            before_last.append(field_before_last_detection(baseline, candidate))
    report = {
        'format': 'aprl-streaming-cpml-onset-diagnostic-v1',
        'method': 'Offline existing CSV analysis only; no new PIC or physics threshold.',
        'script_path': str(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__)),
        'input_sha256': inputs.hashes,
        'static_guard_chronology': static, 'historical_full_control_comparison': history,
        'field_until_common_last_particle_detection': before_last,
        'operation_contract': [
            'One full 424-cell Cowan initial solve, physical domain [-64,320] um, then a copied initial cache; the moving run never reads evolved full fields.',
            'At beta_box=0, integer reference queries return exact stored values; there is no subcell reference translation during this run.',
            'The short grid copies initial total E/B into its CPML as well as its interior. CPML histories begin at zero. Incoming held charge is retained through the original deposition.',
            'Each step: reference E halo(n), native Cowan magnetic update plus original total-field CPML correction, unchanged particle/deposition, native electric update plus original CPML and homogeneous termination, then existing reference override of exterior Ex/Ey and transverse PEC.',
            'After a complete step: WorldPosition field/history remap; newly entering E/B use cached native state n+1, newly unsupported histories use zero. No residual subtraction, projection or source correction.',
        ],
        'interpretation_limits': [
            'Earliest recorded detector disagreement is not the time or location where the full mesh first differs; the plane is sampled every four steps and central E every step.',
            'Disagreement before first detection rules out post-detector ballistic motion and removal as the initial cause, but does not bound their later contributions.',
            'For static beta0, reference subcell translation/CIC mismatch cannot explain the first discrepancy. This does not isolate initial total-field CPML response, moving-interface/history effects, and exterior reference termination from each other.',
            'The CPML corrects z derivative terms even when an unmodified static Cowan curl cancels between axes. With zero initial history and alpha=0, its first-step correction factor is decay/kappa-1. This is an explicit possible source of local response, not a measured attribution of the detector error.',
            'Before/after remap Gauss diagnostics cover different physical sets: a formerly absorbing/interface layer can enter the reported nonabsorbing set. A jump alone does not establish a wrong remap or lost particle charge.',
            'The aggregate history has no first ballistic event timestamp. A geometry estimate must not be labelled an observed event.',
            'cB differences are reported in V/m for scale; native E/B times are staggered. These diagnostics are not simultaneous lab energy flux or pure radiation.',
            'Historical old full PEC and new full CPML have different boundaries; small observed differences are supplementary consistency evidence only.',
        ],
    }
    if args.plots:
        report['plots'] = plot_onset(next(x for x in static if x['guard_cells'] == 24), args.plots)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(str(args.output.resolve()))
    for item in static:
        print('guard', item['guard_cells'], 'first center', item['first_center_numeric_difference']['step'],
              'first plane', item['first_sampled_plane_numeric_difference']['step'],
              'release', item['first_release']['continuous_step'],
              'detector', item['first_detector']['continuous_step'],
              'first removal', item['first_particle_removal_frame']['step'])


if __name__ == '__main__':
    main()
