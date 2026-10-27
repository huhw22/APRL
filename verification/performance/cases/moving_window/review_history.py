#!/usr/bin/env python3
"""Independently review the two moving-CPML history conventions.

The output accepts provenance/state-copy gates, not boundary physics. Plane
time sums of squared E/cB are amplitude diagnostics, not energy flux or a
reflection coefficient. Full per-pixel wave histories are not available.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


POLICIES = ('world-position', 'boundary-relative-z')
WAVES = {'wave': (8, POLICIES[0]), 'wave20': (20, POLICIES[0]),
         'wave-boundary': (8, POLICIES[1]), 'wave-boundary20': (20, POLICIES[1])}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def near(actual, expected, label, abs_tol=2e-15):
    require(math.isclose(actual, expected, rel_tol=2e-12, abs_tol=abs_tol),
            f'Independent arithmetic mismatch {label}: {actual} != {expected}')


def audit_campaign(root):
    require(not (root / 'FAILED.json').exists(), f'Failure marker in {root}')
    plan, completed = read_json(root / 'plan.json'), read_json(root / 'COMPLETED.json')
    require(completed['state_and_cost_probe_gates'] == 'COMPLETED', 'Incomplete campaign')
    for path, expected in read_json(root / 'frozen.json').items():
        require(digest(Path(path)) == expected, 'Changed frozen input: ' + path)
    require(digest(Path(plan['binary'])) == plan['binary_sha256'], 'Executable differs from plan')
    require(plan['ranks'] == plan['threads'] == 1 and plan['case_timeout_s'] == 180,
            'Unexpected resource scope')
    require('[binding] result=PASS failures=0' in (root / 'binding.log').read_text(),
            'Physical-core binding failed')
    names = [case['name'] for case in plan['cases']]
    require(len(names) == len(set(names)), 'Duplicate planned cases')
    require([case['case'] for case in completed['cases']] == names, 'Completion does not match plan')
    for case in completed['cases']:
        directory = root / case['case']
        require(digest(directory / 'metrics.json') == case['metrics_sha256'], 'Changed completed metrics')
        launch = read_json(directory / 'launch.json')
        require(launch['exit_code'] == 0 and not launch.get('timed_out') and not launch.get('error'),
                f'Unsuccessful launch: {directory}')
        require(read_json(directory / 'metrics.json')['completed'] is True, 'Case did not complete')
    return plan


def audit_verification(directory):
    metrics = read_json(directory / 'metrics.json')
    require(metrics['bitwise_failures'] == 0 and metrics['checks'] > 0, 'State-copy gate failed')
    require(metrics['tagged_shift_cases'] == 252 and metrics['particle_schedules'] == 8
            and metrics['particles'] == 256, 'Incomplete dual-policy coverage')
    policies = metrics['history_policies']
    require(len(policies) == 2 and {p['policy'] for p in policies} == set(POLICIES), 'Missing policy')
    for policy in policies:
        require(policy['tagged_shift_cases'] == 126 and policy['particle_schedules'] == 4,
                'Incomplete per-policy coverage')
    for key in ('zero_shift_policy_identity', 'excited_z_history_policy_difference',
                'next_step_policy_field_difference'):
        require(metrics[key] is True, 'Policy-specific gate failed: ' + key)
    rows = read_csv(directory / 'checks.csv')
    expected = {'tagged_native_fields_and_all_history_supports', 'history_omission_negative_control'}
    expected.update(f'PIC_schedule_{n}' for n in range(4))
    require(len(rows) == 15, 'Incomplete labelled verification records')
    for policy in POLICIES:
        observed = [row['test'] for row in rows if row['history_policy'] == policy]
        require(len(observed) == 6 and set(observed) == expected, 'Missing labelled policy tests')
    shared = [row['test'] for row in rows if row['history_policy'] == 'both']
    require(set(shared) == {'vacuum_zero_shift_policy_identity',
                           'nonzero_z_history_policy_distinction_shift_-1',
                           'nonzero_z_history_policy_distinction_shift_1'}, 'Missing policy distinction tests')
    counts = [int(row['checks_after']) for row in rows]
    require(counts == sorted(counts) and counts[-1] <= metrics['checks'], 'Invalid cumulative check count')
    return metrics


def summarize_wave(rows, metrics, steps=176):
    """Recompute observable norms without using precomputed error metrics."""
    require(steps >= 100 and len(rows) == steps + 1
            and [int(row['step']) for row in rows] == list(range(steps + 1)),
            'Missing, repeated, or unordered wave frames')
    for row in rows:
        require(all(math.isfinite(float(value)) for value in row.values()), 'Nonfinite wave sample')
        for signal in ('reference', 'candidate', 'difference'):
            e, b, combined = (float(row[f'{group}_{signal}_rms_V_m'])
                              for group in ('E', 'cB', 'combined'))
            require(min(e, b, combined) >= 0, 'Negative RMS')
            near(combined * combined, e * e + b * b, 'combined plane norm')
    dt = metrics['dt_s']
    periods = {'all': rows, 'main': rows[:100], 'late': rows[100:]}
    result = {'intervals': {'all': [0, steps], 'main': [0, 99], 'late': [100, steps]},
              'groups': {}, 'csv_frames': len(rows),
              'square_sum_contract': 'Equal sample weights; sum of plane RMS squared, not physical energy or flux',
              'time_square_integral_contract': 'dt times the sample square sum; V^2 s / m^2; no temporal-stagger correction'}
    for group in ('E', 'cB', 'combined'):
        values = {}
        for interval, samples in periods.items():
            by_signal = {}
            for signal in ('reference', 'candidate', 'difference'):
                series = [float(row[f'{group}_{signal}_rms_V_m']) for row in samples]
                square_sum = math.fsum(v * v for v in series)
                peak = max(series)
                by_signal[signal] = {'sample_count': len(series), 'sum_rms_squared_V2_m2': square_sum,
                                     'dt_sum_rms_squared_V2_s_m2': dt * square_sum,
                                     'time_rms_V_m': math.sqrt(square_sum / len(series)),
                                     'peak_plane_rms_V_m': peak,
                                     'peak_step': int(samples[series.index(peak)]['step'])}
            values[interval] = by_signal
        reference_sum = values['all']['reference']['sum_rms_squared_V2_m2']
        reference_peak = values['all']['reference']['peak_plane_rms_V_m']
        require(reference_sum > 0 and reference_peak > 0, 'No measured reference pulse')
        values['relative_l2_error'] = math.sqrt(values['all']['difference']['sum_rms_squared_V2_m2'] / reference_sum)
        for interval in periods:
            for signal in ('reference', 'candidate', 'difference'):
                item = values[interval][signal]
                item['square_sum_over_all_reference_square_sum'] = item['sum_rms_squared_V2_m2'] / reference_sum
                item['l2_over_all_reference_l2'] = math.sqrt(item['square_sum_over_all_reference_square_sum'])
                item['peak_over_all_reference_peak'] = item['peak_plane_rms_V_m'] / reference_peak
                item['time_rms_over_all_reference_peak'] = item['time_rms_V_m'] / reference_peak
        result['groups'][group] = values
    for group, name in [('E', 'E_relative_l2_error'), ('cB', 'B_relative_l2_error'),
                        ('combined', 'combined_relative_l2_error')]:
        near(result['groups'][group]['relative_l2_error'], metrics[name], name)
    combined = result['groups']['combined']
    for signal in ('reference', 'candidate', 'difference'):
        for suffix, key in [('l2_over_full_reference_l2', 'l2_over_all_reference_l2'),
                            ('peak_over_reference_peak', 'peak_over_all_reference_peak')]:
            name = f'late_{signal}_{suffix}'
            near(combined['late'][signal][key], metrics[name], name)
    near(combined['all']['difference']['peak_over_all_reference_peak'],
         metrics['peak_plane_difference_over_reference_peak_plane'], 'all peak ratio')
    result['main_frames_all_measured_components_identical'] = all(
        float(row['combined_difference_rms_V_m']) == 0 for row in rows[:100])
    result['first_nonzero_difference_step'] = next(
        (int(row['step']) for row in rows if float(row['combined_difference_rms_V_m']) != 0), None)
    return result


def review_wave(directory, depth, policy):
    metrics = read_json(directory / 'metrics.json')
    require(metrics['mode'] == directory.name and metrics['history_policy'] == policy, 'Wave policy mismatch')
    require(metrics['finite'] and metrics['initial_overlap_bitwise_failures'] == 0
            and metrics['initial_overlap_values_checked'] > 0, 'Wave input/finite gate failed')
    require(metrics['candidate_cpml_nonzero_seen'] and metrics['reference_cpml_nonzero_seen'],
            'Wave did not exercise CPML')
    expected = {'steps': 176, 'plane_frames': 177, 'pixels_per_plane': 80, 'nx': 32, 'ny': 4,
                'full_nz': 512, 'window_nz': 120 + depth, 'full_origin_global_k': -192,
                'window_initial_origin_global_k': 0, 'window_final_origin_global_k': -141,
                'window_beta': .8, 'detector_beta': .8, 'late_start_step': 100, 'late_frames': 77,
                'cpml_cells': [6, 0, depth], 'cpml_polynomial_order': 3, 'cpml_kappa_max': 8,
                'cpml_alpha_fraction': 0, 'cpml_target_reflection': 1e-8, 'pulse_center_global_k': 70,
                'pulse_half_support_cells': 16, 'carrier_wavelength_cells': 8, 'shifts': 141}
    for key, value in expected.items():
        require(metrics[key] == value, f'Changed wave fixture: {key}')
    retained = metrics['history_retained_bytes']
    require((retained == 0 if policy == POLICIES[0] else retained > 0), 'History retention policy not exercised')
    rows = read_csv(directory / 'wave-plane.csv')
    for row in rows:
        n = int(row['step'])
        require(int(row['window_origin_global_k']) == -((4 * n + 4) // 5), 'Incorrect integer origin')
        near(float(row['time_E_box_s']), n * metrics['dt_s'], 'E time', abs_tol=1e-27)
        near(float(row['time_B_box_s']), (n - .5) * metrics['dt_s'], 'B time', abs_tol=1e-27)
        near(float(row['detector_z_box_m']), (100 - .8 * n) * metrics['dz_m'], 'Detector position')
    summary = summarize_wave(rows, metrics)
    summary.update(metrics=metrics, csv_sha256=digest(directory / 'wave-plane.csv'))
    reference_columns = [key for key in rows[0] if 'reference' in key or key in
                         ('step', 'time_E_box_s', 'time_B_box_s', 'detector_z_box_m')]
    reference = [tuple(row[key] for key in reference_columns) for row in rows]
    return summary, reference


def review_cost(directory, width, policy):
    metrics = read_json(directory / 'metrics.json')
    require(metrics['numerical_status'] == 'PASS' and metrics['history_policy'] == policy,
            'Cost policy/numerical gate failed')
    require(metrics['width'] == width and metrics['steps'] == 64 and metrics['repeats'] == 3
            and metrics['bitwise_preflight_steps'] == 8, 'Changed benchmark fixture')
    rows = read_csv(directory / 'timing.csv')
    require(len(rows) == 9, 'Incomplete benchmark repeats')
    numeric = ('repeat', 'width', 'steps', 'kernel_s', 'move_s', 'field_move_s', 'history_move_s',
               'incoming_s', 'total_s', 'setup_s', 'aux_bytes', 'field_bytes', 'cpml_bytes',
               'moved_payload_bytes', 'incoming_bytes', 'history_zero_bytes', 'cell_updates')
    for row in rows:
        for key in numeric:
            require(math.isfinite(float(row[key])) and float(row[key]) >= 0, 'Invalid numeric timing: ' + key)
        require(int(row['width']) == width and int(row['steps']) == 64, 'Timing dimensions differ')
        require(abs(float(row['total_s']) - float(row['kernel_s']) - float(row['move_s'])) < .002,
                'Timing components inconsistent')
    groups = {mode: [row for row in rows if row['mode'] == mode] for mode in ('fixed', 'compact', 'snapshot')}
    modes = {}
    for mode, group in groups.items():
        require(sorted(int(row['repeat']) for row in group) == [0, 1, 2], 'Missing mode/repeat')
        modes[mode] = {key: statistics.median(float(row[key]) for row in group) for key in numeric}
        for bound, fn in [('min', min), ('max', max)]:
            modes[mode]['total_s_' + bound] = fn(float(row['total_s']) for row in group)
    paired = []
    for repeat in range(3):
        fixed = next(row for row in groups['fixed'] if int(row['repeat']) == repeat)
        compact = next(row for row in groups['compact'] if int(row['repeat']) == repeat)
        paired.append(float(compact['total_s']) / float(fixed['total_s']) - 1)
    return {'case': directory.name, 'policy': policy, 'width': width, 'metrics': metrics, 'modes': modes,
            'timing_sha256': digest(directory / 'timing.csv'),
            'same_width_total_overhead_fraction': modes['compact']['total_s'] / modes['fixed']['total_s'] - 1,
            'same_width_native_kernel_difference_fraction': modes['compact']['kernel_s'] / modes['fixed']['kernel_s'] - 1,
            'movement_over_native_kernel_fraction': modes['compact']['move_s'] / modes['compact']['kernel_s'],
            'paired_total_overhead_fractions': paired}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--previous', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root, previous = args.root.resolve(), args.previous.resolve()
    plan, old_plan = audit_campaign(root), audit_campaign(previous)
    expected_cases = {'verify', *WAVES}
    expected_cases.update(f'bench{width}' for width in (64, 128, 256))
    expected_cases.update(f'bench-boundary{width}' for width in (64, 128, 256))
    require({case['name'] for case in plan['cases']} == expected_cases, 'Unexpected history-comparison cases')
    sources, old_sources = read_json(root / 'source-files.json'), read_json(previous / 'source-files.json')
    core = {p: sha for p, sha in sources.items() if '/src/' in p}
    old_core = {p: sha for p, sha in old_sources.items() if '/src/' in p}
    require(core and core == old_core, 'Production kernel sources changed')
    require(digest(root / 'source/CMakeLists.txt') == digest(previous / 'source/CMakeLists.txt'),
            'Top-level build configuration changed')
    result = {'format': 'aprl-moving-window-history-review-v1', 'root': str(root), 'previous': str(previous),
              'source_base': plan['source_base'], 'binary_sha256': plan['binary_sha256'],
              'previous_binary_sha256': old_plan['binary_sha256'], 'provenance_gates': 'PASS',
              'production_core_source_hashes_unchanged': True, 'state_migration': 'PASS',
              'physical_status': 'MEASURED; no boundary-physics acceptance threshold imposed',
              'verification': audit_verification(root / 'verify'), 'waves': {}, 'costs': [],
              'review_script_sha256': digest(Path(__file__).resolve())}
    common_reference = None
    for case, (depth, policy) in WAVES.items():
        summary, reference = review_wave(root / case, depth, policy)
        if common_reference is None:
            common_reference = reference
        require(reference == common_reference, 'Common large-domain reference changed')
        if policy == POLICIES[0]:
            require((root / case / 'wave-plane.csv').read_bytes() == (previous / case / 'wave-plane.csv').read_bytes(),
                    'Default world-position wave no longer reproduces previous CSV exactly')
            summary['previous_csv_byte_identical'] = True
        result['waves'][case] = summary
    result['common_reference_columns_text_identical'] = True
    result['comparisons'] = []
    for depth, old_case, new_case in [(8, 'wave', 'wave-boundary'), (20, 'wave20', 'wave-boundary20')]:
        shared_keys = ('steps', 'plane_frames', 'pixels_per_plane', 'nx', 'ny', 'full_nz', 'window_nz',
                       'dx_m', 'dy_m', 'dz_m', 'dt_s', 'full_origin_global_k',
                       'window_initial_origin_global_k', 'window_final_origin_global_k',
                       'window_beta', 'detector_beta', 'window_initial_cpml_inner_lower_global_k',
                       'window_initial_cpml_inner_upper_global_k', 'cpml_cells', 'cpml_polynomial_order',
                       'cpml_target_reflection', 'cpml_kappa_max', 'cpml_alpha_fraction',
                       'pulse_center_global_k', 'pulse_half_support_cells', 'carrier_wavelength_cells',
                       'field_time_contract', 'plane_contract', 'source_contract')
        for key in shared_keys:
            require(result['waves'][old_case]['metrics'][key] == result['waves'][new_case]['metrics'][key],
                    f'Policy pair changed more than history handling: {key}')
        a, b = (result['waves'][case]['groups']['combined'] for case in (old_case, new_case))
        result['comparisons'].append({'z_cpml_cells': depth,
            'boundary_over_world_combined_error_l2_ratio': b['relative_l2_error'] / a['relative_l2_error'],
            'boundary_over_world_late_difference_peak_ratio':
                b['late']['difference']['peak_plane_rms_V_m'] / a['late']['difference']['peak_plane_rms_V_m'],
            'boundary_over_world_late_difference_square_sum_ratio':
                b['late']['difference']['sum_rms_squared_V2_m2'] / a['late']['difference']['sum_rms_squared_V2_m2'],
            'interpretation': 'Ratios of separate errors against the same reference; not a norm of the direct difference between policies'})
    for policy, stem in [(POLICIES[0], 'bench'), (POLICIES[1], 'bench-boundary')]:
        for width in (64, 128, 256):
            result['costs'].append(review_cost(root / f'{stem}{width}', width, policy))
    result['limitations'] = [
        'Only four longitudinal CPML history groups differ; the original recurrence and production kernels remain unchanged.',
        'A smaller error in one finite vacuum pulse is not general CPML, charged-bunch, weak-radiation, or production acceptance.',
        'Plane RMS and center Ey do not reconstruct per-pixel mutual intensity or direct inter-policy field-error norms.',
        'The late interval is steps 100 through 176; squared-amplitude sums are finite-interval diagnostics, not energy or reflection coefficients.',
        'Native E and B retain different time levels; no simultaneous lab-frame fields or Poynting flux are claimed.',
        'Cost overhead compares each policy with its own same-width fixed-window run; separate policy campaigns are not a paired timing experiment.',
        'Synthetic new-layer supply and one CPU rank omit production reference supply, MPI, and full Simulation workflows.'
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
