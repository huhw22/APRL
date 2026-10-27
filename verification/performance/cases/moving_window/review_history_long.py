#!/usr/bin/env python3
"""Review longer observation of the same two CPML history conventions.

The candidate window is unchanged. The large reference domain and observation
duration grow to exclude reference-boundary contamination and observe delayed
tails. Prefix reproduction is checked explicitly, without loosening tolerances.
"""
import argparse
import json
import math
from pathlib import Path

from review_history import (POLICIES, audit_campaign, digest, near, read_csv,
                            read_json, require, summarize_wave)


CASES = {'wave-long': ('wave', 8, POLICIES[0]),
         'wave-boundary-long': ('wave-boundary', 8, POLICIES[1]),
         'wave-long20': ('wave20', 20, POLICIES[0]),
         'wave-boundary-long20': ('wave-boundary20', 20, POLICIES[1])}
SHARED = ('pixels_per_plane', 'nx', 'ny', 'window_nz', 'dx_m', 'dy_m', 'dz_m', 'dt_s',
          'window_initial_origin_global_k', 'window_beta', 'detector_beta',
          'window_initial_cpml_inner_lower_global_k', 'window_initial_cpml_inner_upper_global_k',
          'cpml_cells', 'cpml_polynomial_order', 'cpml_target_reflection', 'cpml_kappa_max',
          'cpml_alpha_fraction', 'pulse_center_global_k', 'pulse_half_support_cells',
          'carrier_wavelength_cells', 'field_time_contract', 'plane_contract', 'source_contract')


def require_identical_prefix(long_path, short_path):
    actual = long_path.read_bytes().splitlines(keepends=True)[:178]
    expected = short_path.read_bytes().splitlines(keepends=True)
    if actual == expected:
        return
    a, b = read_csv(long_path)[:177], read_csv(short_path)
    details = {'long': str(long_path), 'short': str(short_path),
               'first_differing_line_including_header': next(
                   (n for n, (x, y) in enumerate(zip(actual, expected)) if x != y), None)}
    if a and b and list(a[0]) == list(b[0]) and len(a) == len(b):
        details['max_absolute_difference_by_column'] = {
            key: max(abs(float(x[key]) - float(y[key])) for x, y in zip(a, b)) for key in a[0]}
    raise RuntimeError('Long observation does not reproduce short CSV prefix byte for byte: '
                       + json.dumps(details, sort_keys=True))


def endpoint_observations(rows):
    """Report truncation/delay evidence as measurements, not acceptance gates."""
    result = {}
    for group in ('E', 'cB', 'combined'):
        values = [float(row[f'{group}_difference_rms_V_m']) for row in rows]
        main_reference_peak = max(float(row[f'{group}_reference_rms_V_m']) for row in rows)
        late, added = values[100:], values[177:]
        short_peak, long_peak = max(values[:177]), max(values)
        total_square = math.fsum(value * value for value in values)
        result[group] = {
            'short_observation_peak_step': values[:177].index(short_peak),
            'long_observation_peak_step': values.index(long_peak),
            'short_observation_peak_V_m': short_peak,
            'long_observation_peak_V_m': long_peak,
            'long_peak_over_short_peak': long_peak / short_peak if short_peak else None,
            'long_peak_occurs_after_short_end': values.index(long_peak) > 176,
            'last_frame_is_peak': values[-1] == long_peak,
            'last_frame_exceeds_previous_frame': values[-1] > values[-2],
            'last_frame_rms_V_m': values[-1],
            'last_frame_over_long_peak': values[-1] / long_peak if long_peak else None,
            'last_16_frames_time_rms_V_m': math.sqrt(math.fsum(v * v for v in values[-16:]) / 16),
            'added_interval': [177, len(rows) - 1],
            'added_interval_peak_V_m': max(added),
            'added_interval_peak_over_reference_peak': max(added) / main_reference_peak,
            'added_interval_sum_rms_squared_V2_m2': math.fsum(v * v for v in added),
            'added_interval_fraction_of_all_difference_square_sum':
                math.fsum(v * v for v in added) / total_square if total_square else None,
            'late_sample_count': len(late),
            'interpretation': 'A later peak can invalidate a short-window suppression inference; end decay does not prove infinite-time stability'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--short-review', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    short_review = read_json(args.short_review)
    require(short_review['format'] == 'aprl-moving-window-history-review-v1', 'Wrong short review format')
    short_root = Path(short_review['root'])
    plan, short_plan = audit_campaign(root), audit_campaign(short_root)
    require(short_review['binary_sha256'] == short_plan['binary_sha256'], 'Short review binary mismatch')
    require({case['name'] for case in plan['cases']} == set(CASES), 'Unexpected extended-observation cases')
    core = lambda p: {key: sha for key, sha in read_json(p / 'source-files.json').items() if '/src/' in key}
    require(core(root) and core(root) == core(short_root), 'Production core changed for extended observation')
    for name in ('moving_window_remap.h', 'moving_window_remap.cpp'):
        require(digest(root / 'source/verification' / name) == digest(short_root / 'source/verification' / name),
                'History-policy implementation changed during extended observation')
    require(digest(root / 'source/CMakeLists.txt') == digest(short_root / 'source/CMakeLists.txt'),
            'Top-level build configuration changed')
    result = {'format': 'aprl-moving-window-history-long-review-v1', 'root': str(root),
              'short_root': str(short_root), 'short_review_path': str(args.short_review.resolve()),
              'short_review_sha256': digest(args.short_review), 'source_base': plan['source_base'],
              'binary_sha256': plan['binary_sha256'], 'short_binary_sha256': short_plan['binary_sha256'],
              'provenance_gates': 'PASS', 'original_kernels_and_history_policies_unchanged': True,
              'physical_status': 'MEASURED; extended finite observation, no physical acceptance threshold',
              'review_script_sha256': digest(Path(__file__).resolve()),
              'shared_review_script_sha256': digest(Path(__file__).with_name('review_history.py')),
              'waves': {}, 'comparisons': []}
    common_reference = None
    for case, (short_case, depth, policy) in CASES.items():
        directory = root / case
        metrics = read_json(directory / 'metrics.json')
        short_metrics = read_json(short_root / short_case / 'metrics.json')
        require(metrics['mode'] == case and metrics['history_policy'] == policy, 'Long-case policy mismatch')
        require(metrics['finite'] and metrics['initial_overlap_bitwise_failures'] == 0
                and metrics['initial_overlap_values_checked'] > 0, 'Extended-wave input/finite gate failed')
        require(metrics['candidate_cpml_nonzero_seen'] and metrics['reference_cpml_nonzero_seen'], 'History not exercised')
        expected = {'extended_observation': True, 'steps': 384, 'plane_frames': 385,
                    'full_nz': 1024, 'full_origin_global_k': -512,
                    'window_nz': 120 + depth, 'window_final_origin_global_k': -308, 'shifts': 308,
                    'late_start_step': 100, 'late_frames': 285, 'cpml_cells': [6, 0, depth]}
        for key, value in expected.items():
            require(metrics[key] == value, 'Unexpected extended-wave geometry: ' + key)
        for key in SHARED:
            require(metrics[key] == short_metrics[key], 'Candidate physical setup changed: ' + key)
        require(metrics['history_retained_bytes'] == 0 if policy == POLICIES[0]
                else metrics['history_retained_bytes'] > 0, 'Unexpected history retention')
        csv_path = directory / 'wave-plane.csv'
        short_csv = short_root / short_case / 'wave-plane.csv'
        require(digest(short_csv) == short_review['waves'][short_case]['csv_sha256'], 'Reviewed short wave changed')
        require_identical_prefix(csv_path, short_csv)
        rows = read_csv(csv_path)
        for row in rows:
            n = int(row['step'])
            require(int(row['window_origin_global_k']) == -((4 * n + 4) // 5), 'Wrong extended integer origin')
            near(float(row['time_E_box_s']), n * metrics['dt_s'], 'E time', abs_tol=1e-27)
            near(float(row['time_B_box_s']), (n - .5) * metrics['dt_s'], 'B time', abs_tol=1e-27)
            near(float(row['detector_z_box_m']), (100 - .8 * n) * metrics['dz_m'], 'Detector position')
        reference_keys = [key for key in rows[0] if 'reference' in key or key in
                          ('step', 'time_E_box_s', 'time_B_box_s', 'detector_z_box_m')]
        reference = [tuple(row[key] for key in reference_keys) for row in rows]
        if common_reference is None:
            common_reference = reference
        require(reference == common_reference, 'Extended common reference is not text-identical')
        summary = summarize_wave(rows, metrics, steps=384)
        summary.update(metrics=metrics, csv_sha256=digest(csv_path), short_case=short_case,
                       first_177_frames_csv_byte_identical=True, tail_observation=endpoint_observations(rows))
        result['waves'][case] = summary
    result['common_reference_columns_text_identical'] = True
    for depth, world, boundary, short_world, short_boundary in [
            (8, 'wave-long', 'wave-boundary-long', 'wave', 'wave-boundary'),
            (20, 'wave-long20', 'wave-boundary-long20', 'wave20', 'wave-boundary20')]:
        a, b = (result['waves'][case]['groups']['combined'] for case in (world, boundary))
        sa, sb = (short_review['waves'][case]['groups']['combined'] for case in (short_world, short_boundary))
        result['comparisons'].append({'z_cpml_cells': depth,
            'short_boundary_over_world_error_l2_ratio': sb['relative_l2_error'] / sa['relative_l2_error'],
            'long_boundary_over_world_error_l2_ratio': b['relative_l2_error'] / a['relative_l2_error'],
            'short_boundary_over_world_late_difference_peak_ratio':
                sb['late']['difference']['peak_plane_rms_V_m'] / sa['late']['difference']['peak_plane_rms_V_m'],
            'long_boundary_over_world_late_difference_peak_ratio':
                b['late']['difference']['peak_plane_rms_V_m'] / a['late']['difference']['peak_plane_rms_V_m'],
            'long_boundary_over_world_late_difference_square_sum_ratio':
                b['late']['difference']['sum_rms_squared_V2_m2'] / a['late']['difference']['sum_rms_squared_V2_m2'],
            'world_late_peak_step': a['late']['difference']['peak_step'],
            'boundary_late_peak_step': b['late']['difference']['peak_step'],
            'boundary_minus_world_peak_step': b['late']['difference']['peak_step'] - a['late']['difference']['peak_step'],
            'interpretation': 'Compare both amplitude and arrival time; a smaller short-interval error alone does not establish suppression'})
    result['limitations'] = [
        'Only the observation duration and large reference geometry change; the moving candidate window and history rules are unchanged.',
        'A finite 384-step pulse observation does not prove stability or complete integrated reflection over all time.',
        'CPML-depth memory may delay a response; extended peak and square-sum comparisons must accompany the original short-interval results.',
        'Reported squared-amplitude sums are not physical energy, energy flux, or isolated reflection coefficients.',
        'Native E/B are staggered in time, and only plane RMS plus center Ey are stored.',
        'This is not a new benchmark or production/MPI/charged-bunch acceptance.'
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
