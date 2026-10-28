#!/usr/bin/env python3
"""Independent, read-only review of the completed full-r2 result files.

Only PERFORMANCE_REVIEW.json/txt are written. No simulation or Slurm commands.
The driver science gates are consumed, not redefined here; field amplitudes are
streamed separately to distinguish nonzero replay from a trivial zero result.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import h5py
import numpy as np

BASE = None
ROOT = None
SOURCES = {}


def read(path):
    raw = path.read_bytes()
    SOURCES[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
    return raw.decode()


def load(path):
    return json.loads(read(path))


def records(text, prefix):
    answer = []
    for line in text.splitlines():
        if not line.startswith(prefix):
            continue
        row = dict(re.findall(r'(\w+)=([^ ]+)', line))
        for key, value in row.items():
            try:
                row[key] = float(value) if any(c in value for c in '.eE') else int(value)
            except ValueError:
                pass
        answer.append(row)
    return answer


def pair(a, b):
    return {'baseline': a, 'candidate': b, 'candidate_over_baseline': b / a,
            'reduction_fraction': 1 - b / a, 'speedup_factor': a / b}


def field_amplitude(path):
    before = path.stat()
    stats = {name: {'max_abs': 0., 'scaled_sumsq': 0., 'nonzero_components': 0}
             for name in ('electric_V_per_m', 'magnetic_T')}
    with h5py.File(path, 'r') as handle:
        dataset = handle['field_plane/fields']
        shape = list(dataset.shape)
        for index in range(len(dataset)):
            frame = dataset[index]
            for name, stat in stats.items():
                values = frame[name]
                assert np.isfinite(values).all()
                maximum = float(np.max(np.abs(values)))
                stat['nonzero_components'] += int(np.count_nonzero(values))
                if maximum:
                    scale = max(maximum, stat['max_abs'])
                    stat['scaled_sumsq'] = (stat['scaled_sumsq'] *
                        (stat['max_abs'] / scale) ** 2 +
                        float(np.sum((values / scale) ** 2)))
                    stat['max_abs'] = scale
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    for stat in stats.values():
        stat['sample_component_l2_norm'] = stat['max_abs'] * math.sqrt(stat.pop('scaled_sumsq'))
        stat['all_zero'] = stat['nonzero_components'] == 0
    return {'shape': shape, 'fields': stats, 'stat_unchanged': True,
            'size_bytes': after.st_size, 'mtime_ns': after.st_mtime_ns,
            'scope': 'Unweighted raw sample/component norm; not energy or spectral power.'}


def main():
    global BASE, ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True,
                        help='Completed paired campaign containing TIMING_REVIEW.json')
    parser.add_argument('--output-directory', type=Path,
                        help='Defaults to the campaign parent; writes only PERFORMANCE_REVIEW.json/txt')
    args = parser.parse_args()
    ROOT = args.root.resolve()
    BASE = args.output_directory.resolve() if args.output_directory else ROOT.parent
    assert BASE.is_dir(), 'Output directory must already exist'
    driver = load(ROOT / 'TIMING_REVIEW.json')
    plan = load(ROOT / 'plan.json')
    assert driver['pass'] and driver['simulation_complete'] is False
    assert driver['status'] == 'TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED'
    assert all(case['pass'] for case in driver['cases'].values())
    assert all(comparison['pass'] for comparison in driver['comparisons'].values())
    result = {
        'review_pass': True, 'status': driver['status'], 'simulation_complete': False,
        'scope': 'Independent post-completion performance arithmetic; driver output/science gates consumed unchanged.',
        'binary_sha256': driver['binary_sha256'], 'binding': driver['bindings'], 'cases': {},
        'limits': [
            'One paired early 100-step sample on node3, 16 MPI x 8 threads, all 41,732,775 source particles; no repeat uncertainty or full-run speed claim.',
            'Each case ended by expected bounded-duration exception after diagnostics closed; case exit 1 is not physical completion.',
            'Pure progress endpoints exclude final teardown; final resource elapsed includes final sampling, close, and boundary reports.',
            'Profile 20-step cases diagnose stages only. Throughput claims use separate 100-step cases with profiling disabled.',
            'Stage maxima can belong to different ranks and cannot be summed into a wall-time budget. Communication, moving_window_total, and pending_statistics overlap inclusive timers.',
            'The first pending_statistics interval also includes the step-0 energy sample before loopWallStart; its accumulator is not reset at loop entry. It is diagnostic, not a pure per-step isolated timer.',
            'Process current RSS sum is instantaneous across the resource reductions; sum of process high-water RSS is not a simultaneous node peak.',
            'Initialization includes input/scatter/reference/queue/statistics/diagnostic setup. Reference preparation is only a subset; no standalone queue-build timer exists.',
        ]}
    for item in plan['cases']['full']:
        name = item['name']
        directory = ROOT / 'full' / name
        text = read(directory / 'aprl.log')
        launch = load(directory / 'aprl.log.launch.json')
        preparation = load(directory / 'reference' / 'preparation.json')
        resources = [r for r in records(text, '[resource]') if r.get('phase') in ('progress', 'final')]
        profiles = records(text, '[performance]')
        assert bool(profiles) == bool(item['profile_interval'])
        assert launch['environment']['APRL_PROFILE_INTERVAL'] == str(item['profile_interval'])
        assert launch['exit_code'] == driver['cases'][name]['exit_code'] == 1
        final = [r for r in resources if r['phase'] == 'final'][-1]
        progress = {r['step']: r for r in resources if r['phase'] == 'progress'}
        case = {'steps': item['steps'], 'profile_enabled': bool(item['profile_interval']),
                'resources': resources, 'launcher_seconds': launch['seconds'],
                'reference_preparation_seconds': preparation['preparationSeconds'],
                'launcher_minus_final_loop_seconds_not_pure_initialization': launch['seconds'] - final['elapsed_s'],
                'field_amplitude': field_amplitude(directory / 'detectors' / 'exit-field.h5')}
        assert math.isclose(final['elapsed_s'] / item['steps'], driver['cases'][name]['stepping_average_s'], rel_tol=1e-14)
        if item['steps'] == 100:
            assert not profiles and 50 in progress and 100 in progress
            case['pure_progress_0_100_s_per_step'] = progress[100]['elapsed_s'] / 100
            case['pure_progress_50_100_s_per_step'] = (progress[100]['elapsed_s'] - progress[50]['elapsed_s']) / 50
            case['final_minus_progress100_seconds'] = final['elapsed_s'] - progress[100]['elapsed_s']
            case['driver_final_endpoint_50_100_s_per_step'] = (final['elapsed_s'] - progress[50]['elapsed_s']) / 50
            assert math.isclose(case['driver_final_endpoint_50_100_s_per_step'], driver['cases'][name]['steps_50_100_average_s'], rel_tol=1e-14)
        else:
            init = [r for r in profiles if r.get('phase') == 'initialization'][0]
            intervals = [r for r in profiles if 'stage' in r]
            case['initialization_seconds'] = init
            case['initialization_max_minus_reference_seconds_not_queue_only'] = init['max_s'] - preparation['preparationSeconds']
            case['workload'] = [r for r in profiles if r.get('phase') == 'workload']
            case['profile_stages'] = {}
            for stage in dict.fromkeys(r['stage'] for r in intervals):
                rows = [r for r in intervals if r['stage'] == stage]
                assert sum(r['steps'] for r in rows) == 20
                case['profile_stages'][stage] = {
                    'rank_mean_seconds_per_step': sum(r['mean_s'] for r in rows) / 20,
                    'intervals': [{**r, 'max_s_per_step': r['max_s'] / r['steps']} for r in rows]}
            unique = {r['step']: r for r in intervals}
            case['profile_interval_wall_seconds'] = sum(r['elapsed_s'] for r in unique.values())
        result['cases'][name] = case

    a, b = (result['cases']['throughput-' + x] for x in ('baseline', 'candidate'))
    result['throughput_without_profiling'] = {key: pair(a[key], b[key]) for key in (
        'pure_progress_0_100_s_per_step', 'pure_progress_50_100_s_per_step',
        'driver_final_endpoint_50_100_s_per_step', 'launcher_seconds')}
    def progress100(case):
        return next(r for r in case['resources'] if r['phase'] == 'progress' and r['step'] == 100)
    ra, rb = progress100(a), progress100(b)
    result['memory_step100_without_profiling_mib'] = {key: pair(ra[key], rb[key]) for key in (
        'current_rss_rank_max_mib', 'current_rss_total_mib',
        'peak_rss_rank_max_mib', 'peak_rss_total_mib')}
    pa, pb = (result['cases']['profile-' + x] for x in ('baseline', 'candidate'))
    result['profile_initialization'] = {
        'total_initialization_max_seconds': pair(pa['initialization_seconds']['max_s'], pb['initialization_seconds']['max_s']),
        'reference_preparation_seconds': pair(pa['reference_preparation_seconds'], pb['reference_preparation_seconds']),
        'nonreference_initialization_seconds_not_queue_only': pair(
            pa['initialization_max_minus_reference_seconds_not_queue_only'], pb['initialization_max_minus_reference_seconds_not_queue_only'])}
    result['profile_common_stage_rank_means_s_per_step'] = {name: pair(value['rank_mean_seconds_per_step'],
        pb['profile_stages'][name]['rank_mean_seconds_per_step']) for name, value in pa['profile_stages'].items()}
    result['all_early_detector_fields_zero'] = all(field['all_zero'] for case in result['cases'].values()
        for field in case['field_amplitude']['fields'].values())
    result['sources_sha256'] = SOURCES
    (BASE / 'PERFORMANCE_REVIEW.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')

    t = result['throughput_without_profiling']['pure_progress_50_100_s_per_step']
    lines = [
        f"Independent performance review: {ROOT}; binding={driver['bindings']}",
        'PASS: arithmetic, profiling mode, completed-case and unchanged-file checks. All four cases are TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED, not physically completed simulations.',
        '',
        f"Unprofiled pure progress 50→100: {t['baseline']:.10f} → {t['candidate']:.10f} s/step; {100*t['reduction_fraction']:.5f}% less time, {t['speedup_factor']:.5f}x throughput.",
        f"Pure progress 0→100: {a['pure_progress_0_100_s_per_step']:.10f} → {b['pure_progress_0_100_s_per_step']:.10f} s/step.",
        f"Driver final endpoint adds {a['final_minus_progress100_seconds']:.9f} / {b['final_minus_progress100_seconds']:.9f} s of final control/sampling/close/report work to baseline / candidate. Its 50→100 result is reproduced but not used as the pure-loop result.",
        f"Whole 100-step launcher: {a['launcher_seconds']:.6f} → {b['launcher_seconds']:.6f} s; includes startup, reference, queue construction, stepping, close, srun/timeout teardown.",
        '',
        'RSS at progress step 100 (per-process accounting, GiB):',
    ]
    for key, value in result['memory_step100_without_profiling_mib'].items():
        lines.append(f"  {key}: {value['baseline']/1024:.6f} → {value['candidate']/1024:.6f}; reduction {100*value['reduction_fraction']:.5f}%")
    lines += ['', 'Initialization in separate 20-step PROFILE runs:']
    for key, value in result['profile_initialization'].items():
        lines.append(f"  {key}: {value['baseline']:.9f} → {value['candidate']:.9f} s; delta {value['candidate']-value['baseline']:+.9f} s")
    lines += ['No queue-only timer exists. The nonreference difference is consistent with added one-time queue/statistics work but also contains all other initialization and shared-file/runtime variation.',
        '', '20-step PROFILE stage diagnosis (not the throughput measurement):',
        'stage | baseline rank-mean s/step | candidate rank-mean s/step | per-interval maximum s/step at max_rank (baseline / candidate)']
    for name in pb['profile_stages']:
        x, y = pa['profile_stages'].get(name), pb['profile_stages'][name]
        def maxtext(value):
            return ', '.join(f"{r['max_s_per_step']:.7g}@{r['max_rank']}" for r in value['intervals']) if value else 'unmeasured'
        mean_a = f"{x['rank_mean_seconds_per_step']:.9g}" if x else 'unmeasured'
        lines.append(f"{name} | {mean_a} | {y['rank_mean_seconds_per_step']:.9g} | {maxtext(x)} / {maxtext(y)}")
    lines += ['', 'Never sum stage maxima or overlapping totals. Communication timers include MPI waiting; a large decrease does not itself show fewer bytes or a faster network.',
        'field_communication is inside field_total; particle_communication overlaps particle_migration_total; moving_window_total and pending_statistics overlap diagnostics. Changing the rank0 pending workload can reduce waits in other ranks.',
        'The first pending_statistics interval includes the step-0 energy sample before the loop starts; the existing accumulator is not reset at loop entry. No isolated first-interval per-step claim is made.',
        '', 'Early detector nonzero check:']
    for name, case in result['cases'].items():
        lines.append(name + ': ' + json.dumps(case['field_amplitude']['fields'], sort_keys=True))
    lines += [('All early detector fields are zero; this is only a zero-output consistency check.'
               if result['all_early_detector_fields_zero'] else
               'Early detector fields are nonzero in every case; their amplitudes are very small. This short replay does not establish late radiation equivalence.'),
        'Use the separate science/full-lifecycle evidence for nontrivial entry, propagation, crossing, and exit behavior.',
        '', 'Scope limits: single early pair, same node/layout/card; no later-time, whole-run or asymptotic O(L) proof. Full pending records and one-time preparation remain. Profile and throughput runs are kept separate.',
        'Sources and exact parsed values: PERFORMANCE_REVIEW.json; replay script: verification/performance/cases/pending_admission_optimization/review_performance.py --root <completed-campaign>. No core changes or jobs performed by this review.']
    (BASE / 'PERFORMANCE_REVIEW.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:18]))


if __name__ == '__main__':
    main()
