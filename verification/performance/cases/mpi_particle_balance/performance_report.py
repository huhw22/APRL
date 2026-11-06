"""Read closed MPI-balancing campaigns and report measured throughput/load.

No HDF5 is opened. Particle counts are actual converged-call counts from the
complete rank CSVs; modeled balance scores are never treated as stage timing.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import sys

import numpy as np
import yaml


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def fields(line):
    return dict(re.findall(r'([A-Za-z0-9_]+)=([^\s]+)', line))


def number(value):
    result = float(value)
    require(math.isfinite(result), 'Nonfinite log number')
    return result


def stats(values):
    values = [float(value) for value in values if value is not None]
    if not values:
        return dict(count=0, median=None, minimum=None, maximum=None)
    return dict(count=len(values), median=float(np.median(values)),
                minimum=min(values), maximum=max(values))


def resource_times(text, expected_steps):
    rows = [fields(line) for line in text.splitlines() if line.startswith('[resource]')
            and fields(line).get('phase') in ['progress', 'final']]
    require(rows and rows[-1]['phase'] == 'final', 'Missing final resource report')
    require(int(rows[-1]['step']) == expected_steps, 'Resource final step differs')
    points = []
    for row in rows:
        item = dict(step=int(row['step']), elapsed_s=number(row['elapsed_s']),
                    reported_cumulative_mean_step_s=number(row['measured_step_s']),
                    phase=row['phase'])
        require(item['step'] > 0, 'Resource step must be positive')
        if points and item['step'] == points[-1]['step']:
            require(row['phase'] == 'final', 'Repeated nonfinal progress step')
            points[-1] = item
        else:
            require(not points or (item['step'] > points[-1]['step'] and
                    item['elapsed_s'] >= points[-1]['elapsed_s']), 'Resource time axis is not monotonic')
            points.append(item)
    previous_step, previous_time = 0, 0.
    windows = []
    for point in points:
        count = point['step'] - previous_step
        elapsed = point['elapsed_s'] - previous_time
        require(count > 0 and elapsed >= 0., 'Invalid elapsed-time interval')
        windows.append(dict(first_step=previous_step + 1, last_step=point['step'],
            steps=count, elapsed_s=elapsed, measured_window_mean_step_s=elapsed / count,
            cumulative_elapsed_s=point['elapsed_s'],
            reported_cumulative_mean_step_s=point['reported_cumulative_mean_step_s'],
            final_window=point['phase'] == 'final'))
        previous_step, previous_time = point['step'], point['elapsed_s']
    memory = {}
    for key, label in [('peak_rss_rank_max_mib', 'largest_rank_peak_RSS_MiB'),
                        ('peak_rss_total_mib', 'sum_of_rank_peak_RSS_MiB'),
                        ('current_rss_total_mib', 'largest_observed_simultaneous_total_RSS_MiB')]:
        present = [number(row[key]) for row in rows if key in row]
        memory[label] = max(present) if present else None
    memory['interpretation'] = ('Sum of rank-local peaks is not necessarily a simultaneous job RSS peak; '
                                 'simultaneous total is sampled only at resource-report times.')
    return windows, points[-1]['elapsed_s'], memory


def coupling_steps(directory, ranks, steps):
    paths = sorted((directory / 'reference').glob('source-coupling-rank-*.csv'))
    require([path.name for path in paths] ==
            [f'source-coupling-rank-{rank:05d}.csv' for rank in range(ranks)], 'Missing/extra rank histories')
    committed = np.zeros((steps, ranks), dtype=np.int64)
    trials = np.zeros_like(committed)
    iterations = np.zeros_like(committed)
    maximum_force = maximum_work = 0.
    for rank, path in enumerate(paths):
        require(path.with_suffix('.complete').read_text().split() == ['1', str(steps)],
                'Unclosed rank history: ' + str(path))
        data = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        require(np.array_equal(data['step'], np.arange(1, steps + 1)), 'Missing/repeated rank step')
        require(all(np.isfinite(data[name]).all() for name in data.dtype.names), 'Nonfinite rank history')
        for name, target in [('committed_pushes', committed), ('trial_pushes', trials),
                              ('iterations', iterations)]:
            require(np.all((data[name] >= 0) & (data[name] == np.floor(data[name]))),
                    'Noninteger/negative particle-call count')
            target[:, rank] = data[name].astype(np.int64)
        maximum_force = max(maximum_force, float(np.max(data['self_force_defect'])))
        maximum_work = max(maximum_work, float(np.max(data['relative_fixed_point_work_defect'])))
    require(np.all(iterations == iterations[:, :1]), 'Conservative trial iteration differs across ranks')
    total = committed.sum(axis=1)
    trial_total = trials.sum(axis=1)
    require(np.all(trial_total >= total), 'Trial calls fewer than converged committed calls')
    maximum = committed.max(axis=1)
    active = np.count_nonzero(committed, axis=1)
    mean = total / ranks
    rows = []
    for k in range(steps):
        rows.append(dict(step=k + 1, committed_total=int(total[k]), trial_total=int(trial_total[k]),
            ranks_with_committed_pushes=int(active[k]), maximum_rank_committed=int(maximum[k]),
            all_rank_mean_committed=float(mean[k]),
            actual_committed_max_over_mean=float(maximum[k] / mean[k]) if mean[k] else None,
            trial_over_committed=float(trial_total[k] / total[k]) if total[k] else None,
            conservative_iterations=int(iterations[k, 0])))
    summary = dict(committed_total=int(total.sum()), trial_total=int(trial_total.sum()),
        trial_over_committed=float(trial_total.sum() / total.sum()) if total.sum() else None,
        actual_committed_max_over_mean=stats(row['actual_committed_max_over_mean'] for row in rows),
        ranks_with_committed_pushes=stats(active), maximum_rank_committed=stats(maximum),
        maximum_conservative_iterations=int(iterations.max()),
        maximum_self_force_defect=maximum_force, maximum_fixed_point_work_defect=maximum_work,
        definition='Counts of full particle calls in the conservative-source method; internal Boris substeps are not counted here. A zero-committed rank may still perform fields, frozen-source work and communication.')
    return rows, summary


def balance_events(text, case, card):
    events = []
    for line in text.splitlines():
        if not line.startswith('[particle-balance] step='):
            continue
        row = fields(line)
        cuts = [int(value) for value in row['cuts'].split(':') if value]
        require(len(cuts) == case['ranks'] + 1 and cuts[0] == 0 and
                cuts[-1] == case['active_grid'][2] and all(a < b for a, b in zip(cuts, cuts[1:])),
                'Invalid logged integer partition')
        event = dict(step=int(row['step']), event=int(row['event']),
            logged_successful_call_s=number(row['repartition_s']),
            particle_packet_bytes=int(row['particle_packet_bytes']), cuts=cuts,
            modeled_max_before=number(row['modeled_max_before']),
            modeled_max_after=number(row['modeled_max_after']),
            modeled_imbalance_before=number(row['modeled_imbalance_before']),
            modeled_imbalance_after=number(row['modeled_imbalance_after']))
        require(event['logged_successful_call_s'] >= 0 and event['event'] == len(events) + 1,
                'Invalid repartition timing or event sequence')
        events.append(event)
    interval = card['mesh'].get('mpi_particle_balance', {}).get('interval_steps')
    if interval:
        require(all(event['step'] % int(interval) == 0 for event in events), 'Repartition occurs outside configured cadence')
    if case['mode'] == 'off':
        require(not events, 'Static control unexpectedly repartitioned')
    elif case['mode'] == 'on':
        require(events, 'Balanced control did not repartition')
    times = [event['logged_successful_call_s'] for event in events]
    summary = dict(count=len(events), successful_logged_seconds_sum=sum(times),
        successful_logged_seconds_mean=float(np.mean(times)) if times else 0.,
        successful_logged_seconds_range=[min(times), max(times)] if times else [0., 0.],
        configured_interval_steps=interval,
        interpretation='Log timers cover accepted balance calls, including their planning/remap work. Rejected/unchanged planner calls are not separately timed in this log. This is not all scheduling overhead or a disjoint stage percentage.',
        modeled_scores_are_timing=False)
    return events, summary


def profile_records(text):
    records = []
    for line in text.splitlines():
        if not line.startswith('[performance]'):
            continue
        row = fields(line)
        if 'stage' in row:
            records.append(dict(phase=row['phase'], step=int(row['step']), steps=int(row['steps']),
                elapsed_s=number(row['elapsed_s']), stage=row['stage'],
                mean_rank_s=number(row['mean_s']), max_rank_s=number(row['max_s']),
                maximum_rank=int(row['max_rank'])))
    return records


def read_case(case, plan):
    directory = Path(case['directory'])
    prefix = case.get('bounded_duration_prefix') is True
    marker = directory / ('PREFIX_COMPLETE.json' if prefix else 'PIC_COMPLETE.json')
    require(marker.is_file(), 'Refusing to read a running/incomplete case: ' + case['name'])
    closed = json.loads(marker.read_text())
    if prefix:
        require(closed.get('prefix_complete') is True and closed.get('physical_complete') is False,
                'Unsuccessful or promoted duration prefix')
    else:
        require(closed.get('complete') is True, 'Unsuccessful PIC completion')
    runtime = closed['runtime']
    require(runtime.get('exit_code') == (1 if prefix else 0) and number(runtime['seconds']) > 0,
            'Invalid closed launch timing/termination')
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    text = (directory / 'aprl.log').read_text()
    if prefix:
        # Standalone reports use the campaign's frozen close-proof reader and
        # its helpers, just as the original prefix runner did.
        scripts = Path(plan['directory']) / 'scripts'
        require(all((scripts / name).is_file() for name in
                    ['prefix_closed.py', 'common.py', 'duration_review.py', 'report_discrete.py']),
                'Missing frozen duration-prefix readers')
        if str(scripts) not in sys.path:
            sys.path.insert(0, str(scripts))
        from prefix_closed import duration_status
        duration_status(directory, case, runtime)
    else:
        require(all(value in text for value in ['Configured stop reached:', 'phase=solve-complete',
                                               'phase=mpi-finalize-complete']), 'Missing normal completion log')
    steps = int(case['expected_steps'])
    windows, loop_seconds, memory = resource_times(text, steps)
    rows, calls = coupling_steps(directory, int(case['ranks']), steps)
    events, balance = balance_events(text, case, card)
    result_path = directory / 'RESULT.json'
    science = json.loads(result_path.read_text()) if result_path.exists() else {}
    science_pass = bool(science.get('program_integrity_pass', science.get('passed', False)))
    profile_interval = int(runtime.get('environment', {}).get('APRL_PROFILE_INTERVAL',
                           case.get('profile_interval', plan.get('profile_interval', 0))))
    profile = profile_records(text)
    require(profile_interval != 0 or not profile, 'Profile0 case contains stage profiling records')
    summary = dict(name=case['name'], family=case['family'], mode=case['mode'],
        directory=str(directory), nodes=case['nodes'], ranks=case['ranks'], threads=case['threads'],
        steps=steps, whole_launch_seconds=number(runtime['seconds']),
        solver_loop_seconds=loop_seconds, full_loop_mean_step_seconds=loop_seconds / steps,
        launch_minus_loop_seconds=number(runtime['seconds']) - loop_seconds,
        profile_interval=profile_interval, memory=memory, calls=calls, balance=balance,
        particle_steps=rows, resource_windows=windows, balance_events=events,
        stage_profile=profile, scientific_result_present=result_path.exists(),
        scientific_program_integrity_pass=science_pass, scientific_result=str(result_path),
        completion_record=str(marker), bounded_duration_prefix=prefix,
        PIC_complete=None if prefix else str(marker), physical_acceptance=False)
    return summary


def comparisons(cases):
    grouped = {}
    for case in cases:
        require(case['mode'] in ['off', 'on'], 'Performance report requires static/off and balanced/on cases')
        group = grouped.setdefault(case['family'], {})
        require(case['mode'] not in group, 'Repeated mode within one pair family: ' + case['family'])
        group[case['mode']] = case
    pairs = []
    for family, values in grouped.items():
        require(set(values) == {'off', 'on'}, 'Missing static/balanced pair: ' + family)
        first, second = values['off'], values['on']
        require(first['steps'] == second['steps'] and first['nodes'] == second['nodes'] and
                first['ranks'] == second['ranks'] and first['threads'] == second['threads'],
                'Paired performance duration or resource layout differs')
        pairs.append(dict(family=family, static=first['name'], balanced=second['name'],
            whole_launch_speedup=first['whole_launch_seconds'] / second['whole_launch_seconds'],
            whole_launch_time_reduction_percent=100 * (1 - second['whole_launch_seconds'] / first['whole_launch_seconds']),
            solver_loop_speedup=first['solver_loop_seconds'] / second['solver_loop_seconds'],
            solver_loop_time_reduction_percent=100 * (1 - second['solver_loop_seconds'] / first['solver_loop_seconds'])))
    grouped_summary = {}
    for mode in ['off', 'on']:
        values = [case for case in cases if case['mode'] == mode]
        grouped_summary[mode] = {key:stats(case[key] for case in values) for key in
            ['whole_launch_seconds', 'solver_loop_seconds', 'full_loop_mean_step_seconds']}
        grouped_summary[mode]['successful_balance_calls'] = stats(case['balance']['count'] for case in values)
        grouped_summary[mode]['logged_successful_balance_call_seconds'] = stats(
            case['balance']['successful_logged_seconds_sum'] for case in values)
    pair_count = len(pairs)
    interpretation = (
        'One matched static/balanced pair. Median equals that single paired value and minimum=maximum. '
        'These are not repeat statistics or an estimate of run-to-run variability.' if pair_count == 1 else
        f'{pair_count} separately executed static/balanced pairs. Median and full observed range describe '
        'these repetitions; no confidence interval or statistical significance is claimed.')
    return dict(pairs=pairs, pair_count=pair_count, repeated_pairs=pair_count > 1,
        modes=grouped_summary,
        paired_whole_launch_speedup=stats(pair['whole_launch_speedup'] for pair in pairs),
        paired_whole_launch_time_reduction_percent=stats(pair['whole_launch_time_reduction_percent'] for pair in pairs),
        paired_solver_loop_speedup=stats(pair['solver_loop_speedup'] for pair in pairs),
        paired_solver_loop_time_reduction_percent=stats(pair['solver_loop_time_reduction_percent'] for pair in pairs),
        interpretation=interpretation)


def campaign_scope(plan):
    """Describe only the declared observation scope, never complete an FEL claim."""
    cases = plan.get('cases', [])
    declarations = [str(plan.get('scope', '')), str(plan.get('verification_scope', ''))]
    declarations.extend(str(case.get('verification_scope', '')) for case in cases)
    text = ' '.join(declarations).lower()
    prefix = ('prefix' in text or plan.get('require_all_particles_through_planes') is False or
              any(case.get('require_all_particles_through_planes') is False for case in cases))
    if prefix:
        steps = sorted({int(case['expected_steps']) for case in cases})
        label = '/'.join(str(value) for value in steps)
        return dict(kind='partial_prefix', declared_scope=declarations,
            title=f'Full-input closed{label}-step prefix',
            description=f'Closed{label}-step prefix of the original full input; partial detector cohorts. '
                        'This does not establish full-beam or complete-FEL completion.',
            plot_note='Partial detector cohorts; no claim of full-beam or complete-FEL completion',
            comparison_note='Both controls use the same prefix duration, input and observation geometry. '
                            'Performance does not establish full-beam or complete-FEL completion.')
    if 'long-window' in text or 'long window' in text:
        return dict(kind='short_source_long_window_stress', declared_scope=declarations,
            title='Short-source / long-window performance stress',
            description='The same short-source, long-window geometry is used by both controls. '
                        'This stress test is not the37-period physical FEL acceptance.',
            plot_note='Identical short-source / long-window stress geometry; modeled scores are not timings',
            comparison_note='No gain is credited to changing the window. This stress test does not '
                            'substitute for the37-period physical FEL test.')
    return dict(kind='paired_ownership', declared_scope=declarations,
        title='Paired MPI ownership performance',
        description='Measured performance of the declared matched static/balanced observation.',
        plot_note='Matched observation; modeled load scores are not timings',
        comparison_note='Performance does not establish independent physical acceptance.')


def write_csv(output, cases):
    rows, windows, case_rows = [], [], []
    for case in cases:
        intervals = case['resource_windows']
        cursor = 0
        for step in case['particle_steps']:
            while step['step'] > intervals[cursor]['last_step']:
                cursor += 1
            window = intervals[cursor]
            rows.append(dict(case=case['name'], family=case['family'], mode=case['mode'], **step,
                resource_window_first_step=window['first_step'],
                resource_window_last_step=window['last_step'],
                resource_window_mean_step_s=window['measured_window_mean_step_s']))
        windows.extend(dict(case=case['name'], family=case['family'], mode=case['mode'], **row)
                       for row in intervals)
        case_rows.append(dict(case=case['name'], family=case['family'], mode=case['mode'],
            whole_launch_s=case['whole_launch_seconds'], solver_loop_s=case['solver_loop_seconds'],
            loop_mean_step_s=case['full_loop_mean_step_seconds'],
            successful_balance_calls=case['balance']['count'],
            successful_balance_log_s=case['balance']['successful_logged_seconds_sum'],
            successful_balance_log_mean_s=case['balance']['successful_logged_seconds_mean'],
            largest_rank_peak_RSS_MiB=case['memory']['largest_rank_peak_RSS_MiB'],
            sum_of_rank_peak_RSS_MiB=case['memory']['sum_of_rank_peak_RSS_MiB']))
    for name, values in [('performance.csv', rows), ('time_windows.csv', windows),
                          ('case_summary.csv', case_rows)]:
        with (output / name).open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)


def plot(output, cases, comparison, scope):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-mpi-balance-performance-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10.5, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.2), layout='constrained')
    colors = {'off':'#4D4D4D', 'on':'#0072B2'}
    families = list(dict.fromkeys(case['family'] for case in cases))
    for case in cases:
        style = '-' if families.index(case['family']) % 2 == 0 else '--'
        color = colors[case['mode']]
        label = ('Static' if case['mode'] == 'off' else 'Balanced') + ' / ' + case['family']
        intervals = case['resource_windows']
        axes[0, 0].plot([row['last_step'] for row in intervals],
            [row['measured_window_mean_step_s'] for row in intervals], style + 'o',
            color=color, markersize=3.5, linewidth=1.5, label=label)
        steps = case['particle_steps']
        axes[0, 1].plot([row['step'] for row in steps],
            [row['actual_committed_max_over_mean'] if row['actual_committed_max_over_mean'] is not None
             else np.nan for row in steps], style, color=color, linewidth=1.1, label=label)
        axes[1, 0].plot([row['step'] for row in steps],
            [row['ranks_with_committed_pushes'] for row in steps], style,
            color=color, linewidth=1.3, label=label)
    for ax, title, ylabel in [(axes[0, 0], 'Time per log interval', 'Interval mean (s / field step)'),
                              (axes[0, 1], 'Actual committed particle calls', 'Maximum rank / all-rank mean'),
                              (axes[1, 0], 'Ranks with committed pushes', 'Number of ranks')]:
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel('Maxwell step')
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.25)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].axhline(1., color='.65', linewidth=.7)
    axes[1, 0].set_ylim(0, max(case['ranks'] for case in cases) * 1.05)
    ax = axes[1, 1]
    positions = np.arange(len(cases), dtype=float)
    launch = np.asarray([case['whole_launch_seconds'] for case in cases])
    loop = np.asarray([case['solver_loop_seconds'] for case in cases])
    ax.bar(positions - .17, launch, width=.32, color=[colors[case['mode']] for case in cases],
           alpha=.85, label='Whole launch')
    ax.bar(positions + .17, loop, width=.32, color='white', edgecolor=[colors[case['mode']] for case in cases],
           linewidth=1.5, label='Solver loop')
    ax.set_xticks(positions, [case['name'].replace('perf-', '') for case in cases], rotation=20)
    ax.set_ylabel('Measured wall time (s)')
    ax.set_title('Runtime and accepted balancing cost', fontsize=10.5)
    ax.grid(axis='y', alpha=.25)
    secondary = ax.twinx()
    secondary.spines['top'].set_visible(False)
    balance = [case['balance']['successful_logged_seconds_sum'] for case in cases]
    secondary.plot(positions, balance, 'D', color='#B35806', markersize=5,
                   label='Accepted balance calls only')
    secondary.set_ylim(bottom=0)
    secondary.set_ylabel('Logged balance-call time (s)', color='#8C4605')
    handles, labels = ax.get_legend_handles_labels()
    extra_handles, extra_labels = secondary.get_legend_handles_labels()
    ax.legend(handles + extra_handles, labels + extra_labels, fontsize=7.5, loc='upper left')
    ranks = '/'.join(str(value) for value in sorted({case['ranks'] for case in cases}))
    pair_count = comparison['pair_count']
    pairs = 'one matched pair' if pair_count == 1 else f'{pair_count} repeated pairs'
    fig.suptitle(f'{scope["title"]}: {pairs} ({ranks} MPI ranks, one thread each)\n'
                 + scope['plot_note'], fontsize=10.5)
    fig.savefig(output / 'performance.png', dpi=240, bbox_inches='tight', pad_inches=.12)
    fig.savefig(output / 'performance.pdf', bbox_inches='tight', pad_inches=.12)
    plt.close(fig)


def main(root, output=None, pilots=()):
    root = Path(root).resolve()
    output = Path(output).resolve() if output else root / 'performance-report'
    plan = json.loads((root / 'plan.json').read_text())
    require(output != root and not output.exists(), 'Use a new separate report directory')
    require(len(plan['cases']) in [2, 4], 'Expected two cases for one pair or four cases for two pairs')
    cases = [read_case(case, plan) for case in plan['cases']]
    comparison = comparisons(cases)
    require(comparison['pair_count'] == len(cases) // 2, 'Incomplete or duplicated performance pairs')
    scope = campaign_scope(plan)
    complete_path = root / ('PREFIX_COMPLETE.json' if plan.get('bounded_duration_prefix') is True else 'COMPLETE.json')
    complete = json.loads(complete_path.read_text()) if complete_path.exists() else {}
    science_pass = bool(complete.get('program_integrity_pass', complete.get('passed', False)))
    science_pass = science_pass and all(case['scientific_program_integrity_pass'] for case in cases)
    pilot_records = []
    for directory in pilots:
        directory = Path(directory).resolve()
        marker = json.loads((directory / 'PIC_COMPLETE.json').read_text())
        require(marker.get('complete') is True, 'Stage pilot is incomplete')
        pilot_records.append(dict(directory=str(directory),
            stages=profile_records((directory / 'aprl.log').read_text()),
            interpretation='Measured pilot stages can overlap; their maxima cannot be added or transferred as percentages of unprofiled main-run wall time.'))
    output.mkdir(parents=True)
    summary = dict(read_only_report_complete=True, physical_acceptance=False, campaign=str(root),
        cases=cases, pair_count=comparison['pair_count'], analysis_scope=scope,
        comparisons=comparison, scientific_comparison_complete=complete_path.exists(),
        scientific_comparison_program_integrity_pass=science_pass,
        scientific_comparison_record=str(complete_path), stage_pilots=pilot_records,
        limitations=['All data were read after verified PIC or explicit duration-prefix closure; no running HDF5 was opened.',
            'Whole launch is the timed PIC srun call, including startup and teardown but excluding binding checks, postprocessing, queue wait and build.',
            'Interval means are elapsed-time differences, not reported cumulative means or individual-step timers.',
            'Call imbalance is actual committed particle calls, not model score or CPU utilization.',
            'Logged accepted balance-call cost omits rejected/unchanged planner polling; no timing decomposition percentage is inferred.',
            'Profile0 establishes throughput only. Nested stage maxima are not additive if optional pilots exist.',
            comparison['interpretation'], scope['comparison_note']])
    write_csv(output, cases)
    plot(output, cases, comparison, scope)
    save(output / 'SUMMARY.json', summary)
    speed = comparison['paired_whole_launch_speedup']
    reduction = comparison['paired_whole_launch_time_reduction_percent']
    if comparison['pair_count'] == 1:
        statistics_line = (f'**One matched pair(2 launches)**. Whole-launch speedup: '
            f'**{speed["median"]:.4f}x**; time reduction: **{reduction["median"]:.2f}%**. '
            'The stored median equals this pair and min=max; no repeated-run variability is estimated.')
    else:
        statistics_line = (f'**{comparison["pair_count"]} matched pairs({len(cases)} launches)**. '
            f'Whole-launch speedup: median **{speed["median"]:.4f}x**, '
            f'observed range **{speed["minimum"]:.4f}–{speed["maximum"]:.4f}x**. '
            f'Time reduction: median **{reduction["median"]:.2f}%**, '
            f'range **{reduction["minimum"]:.2f}–{reduction["maximum"]:.2f}%**.')
    lines = ['# MPI balancing: closed-run performance report', '',
             scope['description'], '', statistics_line, '',
             'Scientific comparison status: **' + ('passed' if science_pass else 'not established or failed') + '**. '
             f'The independent result is [{complete_path.name}]({complete_path}); performance does not establish physical acceptance.', '',
             '![Measured throughput and actual particle load](performance.png)', '',
             '| Case | Whole launch (s) | Solver loop (s) | Loop mean (s/step) | Accepted balance calls | Logged accepted-call time (s) |',
             '|---|---:|---:|---:|---:|---:|']
    for case in cases:
        lines.append(f'| {case["name"]} | {case["whole_launch_seconds"]:.3f} | '
                     f'{case["solver_loop_seconds"]:.3f} | {case["full_loop_mean_step_seconds"]:.5f} | '
                     f'{case["balance"]["count"]} | {case["balance"]["successful_logged_seconds_sum"]:.3f} |')
    lines += ['', 'The first panel uses differences in resource elapsed time between log reports. '
        'The repeated window means in performance.csv are not individual-step measurements; '
        'time_windows.csv contains each independent time window once. The second and third panels '
        'use complete per-rank conservative-source CSV counts. Zero committed calls do not imply zero field or communication work.', '',
        'Whole launch is the timed PIC srun call, including startup and teardown; it excludes '
        'binding checks, postprocessing, queue wait and build. The fourth panel separates '
        'whole launch from solver loop. Its second axis shows only '
        'logged successful balance calls; these timers include accepted planning/remap work but '
        'omit unchanged/rejected planner polls. They are a subset of runtime and are not added '
        'to wall time. No model score is treated as a timing component or percentage.', '',
        'RSS summaries distinguish the sum of per-rank peaks from the maximum simultaneously '
        'observed total. Optional stage pilots remain separate; nested stage timers are not additive. '
        + comparison['interpretation'], '', scope['comparison_note'], '',
        'Source records:', '']
    for case in cases:
        completion_label = 'Closed prefix; physical complete=false' if case['bounded_duration_prefix'] else 'PIC completion'
        lines.append(f'- {case["name"]}: [{completion_label}]({case["completion_record"]}), '
                     f'[science result]({case["scientific_result"]}), '
                     f'[log]({Path(case["directory"]) / "aprl.log"}).')
    (output / 'README.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(dict(report=str(output), scientific_comparison_pass=science_pass,
                          pair_count=comparison['pair_count'], analysis_scope=scope['kind'],
                          whole_launch_speedup=speed), indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--pilot', type=Path, action='append', default=[])
    args = parser.parse_args()
    main(args.root, args.output, args.pilot)
