#!/usr/bin/env python3
"""Plot completed energy-transport acceptance data; never launch simulations.

Requires passing SCIENCE_REVIEW.json and PERFORMANCE_REVIEW.json before reading
HDF5. Creates a fresh output directory; does not edit any frozen campaign file.
PNG files require human visual inspection before publication/acceptance.
"""
import os
for variable in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import h5py
import numpy as np
_mpl_cache = tempfile.TemporaryDirectory(prefix='aprl-audit-matplotlib-')
os.environ.setdefault('MPLCONFIGDIR', _mpl_cache.name)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD = np.longdouble
MAGNETS = ('undulator', 'free-drift')
LABELS = {'undulator': 'One-period undulator', 'free-drift': 'Free drift'}
COLORS = {'old': '#657786', 'legacy': '#c97b29', 'boris': '#71429b',
          'particle': '#147e9b', 'field': '#cc503e', 'magnetic': '#4c965b'}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def source_case(plan, stage, magnet, variant):
    selected = [item for item in plan['cases'][stage]
                if item['magnet'] == magnet and item['variant'] == variant and item['repeat'] == 1]
    require(len(selected) == 1, 'Missing/unexpected case: ' + magnet + '/' + variant)
    return Path(selected[0].get('reused_from',selected[0]['directory']))


def read_records(path, group, sources):
    with h5py.File(path, 'r') as stream:
        g = stream[group]
        require(int(g['complete'][()]) == 1, 'Cannot plot unfinished HDF5: ' + str(path))
        committed = int(g['committed_records'][()])
        require(committed == len(g['records']) and committed > 1, 'Incomplete records: ' + str(path))
        rows = g['records'][:]
    sources[str(path)] = {'sha256': digest(path), 'records': len(rows), 'group': group}
    return rows


def scalar(rows, field):
    require(field in rows.dtype.names, 'Missing plot field: ' + field)
    value = rows[field].astype(LD)
    require(np.isfinite(value).all(), 'Nonfinite plotted field: ' + field)
    return value


def precise_terminal(row, field):
    require(field in row.dtype.names and field + '_low' in row.dtype.names,
            'Missing compensated particle energy: ' + field)
    value = LD(row[field]) + LD(row[field + '_low'])
    require(np.isfinite(value), 'Nonfinite compensated particle total')
    return value


def export_csv(path, columns, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def symlog_for_curves(axis, curves):
    peaks = [float(np.max(np.abs(curve))) for curve in curves if np.any(curve)]
    if not peaks:
        return
    threshold = max(min(peaks) / 10, max(peaks) * 1.e-14, 1.e-15)
    axis.set_yscale('symlog', linthresh=threshold, linscale=.6)
    axis.axhline(0, color='#aaaaaa', linewidth=.7, zorder=0)
    axis.grid(True, which='major', color='#dddddd', linewidth=.55)


def save_figure(figure, output, stem, outputs):
    # The default symlog locator can stack zero and its adjacent labels.
    # Retain its numeric ticks, but remove labels closer than24 display pixels.
    figure.canvas.draw()
    for axis in figure.axes:
        if axis.get_yscale() != 'symlog':
            continue
        lower, upper = axis.get_ylim()
        candidates = [value for value in axis.get_yticks() if lower <= value <= upper]
        zero_y = axis.transData.transform((0, 0))[1]
        candidates.sort(key=lambda value: abs(axis.transData.transform((0, value))[1]-zero_y))
        selected, positions = [], []
        for value in candidates:
            position = axis.transData.transform((0, value))[1]
            if all(abs(position-other)>=24 for other in positions):
                selected.append(value); positions.append(position)
        axis.set_yticks(sorted(selected))
    for extension in ('png', 'pdf'):
        path = output / (stem + '.' + extension)
        figure.savefig(path, dpi=180, bbox_inches='tight', facecolor='white')
        outputs.append(str(path))
    plt.close(figure)


def boost_figure(science_plan, output, sources, outputs):
    figure, axes = plt.subplots(2, 2, figsize=(14, 9), sharex='col')
    export_rows, metrics = [], {}
    for index, magnet in enumerate(MAGNETS):
        old = read_records(source_case(science_plan, 'science', magnet, 'old-on') / 'ledger/energy-ledger.h5',
                           'energy_ledger', sources)
        new = read_records(source_case(science_plan, 'science', magnet, 'new-on') / 'ledger/energy-ledger.h5',
                           'energy_ledger', sources)
        require(np.array_equal(old['step'], new['step']), 'Ledger step vectors differ')
        step = new['step']
        old_r = scalar(old, 'balance_residual_J')
        legacy = scalar(new, 'balance_residual_J')
        boris = scalar(new, 'balance_residual_boris_J')
        particle = scalar(new, 'particle_work_balance_residual_J')
        field = scalar(new, 'field_exchange_residual_J')
        magnetic = scalar(new, 'magnetic_rotation_delta_K_J')
        left, right = axes[index]
        curves = (old_r, legacy, boris)
        for values, color, style, label in zip(curves, (COLORS['old'], COLORS['legacy'], COLORS['boris']),
                ('-', '--', '-'), ('Accepted v4 residual', 'Candidate legacy residual', 'Candidate Boris residual')):
            left.plot(step, np.asarray(values * LD(1e15), dtype=float), color=color, linestyle=style,
                      linewidth=1.5 if style != '--' else 1.1, label=label)
        symlog_for_curves(left, [curve * LD(1e15) for curve in curves])
        split = (boris, particle, field, magnetic)
        for values, key, label in zip(split, ('boris', 'particle', 'field', 'magnetic'),
                ('Boris residual', 'Particle work residual', 'Field exchange residual', 'Magnetic rotation delta K')):
            right.plot(step, np.asarray(values * LD(1e15), dtype=float), color=COLORS[key],
                       linewidth=1.4 if key == 'boris' else 1.05, label=label,
                       linestyle='-' if key != 'magnetic' else ':')
        symlog_for_curves(right, [curve * LD(1e15) for curve in split])
        left.set_title(LABELS[magnet] + ' | complete residual history')
        right.set_title(LABELS[magnet] + ' | independent work / field terms')
        for axis in (left, right):
            axis.set_ylabel('Boost-frame energy [fJ], symmetric log scale')
            axis.legend(loc='upper left', fontsize=8, framealpha=.92)
            axis.set_xlim(int(step[0]), int(step[-1]))
        delta_k = scalar(new, 'delta_particle_kinetic_box_J')[-1]
        peak = int(np.argmax(abs(boris)))
        ratio = float(boris[-1] / abs(delta_k)) if delta_k else None
        metrics[magnet] = {'final_delta_K_box_J': float(delta_k),
            'old_final_residual_J': float(old_r[-1]), 'old_max_abs_residual_J': float(max(abs(old_r))),
            'boris_final_residual_J': float(boris[-1]), 'boris_max_abs_residual_J': float(abs(boris[peak])),
            'boris_peak_step': int(step[peak]), 'boris_final_over_abs_delta_K': ratio}
        annotation = ('Boris final: ' + format(float(boris[-1] * LD(1e15)), '.4g') + ' fJ; peak: ' +
                      format(float(abs(boris[peak]) * LD(1e15)), '.4g') + ' fJ at step ' + str(int(step[peak])) +
                      '\nFinal / |net particle change|: ' + ('undefined' if ratio is None else format(ratio, '.3%')))
        right.text(.985, .025, annotation, transform=right.transAxes, ha='right', va='bottom', fontsize=8,
                   bbox={'facecolor': 'white', 'edgecolor': '#cccccc', 'alpha': .94})
        for row_index, current in enumerate(step):
            export_rows.append([magnet, int(current)] + [str(values[row_index]) for values in
                               (old_r, legacy, boris, particle, field, magnetic)])
    for axis in axes[-1]:
        axis.set_xlabel('Maxwell step')
    figure.suptitle('Boost-frame accounting: full history and independently measured terms', fontsize=15)
    figure.text(.5, .015,
        'Residual identities are diagnostic checks, not proof of physical closure. Source supply and moving boundaries remain in the balance.\n'
        'No normalization by the large initial beam energy; signed values and transient maxima are retained.',
        ha='center', fontsize=9)
    figure.tight_layout(rect=(0, .065, 1, .96), h_pad=2)
    save_figure(figure, output, 'boost_energy_accounting', outputs)
    export_csv(output / 'boost_energy_series.csv',
               ['case', 'step', 'accepted_legacy_R_J', 'candidate_legacy_R_J', 'Boris_R_J',
                'particle_work_R_J', 'field_exchange_R_J', 'magnetic_rotation_delta_K_J'], export_rows)
    return metrics


def table(axis, columns, rows, widths):
    axis.axis('off')
    display = axis.table(cellText=rows, colLabels=columns, cellLoc='center', loc='center', colWidths=widths)
    display.auto_set_font_size(False)
    display.set_fontsize(10)
    display.scale(1, 2.4)
    for (row, column), cell in display.get_celld().items():
        cell.set_edgecolor('#c8cdd2')
        if row == 0:
            cell.set_facecolor('#e9eff4'); cell.get_text().set_weight('bold')
            cell.set_height(cell.get_height()*1.85)
        else:
            cell.set_facecolor('#ffffff' if row % 2 else '#f6f8fa')
    return display


def transport_figure(science_plan, science_report, output, sources, outputs):
    figure, axes = plt.subplots(2, 1, figsize=(16, 7))
    particle_rows, field_rows, csv_rows, coverage, metrics = [], [], [], [], {}
    fmt = lambda value: format(float(value), '.7g')
    for magnet in MAGNETS:
        case = source_case(science_plan, 'science', magnet, 'new-audit-on')
        record = {}
        for side in ('entry', 'exit'):
            record[side] = read_records(case / 'detectors' / science_plan['energy_plane_files'][side], 'energy_plane', sources)
        begin, end = record['entry'][-1], record['exit'][-1]
        cohort = [(int(row['particle_downstream_count']), int(row['particle_upstream_count'])) for row in (begin, end)]
        require(cohort[0] == cohort[1] == (science_plan['particles'], 0), 'Entry/exit full particle counts differ')
        k_in, k_out = [precise_terminal(row, 'particle_kinetic_lab_J') for row in (begin, end)]
        loss = k_in - k_out
        flux_in, flux_out = [LD(row['field_energy_lab_J']) for row in (begin, end)]
        flux_net = flux_in - flux_out
        physical_flux = LD(begin['physical_field_energy_lab_J']) - LD(end['physical_field_energy_lab_J'])
        transport_difference = loss + flux_net
        particle_rows.append([LABELS[magnet], str(cohort[0][0]) + ' / ' + str(cohort[1][0]),
            fmt(k_in * LD(1e6)), fmt(k_out * LD(1e6)), fmt(loss * LD(1e12))])
        field_rows.append([LABELS[magnet], fmt(flux_in * LD(1e12)), fmt(flux_out * LD(1e12)),
            fmt(flux_net * LD(1e12)), fmt(physical_flux * LD(1e12)), fmt(transport_difference * LD(1e12))])
        for side, rows in record.items():
            final = rows[-1]
            require(int(final['field_gap_count']) == 0, 'Cannot present uninterrupted field coverage with a gap')
            coverage.append(LABELS[magnet] + ' ' + side + ': ' + str(int(final['field_sample_count'])) +
                ' field samples, ' + str(int(final['outside_field_window_steps'])) + ' outside-window steps')
        csv_rows.append([magnet] + [str(value) for value in
            (k_in, k_out, loss, flux_in, flux_out, flux_net, physical_flux, transport_difference)])
        metrics[magnet] = {'particle_entry_kinetic_lab_J': str(k_in), 'particle_exit_kinetic_lab_J': str(k_out),
            'particle_loss_lab_J': str(loss), 'full_field_entry_minus_exit_J': str(flux_net),
            'physical_field_entry_minus_exit_J': str(physical_flux),
            'full_particle_plus_field_entry_minus_exit_J': str(transport_difference),
            'scope': 'Transport comparison; not a closed control-volume balance.'}
    axes[0].set_title('Particle crossing energies | full cohort | lab frame', fontsize=13, pad=10)
    table(axes[0], ['Case', 'Entry / exit count', 'Entry K [microJ]', 'Exit K [microJ]', 'Particle loss [pJ]\nentry minus exit'],
          particle_rows, [.24, .18, .18, .18, .22])
    axes[1].set_title('Signed field transport | positive along +z | lab frame', fontsize=13, pad=10)
    table(axes[1], ['Case', 'Full entry flux\nintegral [pJ]', 'Full exit flux\nintegral [pJ]',
          'Full field\nentry - exit [pJ]', 'Inner-aperture field\nentry - exit [pJ]',
          'Particle loss + full\nfield difference [pJ]'], field_rows, [.20, .15, .15, .16, .17, .17])
    figure.suptitle('Fixed-lab entry / exit transport comparison — not a closed energy balance', fontsize=15)
    footer = ('Small particle differences are formed from high + low components in long double before display; large totals are rounded here.\n'
        'Field flux is integrated only over available samples. Full particles include all crossings; full fields cover the finite mesh aperture.\n'
        'The last column excludes lateral flux, stored-field changes and source work. Total grid fields have no near-field subtraction.\n' +
        '\n'.join(coverage))
    figure.text(.5, .02, footer, ha='center', va='bottom', fontsize=8.7, linespacing=1.5)
    figure.subplots_adjust(top=.86, bottom=.34, hspace=.55)
    save_figure(figure, output, 'lab_plane_transport', outputs)
    export_csv(output / 'lab_transport_values.csv', ['case', 'K_entry_lab_J', 'K_exit_lab_J', 'particle_loss_lab_J',
        'field_entry_lab_J', 'field_exit_lab_J', 'full_field_entry_minus_exit_lab_J',
        'physical_field_entry_minus_exit_lab_J', 'particle_plus_full_field_entry_minus_exit_lab_J'], csv_rows)
    return metrics


def performance_figure(report, output, outputs):
    figure, axes = plt.subplots(2, 2, figsize=(15, 10))
    variants = ('old-on', 'new-on', 'new-off', 'new-audit-on')
    variant_labels = ('Accepted v4\nledger on', 'Candidate\nledger on', 'Candidate\ndiagnostics off', 'Candidate\nledger + planes')
    pairs = ('increment_on_new_vs_old', 'new_ledger_on_vs_off', 'new_planes_on_vs_no_planes', 'new_all_diagnostics_on_vs_off')
    pair_labels = ('Ledger on:\nnew vs old', 'Candidate:\nledger on / off', 'Candidate:\nplanes added', 'Candidate:\nall on / off')
    colors = ('#147e9b', '#cc503e', '#72954f')
    raw_csv, paired_csv = [], []
    for column, magnet in enumerate(MAGNETS):
        values = report['performance'][magnet]
        raw = np.array([values['configurations'][v]['samples_s_per_step'] for v in variants], dtype=float)
        require(raw.shape == (4, 3) and np.isfinite(raw).all() and np.all(raw>0), 'Expected all three timing repetitions')
        x = np.arange(4)
        for repeat in range(3):
            axes[0, column].plot(x, raw[:, repeat] * 1000, 'o-', color=colors[repeat],
                alpha=.85, linewidth=1, label='Repeat ' + str(repeat+1))
            for index, value in enumerate(raw[:, repeat]):
                raw_csv.append([magnet, repeat+1, variants[index], value])
        axes[0, column].set_xticks(x, variant_labels)
        axes[0, column].set_ylabel('Time per step [ms], steps 40–240')
        axes[0, column].set_title(LABELS[magnet] + ' | every throughput sample')
        axes[0, column].grid(axis='y', color='#dddddd', linewidth=.6)
        axes[0, column].legend(fontsize=9)
        pair_values = []
        for comparison, name in enumerate(pairs):
            samples = values['comparisons'][name]['paired_percent_changes']
            require(len(samples) == 3 and np.isfinite(samples).all(), 'Missing paired measurements')
            pair_values.append(samples)
            for repeat, value in enumerate(samples):
                paired_csv.append([magnet, name, repeat+1, value])
        paired_rows = []
        for repeat in range(3):
            paired_rows.append(['Repeat '+str(repeat+1)] + [format(v[repeat], '+.2f')+'%' for v in pair_values])
        paired_rows.append(['Paired median'] + [format(float(np.median(v)), '+.2f')+'%' for v in pair_values])
        paired_rows.append(['Ratio of medians'] + [format(values['comparisons'][name]['ratio_of_medians_percent'], '+.2f')+'%' for name in pairs])
        summary = table(axes[1, column], ['Time change'] + list(pair_labels), paired_rows, [.24,.19,.19,.19,.19])
        summary.set_fontsize(8.7)
        for (row, col), cell in summary.get_celld().items():
            if row in (1,2,3):
                cell.get_text().set_color(colors[row-1])
            if row in (4,5):
                cell.get_text().set_weight('bold')
        axes[1, column].set_title('Every paired change [%]; positive costs more')
    figure.suptitle('Small-case diagnostic overhead | fixed 16 MPI × 8 threads | every repetition retained', fontsize=15)
    figure.text(.5, .015, 'Same-node paired cases, profiling disabled; timing excludes preparation and closing.\n'
        'Three repeats characterize this short fixture only. No production-cost extrapolation or strict zero-cost claim.',
        ha='center', fontsize=9)
    figure.tight_layout(rect=(0, .07, 1, .955), h_pad=3)
    save_figure(figure, output, 'diagnostic_performance', outputs)
    export_csv(output / 'performance_samples.csv', ['case', 'repeat', 'variant', 'seconds_per_step'], raw_csv)
    export_csv(output / 'performance_pairs.csv', ['case', 'comparison', 'repeat', 'percent_change'], paired_csv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='Completed scientific campaign root')
    parser.add_argument('--performance-root', type=Path, help='Completed performance campaign; defaults to --root')
    parser.add_argument('--output', type=Path, required=True, help='New directory for figures and inspectable values')
    args = parser.parse_args()
    root = args.root.resolve()
    perf_root = (args.performance_root or root).resolve()
    output = args.output.resolve()
    require(not output.exists(), 'Use a new figure directory; do not overwrite reviewed artifacts')
    science = read_json(root / 'SCIENCE_REVIEW.json')
    performance = read_json(perf_root / 'PERFORMANCE_REVIEW.json')
    require(science['checks_pass'] and performance['checks_pass'], 'Both completed scientific and performance gates must pass')
    require(science['binary_sha256'] == performance['binary_sha256'], 'Science and performance used different binaries')
    plan = read_json(root / 'plan.json')
    require(plan['frozen'] and plan['binary_sha256'] == science['binary_sha256'], 'Scientific plan identity differs')
    perf_plan = read_json(perf_root / 'plan.json')
    require(perf_plan['frozen'] and perf_plan['binary_sha256'] == performance['binary_sha256'], 'Performance plan identity differs')
    output.mkdir(parents=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.titlepad': 11, 'pdf.fonttype': 42})
    sources, outputs = {}, []
    try:
        boost = boost_figure(plan, output, sources, outputs)
        transport = transport_figure(plan, science, output, sources, outputs)
        performance_figure(performance, output, outputs)
        manifest = {'format': 'aprl-energy-audit-figures-v1', 'script': str(Path(__file__).resolve()),
            'script_sha256': digest(__file__), 'science_root': str(root), 'performance_root': str(perf_root),
            'binary_sha256': science['binary_sha256'], 'input_hdf5': sources,
            'science_review_sha256': digest(root / 'SCIENCE_REVIEW.json'),
            'performance_review_sha256': digest(perf_root / 'PERFORMANCE_REVIEW.json'),
            'output_files': {path: digest(path) for path in outputs}, 'boost_metrics': boost,
            'lab_transport_metrics': transport, 'visual_inspection_complete': False,
            'required_next_step': 'Open all three PNG files, check clipping/overlap/legibility and source consistency; record visual review separately.',
            'scope': 'Diagnostic accounting and small-case overhead only; no claim of exact physical closure or isolated radiation.'}
        save_json(output / 'PLOT_MANIFEST.json', manifest)
        (output / 'VISUAL_REVIEW_REQUIRED.txt').write_text('All three PNG files require visual review. PDF files share the same matplotlib source but are not separately rendered.\n')
        print(output)
    except BaseException as error:
        save_json(output / 'INCOMPLETE.json', {'complete': False, 'error': str(error)})
        raise


if __name__ == '__main__':
    main()
