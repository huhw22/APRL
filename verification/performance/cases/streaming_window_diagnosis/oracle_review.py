#!/usr/bin/env python3
"""Offline native-field oracle audit; no solver execution or physical pass gate."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import struct

import numpy as np

CASES = {
    'b023-full-record': 'full-record',
    'b023-residual-g24': 'residual',
    'b023-oracle-feed': 'residual-oracle-feed',
    'b023-oracle-cpml': 'residual-oracle-cpml',
    'b023-oracle-all': 'residual-oracle-all',
    'b023-native-b': 'residual-native-b',
}
LABELS = ('Frozen reference', 'Replay feed only', 'Replay CPML only',
          'Replay both', 'Native-B reconstruction')
CSV_NAMES = ('initial-particles.csv', 'release-particles.csv',
             'detector-crossings.csv', 'history.csv', 'field-plane.csv')
E_KEYS = ('Ex_V_m', 'Ey_V_m', 'Ez_V_m')
B_KEYS = ('Bx_T', 'By_T', 'Bz_T')
Z_STAGGER = (0., 0., .5, .5, .5, 0.)
C = np.longdouble('299792458')
REST = np.longdouble('9.1093837139e-31') * C*C / np.longdouble('1.602176634e-19')
BOOST = np.longdouble(5)/3
BOOST_BETA = np.longdouble('.8')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha256(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def table(path, sort_id=False):
    result = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
    require(result.size and result.dtype.names, 'Empty table: '+str(path))
    require(all(np.isfinite(result[k]).all() for k in result.dtype.names), 'Nonfinite table: '+str(path))
    if sort_id:
        result = np.sort(result, order='id')
        require(len(np.unique(result['id'])) == len(result), 'Duplicate particle ID: '+str(path))
    return result


def vector(table_, keys):
    return np.column_stack([table_[k] for k in keys])


def errors(reference, candidate):
    require(reference.shape == candidate.shape and reference.size, 'Mismatched/empty samples')
    r = np.asarray(reference, dtype=np.longdouble)
    d = np.asarray(candidate, dtype=np.longdouble)-r
    reference_l2 = np.sqrt(np.sum(r*r))
    delta_l2 = np.sqrt(np.sum(d*d))
    result = {'relative_l2': float(delta_l2/reference_l2) if reference_l2 else None,
              'difference_l2': float(delta_l2), 'reference_l2': float(reference_l2),
              'scalar_rms': float(np.sqrt(np.mean(d*d))),
              'maximum_absolute': float(np.max(np.abs(d))), 'scalar_values': int(d.size)}
    if d.ndim == 2 and d.shape[1] == 3:
        result['vector_rms'] = float(np.sqrt(np.mean(np.sum(d*d, axis=1))))
    return result


def particle_state(events, initial):
    require(np.array_equal(events['id'], initial['id']), 'Different initial/event cohorts')
    weights = np.abs(events['charge_C']).astype(np.longdouble)
    require(np.all(weights > 0), 'Zero particle weight')
    weights /= weights.sum()
    u, u0 = (vector(t, ('ux', 'uy', 'uz')).astype(np.longdouble) for t in (events, initial))
    gamma, gamma0 = (np.sqrt(1+np.sum(v*v, axis=1)) for v in (u, u0))
    dg = np.sum((u-u0)*(u+u0), axis=1)/(gamma+gamma0)
    delta = BOOST*(dg+BOOST_BETA*(u[:, 2]-u0[:, 2]))*REST
    box = dg*REST
    energy = BOOST*(gamma+BOOST_BETA*u[:, 2])*REST
    mean = np.sum(weights*energy)
    # This is an algebra/serialization check, not a physics tolerance.
    require(np.max(np.abs(delta-events['delta_energy_lab_eV'])) < 1e-9 and
            np.max(np.abs(box-events['delta_energy_box_eV'])) < 1e-9,
            'Momentum-derived energy disagrees with recorded energy')
    summary = {'count': len(events), 'mean_delta_energy_lab_eV': float(np.sum(weights*delta)),
               'mean_loss_lab_eV': -float(np.sum(weights*delta)),
               'sigma_energy_lab_eV': float(np.sqrt(np.sum(weights*(energy-mean)**2)))}
    return summary, u, gamma, weights


def plane_comparison(full, candidate, full_events, events):
    for key in ('step', 'pixel_i', 'pixel_j', 'time_E_box_s', 'time_B_box_s',
                'time_E_lab_s', 'x_m', 'y_m', 'z_box_m'):
        require(np.array_equal(full[key], candidate[key]), 'Different plane coordinates/times: '+key)
    first = min(np.min(full_events['step']+full_events['fraction']),
                np.min(events['step']+events['fraction']))
    last = min(np.max(full_events['step']+full_events['fraction']),
               np.max(events['step']+events['fraction']))
    masks = {'whole': np.ones(len(full), dtype=bool),
             'before_first_detection': full['step'] < first,
             'until_common_last_detection': full['step'] <= last,
             'step32': full['step'] == 32}
    result = {'first_detection_continuous_step': float(first),
              'common_last_detection_continuous_step': float(last), 'intervals': {}, 'by_step': []}
    for name, mask in masks.items():
        result['intervals'][name] = {'first_sample_step': int(full['step'][mask].min()),
                                    'last_sample_step': int(full['step'][mask].max())}
        for label, keys in (('E', E_KEYS), ('B', B_KEYS)):
            result['intervals'][name][label] = errors(vector(full[mask], keys), vector(candidate[mask], keys))
    for step in np.unique(full['step']):
        mask = full['step'] == step
        row = {'step': int(step)}
        for label, keys in (('E', E_KEYS), ('B', B_KEYS)):
            row[label] = errors(vector(full[mask], keys), vector(candidate[mask], keys))
        row['c_delta_B_vector_rms_V_m'] = float(C*row['B']['vector_rms'])
        result['by_step'].append(row)
    return result


def channel_checks(metrics):
    mode = metrics['mode']
    feed = mode in ('residual-oracle-feed', 'residual-oracle-all')
    cpml = mode in ('residual-oracle-cpml', 'residual-oracle-all')
    native_b = mode == 'residual-native-b'
    recording = mode == 'full-record'
    for key, expected in (('oracle_recording', recording), ('oracle_replay', feed or cpml),
                          ('oracle_feed_channel', feed), ('oracle_cpml_channel', cpml),
                          ('native_B_reference_reconstruction', native_b)):
        require(metrics[key] is expected, 'Incorrect channel flag: '+key)
    exterior = metrics['boundary_values']+metrics['halo_values']+metrics['new_layer_values']
    residual = metrics['residual_reference_input_values']
    expected_feed_calls = exterior if feed else (exterior+residual if native_b else 0)
    expected_cpml_calls = residual if cpml else 0
    require(metrics['oracle_feed_sample_calls'] == expected_feed_calls and
            metrics['oracle_cpml_sample_calls'] == expected_cpml_calls,
            'Reference channels do not have the isolated expected query counts')
    require(metrics['oracle_reader_queries'] == (exterior if feed else 0)+expected_cpml_calls,
            'Oracle reader query accounting differs from selected channels')
    if not (feed or cpml):
        require(all(metrics[k] == 0 for k in ('oracle_reader_frame_loads', 'oracle_reader_bytes',
                                              'oracle_reader_live_bytes', 'oracle_reader_queries')),
                'A non-replay case read full-domain history')
        require(recording or metrics['oracle_directory'] == '', 'Non-oracle case has an oracle directory')
    else:
        require(metrics['oracle_reader_frame_loads'] > 0 and metrics['oracle_reader_bytes'] > 0,
                'Replay did not read any recorded fields')
    require(metrics['oracle_recorded_frames'] == (177 if recording else 0), 'Wrong oracle recording length')
    return {k: metrics[k] for k in metrics if k.startswith('oracle_') or k == 'native_B_reference_reconstruction'}


def stage_comparison(full, candidate, width):
    require(set(full) == set(candidate), 'Different native section coverage')
    selected = []
    for key in ((0, 0),)+tuple((n, s) for n in (1, 2, 32, 40) for s in (1, 2, 3)):
        f, c = full[key], candidate[key]
        require(f['E_step'] == c['E_step'] and f['B_step'] == c['B_step'], 'Different native stage times')
        start = c['offset']-f['offset']
        groups = {name: {'E': [[], []], 'B': [[], []]} for name in
                  ('all', 'low_cpml', 'low_cpml_bulk', 'interior', 'high_cpml',
                   'high_cpml_bulk', 'first_native_plane')}
        for a, values in enumerate(c['fields']):
            reference = f['fields'][a][start:start+len(values)]
            require(start >= 0 and reference.shape == values.shape, 'Section world coordinates differ')
            z = np.arange(len(values))+Z_STAGGER[a]
            masks = {'all': np.ones(len(values), dtype=bool), 'low_cpml': z < 20,
                     'low_cpml_bulk': (z > 0) & (z < 20),
                     'interior': (z >= 20) & (z <= width-20), 'high_cpml': z > width-20,
                     'high_cpml_bulk': (z > width-20) & (z < width),
                     'first_native_plane': np.arange(len(values)) == 0}
            for region, mask in masks.items():
                group = groups[region]['E' if a < 3 else 'B']
                group[0].append(reference[mask].ravel());group[1].append(values[mask].ravel())
        selected.append({'state': key[0], 'stage': key[1], 'offset': c['offset'],
                         'regions': {region: {label: errors(np.concatenate(pair[0]), np.concatenate(pair[1]))
                                              for label, pair in group.items()} for region, group in groups.items()}})
    return selected


def audit_recording(directory, full_sections):
    manifest = json.loads((directory/'manifest.json').read_text())
    require(manifest['completed'] and manifest['frames'] == 177, 'Incomplete oracle recording')
    files = [directory/('state-%06d.bin' % n) for n in range(177)]
    require(sum(p.stat().st_size for p in files) == manifest['frame_bytes'], 'Oracle frame bytes mismatch')
    inspected = (0, 1, 32, 40)
    for state in inspected:
        reference = full_sections[(state, 3) if state else (0, 0)]
        with files[state].open('rb') as stream:
            require(stream.read(8) == b'APRLOF01', 'Wrong oracle frame magic')
            n, first, cells, count = struct.unpack('<4Q', stream.read(32))
            require(n == state and count == 6 and cells == 108, 'Wrong oracle frame geometry/state')
            for a in range(6):
                component, nx, ny, nz = struct.unpack('<4Q', stream.read(32))
                require(component == a, 'Wrong oracle component order')
                data = np.frombuffer(stream.read(nx*ny*nz*8), '<f8').reshape(nz, ny, nx)
                expected = reference['fields'][a][first:first+nz]
                require(data[:, 14, :].shape == expected.shape and
                        np.array_equal(data[:, 14, :].view(np.uint64), expected.view(np.uint64)),
                        'Oracle native values differ from independently written full sections')
            require(not stream.read(1), 'Trailing oracle frame bytes')
    return {'manifest': manifest, 'independently_compared_full_section_states': list(inspected),
            'comparison': 'All six components on recorded j=14 XZ sections, exact binary64 bits'}, files


def plots(directory, reports, full_moments):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = ('#626870', '#d08024', '#2777a6', '#8f57a2', '#20936d')
    styles = ('-', '--', '-.', ':', '-')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
    for row, label, color, style in zip(reports, LABELS, colors, styles):
        steps = np.array([p['step'] for p in row['fields']['by_step']])
        values = np.array([p['B']['relative_l2'] for p in row['fields']['by_step']])*100
        valid = values > 0
        for ax in axes:
            ax.plot(steps[valid], values[valid], linestyle=style, color=color, label=label, linewidth=1.9)
    first = min(r['fields']['first_detection_continuous_step'] for r in reports)
    last = min(r['fields']['common_last_detection_continuous_step'] for r in reports)
    for ax, upper, title in zip(axes, (176, 48), ('Complete observation', 'Early observation')):
        ax.set_yscale('log');ax.set_xlim(0, upper);ax.set_title(title)
        ax.axvline(32, color='#89929a', linestyle=':', linewidth=.9)
        ax.axvspan(first, min(last, upper), color='#a8b8c5', alpha=.13)
        ax.set_xlabel('Completed field state n');ax.grid(True, which='major', alpha=.18)
    axes[0].set_ylabel('Native B vector L2 difference / full B (%)')
    fig.suptitle('Detector-plane magnetic-field discrepancy', y=.99, fontsize=14)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.5, .93), ncol=3, frameon=False)
    fig.text(.04, .02, 'E(n), B(n−½), sampled at the detector position at n. Shading: first through common last detection.\n'
             'Exact zero errors are omitted on the log scale. Replay requires a completed full-domain run.', fontsize=9)
    fig.subplots_adjust(top=.76, bottom=.18, wspace=.09)
    for ext in ('png', 'pdf'):fig.savefig(directory/('native_B_temporal_error.'+ext), dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.3), sharey=True)
    keys = ('mean_delta_difference_eV', 'sigma_difference_eV', 'particle_energy_rms_error_eV')
    titles = ('Mean ΔE difference', 'Energy-spread difference', 'Per-particle energy RMS error')
    for ax, key, title in zip(axes, keys, titles):
        values = [r['particles'][key] for r in reports]
        for i, (value, color) in enumerate(zip(values, colors)):
            ax.scatter(value, i, color=color, s=58, zorder=3)
        if key != keys[-1]:ax.axvline(0, color='#89929a', linewidth=.8)
        elif all(v > 0 for v in values):ax.set_xscale('log')
        ax.set_xlabel('eV');ax.set_title(title, fontsize=11);ax.grid(True, axis='x', alpha=.18)
    axes[0].set_yticks(range(len(LABELS)));axes[0].set_yticklabels(LABELS);axes[0].invert_yaxis()
    fig.suptitle('Particle energy at the common detector, versus full domain', fontsize=14)
    fig.text(.03, .03, 'Same 1,024 particle IDs; charge-weighted moments and matched-particle errors. '
             'Full mean ΔE = %.7g eV; full σE = %.7g eV.' %
             (full_moments['mean_delta_energy_lab_eV'], full_moments['sigma_energy_lab_eV']), fontsize=9)
    fig.subplots_adjust(left=.20, right=.98, bottom=.17, top=.83, wspace=.20)
    for ext in ('png', 'pdf'):fig.savefig(directory/('particle_energy_comparison.'+ext), dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args();root = args.root.resolve()
    require((root/'COMPLETED.json').is_file(), 'Wait for the complete bounded campaign')
    output = root/'oracle-diagnosis';output.mkdir(exist_ok=True)
    inputs = {Path(__file__).resolve(), root/'COMPLETED.json', root/'plan.json'}
    parser_path = root/'analysis/analyze.py'
    spec = importlib.util.spec_from_file_location('frozen_section_parser', parser_path)
    helper = importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);inputs.add(parser_path)
    cases, channels = {}, {}
    for name, mode in CASES.items():
        directory = root/name;metrics = json.loads((directory/'metrics.json').read_text())
        require(metrics['completed'] and metrics['mode'] == mode and metrics['beta_box'] == .23 and
                metrics['guard_cells'] == 24 and metrics['steps'] == 176 and
                metrics['cpml_cells'] == [0, 0, 20] and metrics['history_policy'] == 'world-position',
                'Unexpected physical fixture: '+name)
        channels[name] = channel_checks(metrics)
        section_metrics = dict(metrics)
        if mode == 'full-record':section_metrics['mode'] = 'full'
        sections, manifest = helper.sections(directory, section_metrics)
        require(manifest['maximum_state_step'] == 40, 'Expected 40 native section steps')
        initial, release, events = (table(directory/file, True) for file in CSV_NAMES[:3])
        plane = table(directory/'field-plane.csv')
        require(len(initial) == len(events) == len(release) == 1024 and len(plane) == 18000,
                'Incomplete cohort or plane coverage')
        require(np.array_equal(np.unique(plane['step']), np.arange(0, 177, 4)), 'Different plane cadence')
        for step in np.unique(plane['step']):
            sample = plane[plane['step'] == step]
            require(len(sample) == 400 and len(set(zip(sample['pixel_i'], sample['pixel_j']))) == 400,
                    'Repeated or missing plane pixel')
        require(np.allclose(plane['time_E_box_s'], plane['step']*metrics['dt_s'], rtol=1e-14, atol=1e-27)
                and np.allclose(plane['time_B_box_s'], (plane['step']-.5)*metrics['dt_s'], rtol=1e-14, atol=1e-27),
                'Changed native field time layers')
        summary, u, gamma, weights = particle_state(events, initial)
        cases[name] = dict(metrics=metrics, sections=sections, initial=initial, events=events,
                           release=release, plane=plane, summary=summary, u=u, gamma=gamma, weights=weights)
        inputs.update(directory/file for file in CSV_NAMES)
        inputs.update(directory/file for file in ('metrics.json', 'field-sections.bin', 'field-sections.json'))
    full = cases['b023-full-record'];reports = []
    for name in list(CASES)[1:]:
        case = cases[name]
        require(np.array_equal(case['initial'], full['initial']) and
                np.array_equal(case['release'], full['release']), 'Prescribed initial/release states changed: '+name)
        require(np.array_equal(case['events']['id'], full['events']['id']) and
                np.array_equal(case['events']['charge_C'], full['events']['charge_C']), 'Different detector cohort/weights')
        d = case['u']-full['u']
        dg = np.sum(d*(case['u']+full['u']), axis=1)/(case['gamma']+full['gamma'])
        energy_error = BOOST*(dg+BOOST_BETA*d[:, 2])*REST
        particle = {**case['summary'],
                    'mean_delta_difference_eV': case['summary']['mean_delta_energy_lab_eV']-full['summary']['mean_delta_energy_lab_eV'],
                    'sigma_difference_eV': case['summary']['sigma_energy_lab_eV']-full['summary']['sigma_energy_lab_eV'],
                    'particle_energy_rms_error_eV': float(np.sqrt(np.sum(case['weights']*energy_error**2))),
                    'particle_energy_maximum_abs_error_eV': float(np.max(np.abs(energy_error)))}
        reports.append({'case': name, 'fields': plane_comparison(full['plane'], case['plane'], full['events'], case['events']),
                        'particles': particle, 'native_stages': stage_comparison(full['sections'], case['sections'], 104)})
    plan = json.loads((root/'plan.json').read_text())
    legacy = []
    pairs = [('b023-residual-g24', root.parent/'r1/b023-residual-g24'),
             ('b023-full-record', root.parent/'r1/b023-full-g24'),
             ('b023-full-record', Path(plan['legacy_replay']['b023-full-record']['directory']))]
    for name, old in pairs:
        record = {'case': name, 'previous_directory': str(old), 'csv_sha256': {}}
        for file in CSV_NAMES:
            now, before = root/name/file, old/file
            require(now.read_bytes() == before.read_bytes(), 'Scientific CSV changed: '+str(now)+' vs '+str(before))
            inputs.add(before);record['csv_sha256'][file] = sha256(now)
        record['all_five_csv_byte_identical'] = True;legacy.append(record)
    oracle = root/'b023-full-record/native-oracle'
    recording, frames = audit_recording(oracle, full['sections']);inputs.update(frames)
    inputs.update(oracle/file for file in ('header.bin', 'manifest.json', 'complete.bin'))
    plots(output, reports, full['summary'])
    result = {'format': 'aprl-streaming-window-oracle-review-v1', 'campaign': str(root),
              'data_contract_status': 'PASS', 'physical_acceptance': 'Not assigned; no physical error threshold imposed',
              'binary_sha256': plan['binary_sha256'], 'channel_accounting': channels,
              'legacy_replay': legacy, 'recording_audit': recording, 'full_particle_moments': full['summary'],
              'cases': reports, 'input_sha256': {str(p): sha256(p) for p in sorted(inputs)},
              'limitations': ['Replay modes require completed full-domain history; they are causal diagnostics, not a production reference method or speedup.',
                              'Native-B reconstruction uses only the initial cache but remains a measured reference approximation.',
                              'E and B retain distinct native times; no mixed-time Lorentz transform or radiation-energy inference.',
                              'Before detection excludes every detector crossing; until common last detection may contain earlier outgoing particles.',
                              'Section norms combine native staggered component samples at j=14; plane norms use collocated sampled vectors.',
                              'No relative-error subtraction is interpreted as a separately measured reflection.']}
    (output/'REVIEW.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(output/'REVIEW.json')
    for row in reports:
        print(row['case'], 'B32=', row['fields']['intervals']['step32']['B']['relative_l2'],
              'Bwhole=', row['fields']['intervals']['whole']['B']['relative_l2'],
              'particle_RMS_eV=', row['particles']['particle_energy_rms_error_eV'])


if __name__ == '__main__':
    main()
