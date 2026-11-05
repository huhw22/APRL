"""Summarize closed-run field-origin evidence; no PIC or source modification."""
import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np


def main(root):
    root = Path(root).resolve()
    evidence = root / 'independent_evidence'
    data = {}
    for case in ['energy', 'momentum']:
        with np.load(evidence / (case + '-strict-traces.npz')) as f:
            data[case] = {k: f[k].copy() for k in f.files}
        assert all(np.isfinite(v).all() for v in data[case].values())
    channels = {}
    for case, d in data.items():
        channels[case] = {
            'native_change_J': d['native_field_J'] - d['native_field_J'][0],
            'numerical_time_cowan_J': sum(d[k] for k in [
                'time_electric_J', 'time_magnetic_J',
                'current_time_pairing_J', 'cowan_operator_J']),
            'current_input_J': -d['applied_current_mid_J'],
            'curl_remap_J': -d['discrete_curl_exchange_J'] + d['remap_J'],
            'inlet_cpml_J': sum(d[k] for k in [
                'inlet_overwrite_J', 'magnetic_cpml_J', 'electric_cpml_J']),
            'operator_residual_J': d['field_operator_residual_J']}
        x = channels[case]
        identity = x['native_change_J'] - sum(v for k, v in x.items() if k != 'native_change_J')
        assert np.max(np.abs(identity)) < 1e-14
    excess = {k: channels['momentum'][k] - channels['energy'][k] for k in channels['energy']}
    denominator = float(excess['native_change_J'][-1])
    result = dict(new_pic_simulation=False, frame='boost',
        definition='Native staggered positive field-square norm, not a proof of physical energy.',
        denominator_J=denominator,
        final_excess={k: dict(J=float(v[-1]), percent=100 * float(v[-1]) / denominator)
                      for k, v in excess.items()})
    (root / 'CHANNEL_SUMMARY.json').write_text(json.dumps(result, indent=2) + '\n')
    with (root / 'CHANNEL_SUMMARY.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['channel', 'momentum_minus_energy_J', 'percent_of_final_native_difference'])
        w.writerows((k, v['J'], v['percent']) for k, v in result['final_excess'].items())

    os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-field-origin-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.3), layout='constrained')
    colors = {'energy': '#23678c', 'momentum': '#b17c25'}
    names = {'energy': 'Energy coupled', 'momentum': 'Momentum matched'}
    for case, d in data.items():
        step = np.arange(len(d['native_field_J']))
        axes[0, 0].plot(step, d['native_field_J'] * 1e6, color=colors[case], label=names[case])
        ratio = np.diff(d['time_electric_J']) / d['native_field_J'][1:]
        axes[1, 1].plot(step[1::10], ratio[::10], color=colors[case], label=names[case])
    axes[0, 0].set(xlabel='Maxwell step', ylabel='Native field norm (microjoule)',
                   title='A. The field difference already exists in boost frame')
    axes[0, 0].legend(fontsize=9)
    for k, label in [('native_change_J', 'Native field difference'),
                     ('numerical_time_cowan_J', 'Time + Cowan terms'),
                     ('current_input_J', 'Midpoint current work'),
                     ('curl_remap_J', 'Curl exchange + remap')]:
        axes[0, 1].plot(np.arange(len(excess[k]))[::10],
                        (100 * excess[k] / denominator)[::10], label=label)
    axes[0, 1].set(xlabel='Maxwell step', ylabel='Percent of final native difference',
                   title='B. Independent discrete-channel decomposition')
    axes[0, 1].legend(fontsize=8)
    with (root / 'independent_source/ADMISSION_CONTINUITY.csv').open() as f:
        r = list(csv.DictReader(f))
    node = np.array([int(v['node']) for v in r])
    full = np.array([float(v['full_continuity_residual_over_q']) for v in r])
    clipped = np.array([float(v['clipped_continuity_residual_over_q']) for v in r])
    keep = (node >= 0) & (node <= 6)
    axes[1, 0].bar(node[keep], clipped[keep], width=.55, color='#b17c25', label='Current center clipping')
    axes[1, 0].plot(node[keep], full[keep], 'o', color='#23678c', label='Full source trajectory')
    axes[1, 0].axhline(0, color='#777777', lw=.7)
    axes[1, 0].set(xlabel='Old-window charge node k (retained k >= 1)',
                   ylabel='Continuity defect / particle charge',
                   title='C. Admission unit probe: -0.2 to +0.3148 cells')
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].axhline(2, color='#777777', ls=':', label='Flip-mode indicator under E/B equipartition')
    axes[1, 1].set(xlabel='Maxwell step', ylabel='Per-step electric-square term / native norm',
                   title='D. Strong time-grid-scale content in both versions')
    axes[1, 1].legend(fontsize=8, loc='lower right')
    for ax in axes.flat:
        ax.grid(axis='y', alpha=.18)
    fig.suptitle('Field-origin diagnosis: recorded channels and a source-admission defect\n'
                 'Channel closure identifies bookkeeping; the production excitation chain is not yet replay-verified',
                 fontsize=12)
    fig.savefig(root / 'field_origin_diagnosis.png', dpi=200)
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root')
    main(p.parse_args().root)
