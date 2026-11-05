"""Closed-result BF comparison for legacy and actual-substep plane events."""
import argparse
import csv
import json
import os
from pathlib import Path
import re

import h5py
import numpy as np
from scipy.constants import c

from prepare import require, save
from analyze_closed import verify
import report_discrete


def bunching(records, wavelength, first_slice=None, end_slice=None):
    if not len(records):
        return None
    ids = records['source_id']
    mask = np.ones(len(records), dtype=bool)
    slices = (ids-1)//1024
    if first_slice is not None:
        mask &= slices >= first_slice
        mask &= slices < end_slice
    a = records[mask]
    if not len(a):
        return None
    w = a['weight'].astype(np.longdouble)
    t = a['time_s'].astype(np.longdouble)
    anchor = np.median(t)
    angle = np.asarray(np.remainder((t-anchor)*np.longdouble(c/wavelength), 1)*np.longdouble(2*np.pi), float)
    phase = np.exp(1j*angle)
    b = np.sum(w*phase)/np.sum(w)
    groups = np.asarray((a['source_id']-1)//(1024*10), dtype=int)
    groups -= groups.min()
    q = np.bincount(groups, weights=np.asarray(w, float))
    reals = np.bincount(groups, weights=np.asarray(w, float)*phase.real)
    imags = np.bincount(groups, weights=np.asarray(w, float)*phase.imag)
    good = q > 0
    local = (reals[good]+1j*imags[good])/q[good]
    return dict(records=len(a), represented_electrons=float(w.sum()),
        global_BF=float(abs(b)), local_RMS_BF=float(np.sqrt(np.sum(q[good]*abs(local)**2)/q.sum())),
        phase_anchor_s=float(anchor), local_bin='10 original optical slices; charge-weighted RMS')


def analyze_case(directory, plan):
    directory = Path(directory)
    files = verify(directory)
    rows, sums, ranks, steps = report_discrete.load(directory/'ledger', directory/'ledger/energy-ledger.h5')
    require(abs(steps-plan['expected_steps']) <= 2, 'Configured stop differs unexpectedly')
    require(ranks == plan['ranks'], 'MPI layout differs')
    _, _, ledger = report_discrete.calculate(rows, sums, steps)
    events = {}
    source_input = Path(plan['input'])
    with h5py.File(source_input, 'r') as f:
        input_ids = np.sort(f['particles/records']['source_id'])
    for name in ['entry', 'und-entry', 'mid', 'p1']:
        events[name] = {}
        for method, suffix, groupname in [
                ('legacy', '-particles.h5', 'particle_plane'),
                ('substeps', '-particles-substeps.h5', 'particle_substep_plane')]:
            path = directory/'detectors'/(name+suffix)
            with h5py.File(path, 'r') as f:
                g = f[groupname]
                require(int(g['complete'][()]) == 1, 'Unclosed particle plane: '+str(path))
                a = g['records'][:]
                require(len(a) == int(g['committed_records'][()]), 'Particle commit count differs')
            require(len(a) == len(np.unique(a['source_id'])), 'Duplicate source ID: '+str(path))
            require(np.all(a['crossing_direction'] == 1), 'Reverse crossings in control case')
            if plan['require_all_particles_through_planes']:
                require(np.array_equal(np.sort(a['source_id']), input_ids), 'Incomplete cohort: '+str(path))
            else:
                require(np.all(np.isin(a['source_id'], input_ids)), 'Unknown source ID')
            a = np.sort(a, order='source_id')
            events[name][method] = dict(path=str(path), records=len(a),
                full=bunching(a, plan['wavelength_m']),
                core=bunching(a, plan['wavelength_m'], plan['core_first_slice_inclusive'], plan['core_end_slice_exclusive']))
            if method == 'legacy':
                legacy = a
            else:
                common, ia, ib = np.intersect1d(legacy['source_id'], a['source_id'], return_indices=True)
                if len(common):
                    difference = a['time_s'][ib].astype(np.longdouble)-legacy['time_s'][ia].astype(np.longdouble)
                    phase = difference*np.longdouble(2*np.pi*c/plan['wavelength_m'])
                    events[name]['time_difference_substeps_minus_legacy'] = dict(common_records=len(common),
                        mean_s=float(np.mean(difference)), rms_s=float(np.sqrt(np.mean(difference*difference))),
                        max_abs_s=float(np.max(abs(difference))), rms_optical_phase_rad=float(np.sqrt(np.mean(phase*phase))))
    log = (directory/'aprl.log').read_text()
    samples = [(int(n), float(t)) for n,t in re.findall(r'phase=progress step=(\d+) elapsed_s=([\d.eE+-]+)', log)]
    retained = [(n,t) for n,t in samples if n >= 200]
    speed = None
    if len(retained) >= 2:
        speed = (retained[-1][1]-retained[0][1])/(retained[-1][0]-retained[0][0])
    source_logs = sorted((directory/'reference').glob('source-coupling-rank-*.csv'))
    require(len(source_logs) == ranks, 'Missing source-coupling logs')
    for path in source_logs:
        require(path.with_suffix('.complete').read_text().split() == ['1', str(steps)], 'Unclosed source-coupling log')
        a = np.genfromtxt(path, delimiter=',', names=True)
        require(np.all(a['self_force_defect'] < 1e-13), 'Coupled force did not converge')
        require(np.all(a['relative_fixed_point_work_defect'] < 1e-8), 'Coupled work did not converge')
    result = dict(program_integrity_pass=True, physical_acceptance=False,
        stage=plan['stage'], steps=steps, ranks=ranks, closed_scientific_files=files,
        strict_ledger=ledger, planes=events, progress_seconds_per_step_200_to_end=speed,
        native_field_final_J=float(sums['native_end_J'][-1]),
        initial_field_J=float(sums['native_begin_J'][0]))
    save(directory/'RESULT.json', result)
    return result


def compare(root):
    root = Path(root)
    result = {}
    for mode in ['off', 'on']:
        directory = root/mode
        data = json.loads((directory/'RESULT.json').read_text())
        require(data['program_integrity_pass'] and data['stage'] == 'full', 'Control case not complete')
        result[mode] = data
    table = []
    for plane in ['entry', 'und-entry', 'mid', 'p1']:
        for scope in ['full', 'core']:
            for method in ['legacy', 'substeps']:
                row = dict(plane=plane, scope=scope, time_reconstruction=method)
                for mode in ['off', 'on']:
                    row[mode] = result[mode]['planes'][plane][method][scope]
                table.append(row)
    save(root/'COMPARISON.json', dict(complete=True, physical_acceptance=False, rows=table,
        counterfactual='Same selected actual phase-space events and same grid; conservative-source linear versus linear-notch.',
        caution='Central hard crop and finite record are not an identical replica of the original full bunch. Baseline must reproduce the startup feature before interpreting its suppression.'))
    flat = []
    for row in table:
        for mode in ['off', 'on']:
            values = row[mode]
            flat.append(dict(plane=row['plane'], scope=row['scope'],
                time_reconstruction=row['time_reconstruction'], mode=mode,
                global_BF=values['global_BF'], local_RMS_BF=values['local_RMS_BF'],
                represented_electrons=values['represented_electrons']))
    with (root/'bunching_comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-critical-control-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.6), sharex=True, layout='constrained')
    z = [0., .01, 1.088, 2.166]
    for i, scope in enumerate(['full', 'core']):
        for j, (key, label) in enumerate([('global_BF', 'Global BF'), ('local_RMS_BF', 'Local RMS BF')]):
            for mode, color in [('off', '#24648c'), ('on', '#b18428')]:
                for method, style in [('legacy', '-'), ('substeps', '--')]:
                    values = [next(row[mode][key] for row in table if row['plane']==plane
                              and row['scope']==scope and row['time_reconstruction']==method)
                              for plane in ['entry', 'und-entry', 'mid', 'p1']]
                    axes[i,j].semilogy(z, values, style, marker='o', ms=4, color=color,
                                      label=mode+' / '+method)
            axes[i,j].set_title(('Full 231-slice crop' if scope=='full' else 'Central 77 slices')+': '+label)
            axes[i,j].set_ylabel(label)
            axes[i,j].grid(alpha=.2)
            axes[i,j].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel('Fixed laboratory plane z (m)')
    fig.suptitle('First-segment conservative-source control\nSame selected input; linear versus linear-notch; no near-field subtraction', fontsize=12)
    fig.savefig(root/'bunching_comparison.png', dpi=200)
    fig.savefig(root/'bunching_comparison.pdf')
    plt.close(fig)
    return table


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    a = p.parse_args()
    print(json.dumps(compare(a.root), indent=2))
