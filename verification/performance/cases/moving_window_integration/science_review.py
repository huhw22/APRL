#!/usr/bin/env python3
"""Review completed sparse/dense same-input pairs; never starts a simulation."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('moving_science_helpers', HERE/'analyze.py')
review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(review)
LD = np.longdouble
PHOTON_EV = LD('91.5')
# Exact SI definitions h and e; derive hbar in extended precision.
HBAR_EV_S = LD('6.62607015e-34')/(2*np.arccos(LD(-1))*LD('1.602176634e-19'))


def bunching(reference, candidate):
    weights = np.abs(np.asarray(reference['charge_C'], dtype=LD))
    times = np.asarray(reference['time_s'], dtype=LD)
    origin = np.sum(weights*times)/weights.sum()
    output = {'photon_energy_eV': float(PHOTON_EV), 'shared_time_origin_s': float(origin),
        'definition': '|sum(abs(charge)*exp(i*omega*(t-tref)))|/sum(abs(charge)); shared long-double-centred time within each pair'}
    complex_values = []
    for name, records in [('fixed', reference), ('moving', candidate)]:
        review.require(np.array_equal(reference['source_id'], records['source_id']) and
            np.array_equal(reference['charge_C'], records['charge_C']), 'Bunching cohorts/weights differ')
        phase = PHOTON_EV/HBAR_EV_S*(np.asarray(records['time_s'], dtype=LD)-origin)
        real = np.sum(weights*np.cos(phase))/weights.sum()
        imaginary = np.sum(weights*np.sin(phase))/weights.sum()
        amplitude = np.hypot(real, imaginary)
        review.require(amplitude <= 1+LD('1e-15'), 'Unphysical bunching amplitude')
        complex_values.append((real, imaginary))
        output[name] = {'bf': float(amplitude), 'real': float(real), 'imaginary': float(imaginary)}
    output['moving_minus_fixed_bf'] = output['moving']['bf']-output['fixed']['bf']
    output['complex_difference_magnitude'] = float(np.hypot(
        complex_values[1][0]-complex_values[0][0], complex_values[1][1]-complex_values[0][1]))
    return output


def pair(root):
    plan = json.loads((root/'plan.json').read_text())
    with h5py.File(root/'inputs/held.h5') as stream:
        source = stream['particles/records'][:]
    fixed = review.inspect_case(root/'fixed', plan, source)
    moving = review.inspect_case(root/'moving', plan, source)
    compared = review.compare(fixed, moving, 'same-input fixed-to-moving; no cross-density timing ratio')
    local_traces = {}
    first, attrs = fixed['files']['detectors/exit-field.h5']
    second = moving['files']['detectors/exit-field.h5'][0]
    x = attrs['x_first_m']+np.arange(first['fields'].shape[2])*attrs['dx_m']
    y = attrs['y_first_m']+np.arange(first['fields'].shape[1])*attrs['dy_m']
    i, j = int(np.argmin(abs(x+37.5e-6))), int(np.argmin(abs(y+22.5e-6)))
    for field in ('electric_V_per_m', 'magnetic_T'):
        a, b = first['fields'][field][:, j, i], second['fields'][field][:, j, i]
        norm = np.linalg.norm(a, axis=1)
        error = np.linalg.norm(b-a, axis=1)
        # Preserve weak-signal locations, but explicitly mark the norm floor.
        active = norm>norm.max()*1e-6
        local_traces[field] = {'fixed_vector_norm': norm.tolist(), 'difference_vector_norm': error.tolist(),
            'relative_error': [float(error[k]/norm[k]) if active[k] else None for k in range(len(norm))],
            'norm_floor_relative_to_peak': 1e-6}
    return {'root': str(root), 'ranks': plan['ranks'], 'threads': plan['threads'], 'records': plan['source']['records'],
        'macro_electrons': plan['source']['macro_electrons'], 'geometry': plan['geometry'],
        'fixed': fixed['summary'], 'moving': moving['summary'], 'comparison': compared,
        'bf_91p5eV': {name: bunching(fixed['planes'][name], moving['planes'][name]) for name in ('release', 'aux', 'exit')},
        'historical_sparse_pixel': {'x_m': float(x[i]), 'y_m': float(y[j]), 'traces': local_traces}}


def plot(report, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 8, 'pdf.fonttype': 42})
    figure, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    colours = {'sparse': '#276ba4', 'dense': '#c15a19'}
    energy_error_values = []
    for name in ('sparse', 'dense'):
        item = report[name]
        fields = item['comparison']['field_regions']['fields']
        electric = fields['electric_V_per_m']['physical_xy']
        trace = electric['time_traces']
        times = np.asarray(trace['relative_lab_time_fs'])
        for mode, style in [('fixed', '-'), ('moving', '--')]:
            value = np.asarray(trace[mode+'_vector_rms'])**2*electric['area_m2']
            axes[0, 0].plot(times, value, style, color=colours[name], label=name+' '+mode)
        for field, label, style in [('electric_V_per_m', 'E', '-'), ('magnetic_T', 'B', '--')]:
            region = fields[field]['physical_xy']['time_traces']
            norm = np.asarray(region['fixed_vector_rms'])
            relative = np.asarray(region['difference_vector_rms'])/norm*100
            axes[0, 1].plot(times, np.where(relative>0, relative, np.nan), style,
                color=colours[name], label=name+' '+label)
            local = np.array([np.nan if value is None else value*100
                for value in item['historical_sparse_pixel']['traces'][field]['relative_error']])
            axes[1, 0].plot(times, np.where(local>0, local, np.nan), style,
                color=colours[name], label=name+' '+label)
        for k, plane in enumerate(('aux', 'exit')):
            value = item['comparison']['planes'][plane]['rms_candidate_minus_reference_K_eV']
            energy_error_values.append(value)
            axes[1, 1].bar(k+(-.18 if name=='sparse' else .18), value, width=.32,
                color=colours[name], label=name if k==0 else None)
            axes[1, 1].text(k+(-.18 if name=='sparse' else .18), value*1.12, f'{value:.3g}', ha='center', va='bottom', fontsize=8)
        for axis, index in [('x', 0), ('y', 1)]:
            location=index+(0 if name=='sparse' else 3)
            for mode, marker, offset in [('fixed', 'o', -.10), ('moving', 'x', .10)]:
                value=item[mode]['particle_planes']['aux']['transverse'][axis]['sigma_minus_ballistic_m']*1e12
                axes[2, 0].scatter(location+offset, value, color=colours[name], marker=marker,
                    label=name+' '+mode if axis=='x' else None, zorder=3)
        for mode, marker, style in [('fixed', 'o', '-'), ('moving', 'x', '--')]:
            values=[item['bf_91p5eV'][plane][mode]['bf'] for plane in ('release', 'aux', 'exit')]
            axes[2, 1].plot(range(3), values, style, marker=marker, color=colours[name], label=name+' '+mode)
    axes[0, 0].set(title='(a) Physical aperture: total electric field', ylabel=r'$\int |E|^2\,dA$ (V$^2$)', yscale='log')
    axes[0, 1].set(title='(b) Physical aperture: field difference', ylabel=r'$\|F_{moving}-F_{fixed}\|_2/\|F_{fixed}\|_2$ (%)', yscale='log')
    axes[1, 0].set(title='(c) Same sparse-error pixel: x=-37.5, y=-22.5 um', ylabel='Local vector relative difference (%)', yscale='log')
    axes[1, 1].set(title='(d) Per-ID kinetic-energy difference', ylabel='RMS(moving - fixed) (eV/electron)', yscale='log', xticks=[0,1], xticklabels=['Auxiliary plane', 'Exit plane'])
    axes[2, 0].set(title='(e) Auxiliary-plane expansion after ballistic drift', ylabel=r'$\sigma-\sigma_{ballistic}$ (pm)', xticks=[0,1,3,4], xticklabels=['Sparse x','Sparse y','Dense x','Dense y'])
    axes[2, 1].set(title='(f) 91.5 eV bunching amplitude', ylabel='|weighted phase sum| / total weight', yscale='log', xticks=[0,1,2], xticklabels=['Release','Auxiliary','Exit'])
    axes[1, 1].set_ylim(min(energy_error_values)*.6, max(energy_error_values)*1.8)
    dense_trace = report['dense']['comparison']['field_regions']['fields']['electric_V_per_m']['physical_xy']['time_traces']
    cohort_last = report['dense']['comparison']['fields']['through_cohort_exit_samples']-1
    cohort_time = dense_trace['relative_lab_time_fs'][cohort_last]
    for axis in (axes[0,0],axes[0,1],axes[1,0]):
        axis.set_xlabel('Lab time from first field sample (fs)')
        axis.axvline(cohort_time, color='0.5', linewidth=.8, linestyle=':', label='Dense cohort passed plane')
    for axis in axes.flat:
        axis.grid(alpha=.22, which='major')
        axis.legend(loc='best')
    layout = 'Single-rank' if report['sparse']['ranks']==report['dense']['ranks']==1 else 'Same-layout'
    figure.suptitle(layout+' comparison: sparse N='+str(report['sparse']['records'])+'; dense N='+format(report['dense']['records'], ',')+
        '\nMatched fixed/moving input; native E and staggered B, not radiation-only power', fontsize=12)
    figure.savefig(output/'SPARSE_DENSE_SCIENCE.png', dpi=180)
    figure.savefig(output/'SPARSE_DENSE_SCIENCE.pdf')
    plt.close(figure)


def main(args):
    output=args.output.resolve()
    review.require(not output.exists(), 'Use a fresh analysis directory')
    report={'format':'aprl-moving-sparse-dense-science-v1',
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'helper_sha256':hashlib.sha256((HERE/'analyze.py').read_bytes()).hexdigest(),
        'sparse':pair(args.sparse_root.resolve()), 'dense':pair(args.dense_root.resolve()),
        'limits':['Each fixed/moving pair has the same input. Sparse and dense sampling differs; do not turn cross-N runtimes into acceleration.',
            'Source ID2 remains an explicit tail witness, not a random-tail statistic.',
            'Short one-period diagnostics do not establish radiation-only fields or far-field spectral convergence.',
            'These single-rank comparisons do not establish multi-node scientific equivalence.',
            'All E and B norms are separate in their own units; amplitude-squared area sums are not Poynting power.',
            'The previous sparse maximum is inspected alongside each new global maximum. No automatic physical PASS threshold is invented.']}
    output.mkdir(parents=True)
    review.save(output/'SPARSE_DENSE_SCIENCE.json', report)
    if not args.no_plot:
        plot(report,output)
    print(output/'SPARSE_DENSE_SCIENCE.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sparse-root',type=Path,required=True)
    parser.add_argument('--dense-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--no-plot',action='store_true')
    main(parser.parse_args())
