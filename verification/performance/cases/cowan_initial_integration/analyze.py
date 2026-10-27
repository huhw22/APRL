#!/usr/bin/env python3
"""Plot completed small Cowan integration runs from their original HDF5 data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

import h5py
_MPL_CACHE = tempfile.TemporaryDirectory(prefix='cowan-plot-cache-')
os.environ.setdefault('MPLCONFIGDIR', _MPL_CACHE.name)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml

C = 299792458.0
E = 1.602176634e-19
REST_EV = 510998.9499961642
MODELS = ('legacy', 'cowan')
COLORS = {'legacy': '#325777', 'cowan': '#c47424'}
LABELS = {'legacy': 'Legacy initializer', 'cowan': 'Cowan initializer'}
CASES = ('gaussian-static-pec', 'gaussian-moving-cpml')
TITLES = ('Static in boost frame / PEC', 'Moving in boost frame / CPML')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def records(path):
    with h5py.File(path) as file:
        group = next(iter(file.values()))
        require(int(group['complete'][()]) == 1, 'Incomplete output: ' + str(path))
        array = group['records'][:]
        require(int(group['committed_records'][()]) == len(array), 'Uncommitted records')
    return array


def particles(path):
    with h5py.File(path) as file:
        require(int(file['header'][2]) == 1, 'Incomplete checkpoint')
        result = file['blocks/particles/data'][:]
    return np.sort(result, order='id')


def radius(xy, weights):
    mean = np.average(xy, weights=weights, axis=0)
    return float(np.sqrt(np.sum(np.average((xy-mean)**2, weights=weights, axis=0))))


def pair_summary(old, new):
    old, new = np.asarray(old, dtype=float), np.asarray(new, dtype=float)
    require(old.shape == new.shape, 'Comparison shapes differ')
    relative = np.full_like(old, np.nan)
    np.divide(new-old, old, out=relative, where=old != 0)
    return {'legacy_final': float(old[-1]), 'cowan_final': float(new[-1]),
            'new_minus_old_final': float(new[-1]-old[-1]),
            'new_minus_old_max_abs': float(np.max(np.abs(new-old))),
            'relative_difference_final': None if old[-1] == 0 else float(relative[-1]),
            'relative_difference_max_abs_nonzero_reference':
                float(np.nanmax(np.abs(relative))) if np.any(np.isfinite(relative)) else None,
            'undefined_relative_samples': int(np.count_nonzero(old == 0))}


def load_model(directory, expected_particles):
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    manifest = yaml.safe_load((directory / 'run-manifest.yaml').read_text())
    require(yaml.safe_load(manifest['configuration_yaml']) == card, 'Executed configuration mismatch')
    ledger_path = directory / 'ledger/energy-ledger.h5'
    ledger = records(ledger_path)
    require(np.all(ledger['active_macroparticles'] == expected_particles) and
            np.all(ledger['removed_macroparticles'] == 0), 'Particle cohort changed')
    require(np.all(ledger['prescribed_work_J'] == 0), 'Unexpected prescribed work')
    require(np.all(ledger['active_represented_electrons'] == ledger['active_represented_electrons'][0]),
            'Represented charge changed')
    gamma = card['mesh']['boost_gamma']
    beta = np.sqrt(1-1/gamma**2)
    dt = card['mesh']['cell_size'][2] * 1e-6 / C
    origin = re.search(r'Lorentz synchronization: anchor_lab_z=([^ ]+) m, origin_lab_t=([^ ]+) s',
                       (directory / 'aprl.log').read_text())
    require(origin is not None, 'Missing Lorentz origin')
    z0, t0 = map(float, origin.groups())
    trajectory_paths = sorted((directory / 'trajectory').glob('*.h5'))
    trajectory = np.concatenate([records(path) for path in trajectory_paths])
    trajectory = trajectory[trajectory['event_type'] == 0]
    box_time = gamma*(trajectory['time_s']-t0) - gamma*beta*(trajectory['position_m'][:, 2]-z0)/C
    steps = np.rint(box_time/dt).astype(np.int64)
    require(np.max(np.abs(box_time/dt-steps)) < 1e-6, 'Cannot recover common boost sampling time')
    initial_path = directory / 'checkpoints/step-000000000000.h5'
    initial = particles(initial_path)
    weights0 = np.abs(initial['charge'])
    initial_radius = radius(np.column_stack([initial['x'], initial['y']]), weights0)
    rms_steps, rms, excluded = [], [], []
    for step in np.unique(steps):
        sample = trajectory[steps == step]
        sample = np.sort(sample, order='particle_id')
        if len(sample) != expected_particles:
            excluded.append({'step': int(step), 'particles': len(sample)})
            continue
        require(np.array_equal(sample['particle_id'], initial['id']), 'Trajectory particle IDs changed')
        rms_steps.append(int(step))
        rms.append(radius(sample['position_m'][:, :2], np.abs(sample['charge_C'])))
    require(len(rms) > 2, 'Insufficient complete-cohort trajectory samples')
    snapshots = {}
    initial_u = np.column_stack([initial[key] for key in ('ux', 'uy', 'uz')]).astype(np.longdouble)
    initial_gamma = np.sqrt(1 + np.sum(initial_u**2, axis=1))
    for path in sorted((directory / 'checkpoints').glob('step-*.h5')):
        p = particles(path)
        require(np.array_equal(initial['id'], p['id']), 'Snapshot particle IDs changed')
        u = np.column_stack([p[key] for key in ('ux', 'uy', 'uz')]).astype(np.longdouble)
        g = np.sqrt(1+np.sum(u**2, axis=1))
        delta = np.sum((u-initial_u)*(u+initial_u), axis=1)/(g+initial_gamma)
        delta_lab = gamma*(delta+beta*(u[:, 2]-initial_u[:, 2]))
        step = int(path.stem.split('-')[-1])
        snapshots[str(step)] = {'held': int(np.count_nonzero(p['frozen'])), 'particles': len(p),
            'rms_radius_m': radius(np.column_stack([p['x'], p['y']]), np.abs(p['charge'])),
            'stable_delta_mean_box_energy_eV': float(np.average(delta, weights=weights0)*REST_EV),
            'stable_delta_mean_lab_energy_eV': float(np.average(delta_lab, weights=weights0)*REST_EV)}
    electrons = float(ledger['active_represented_electrons'][0])
    mean_delta_box = ((ledger['particle_kinetic_active_J'].astype(np.longdouble)-
                       np.longdouble(ledger['particle_kinetic_active_J'][0]))/(electrons*E)).astype(float)
    mean_delta_lab = ((ledger['mean_gamma_lab']-ledger['mean_gamma_lab'][0])*REST_EV).astype(float)
    for step, snapshot in snapshots.items():
        index = np.flatnonzero(ledger['step'] == int(step))
        require(len(index) == 1, 'Snapshot missing ledger sample')
        require(abs(mean_delta_box[index[0]]-snapshot['stable_delta_mean_box_energy_eV']) < 1e-5 and
                abs(mean_delta_lab[index[0]]-snapshot['stable_delta_mean_lab_energy_eV']) < 1e-5,
                'Ledger disagrees with stable fixed-cohort energy change')
    result = {
        'steps': ledger['step'].astype(int), 'time_box_fs': ledger['time_box_s']*1e15,
        'delta_mean_box_energy_eV': mean_delta_box, 'delta_mean_lab_energy_eV': mean_delta_lab,
        'delta_boost_beam_kinetic_J': ledger['particle_kinetic_active_J']-ledger['particle_kinetic_active_J'][0],
        'sigma_lab_eV': ledger['sigma_gamma_lab']*REST_EV,
        'sigma_box_eV': ledger['sigma_gamma_box']*REST_EV,
        'relative_lab_spread': ledger['relative_energy_spread_lab'],
        'field_energy_J': ledger['field_energy_interior_J'], 'residual_J': ledger['balance_residual_J'],
        'closure_valid': ledger['closure_valid'].astype(int),
        'rms_steps': np.array(rms_steps), 'rms_time_fs': np.array(rms_steps)*dt*1e15,
        'rms_radius_m': np.array(rms), 'rms_growth_percent': (np.array(rms)/initial_radius-1)*100,
        'initial_radius_m': initial_radius, 'excluded_partial_trajectory_samples': excluded,
        'snapshots': snapshots, 'particles': expected_particles, 'represented_electrons': electrons,
        'source_revision': manifest['source_revision'],
        'sources': {str(path): digest(path) for path in
                    [ledger_path, initial_path, directory/'run.yaml', directory/'run-manifest.yaml']+trajectory_paths}}
    return result


def summarize(old, new):
    for key in ('steps', 'rms_steps'):
        require(np.array_equal(old[key], new[key]), 'Different physical sampling steps')
    names = ('rms_radius_m', 'rms_growth_percent', 'delta_mean_box_energy_eV',
             'delta_mean_lab_energy_eV', 'delta_boost_beam_kinetic_J', 'sigma_lab_eV',
             'sigma_box_eV', 'relative_lab_spread', 'field_energy_J', 'residual_J')
    return {'comparisons': {name: pair_summary(old[name], new[name]) for name in names},
            'ledger_final_time_box_fs': float(old['time_box_fs'][-1]),
            'trajectory_final_time_box_fs': float(old['rms_time_fs'][-1]),
            'legacy': old, 'cowan': new}


def draw_main(cases, output):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.alpha': .18, 'grid.linewidth': .6})
    fig, axes = plt.subplots(4, 4, figsize=(17, 11), gridspec_kw={'height_ratios': [1.6, 1, 1.6, 1]})
    columns = [('rms_growth_percent', 'Transverse RMS growth', 'Growth (%)', 'rms_time_fs'),
               ('delta_mean_box_energy_eV', 'Boost mean energy change', r'$\Delta\langle K\rangle_{box}$ (eV/electron)', 'time_box_fs'),
               ('delta_mean_lab_energy_eV', 'Lab mean energy change', r'$\Delta\langle E\rangle_{lab}$ (eV/electron)', 'time_box_fs'),
               ('sigma_lab_eV', 'Lab energy spread', r'$\sigma(E)_{lab}$ (eV)', 'time_box_fs')]
    for case_index, (name, title) in enumerate(zip(CASES, TITLES)):
        row = case_index*2
        data = cases[name]
        for col, (key, label, unit, tkey) in enumerate(columns):
            ax, diff = axes[row, col], axes[row+1, col]
            for model in MODELS:
                m = data[model]
                ax.plot(m[tkey], m[key], color=COLORS[model],
                        ls='--' if model == 'legacy' else '-', lw=1.7, label=LABELS[model])
            ax.set_title(label if row == 0 else '', loc='left', fontsize=11)
            ax.set_ylabel(unit)
            ax.set_xlim(0, 216)
            ax.set_ylim(bottom=0)
            old, new = data['legacy'], data['cowan']
            diff.plot(new[tkey], new[key]-old[key], color='#58505f', lw=1.5)
            diff.axhline(0, color='#888888', lw=.7)
            diff.set_ylabel('New − old\n' + ('percentage points' if col == 0 else 'eV'))
            diff.set_xlim(0, 216)
            diff.set_xlabel('Boost time (fs)')
            diff.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useMathText=True)
        axes[row, 0].text(0, 1.12, title, transform=axes[row, 0].transAxes,
                          fontsize=12, fontweight='bold')
        axes[row, 3].legend(loc='upper left', frameon=False, fontsize=9)
    for col in range(4):
        maximum = max(axes[0, col].get_ylim()[1], axes[2, col].get_ylim()[1])
        axes[0, col].set_ylim(0, maximum)
        axes[2, col].set_ylim(0, maximum)
    fig.suptitle('Full-program initializer comparison · 1 MPI rank · fixed 4096-particle cohort',
                 x=.06, ha='left', fontsize=15, fontweight='bold')
    fig.text(.06, .942, '1 pC Gaussian, fixed seed · 32 × 32 × 192 cells · no prescribed fields · same physical settings', fontsize=10)
    fig.text(.06, .014, 'Absolute panels share scales between cases; difference panels use their labeled scales. RMS ends at 210.15 fs; ledger at 215.15 fs.\n'
             'This compares two finite-domain initializers; it does not isolate radiation or measure a production speedup.', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.07, right=.985, bottom=.085, top=.88, hspace=.7, wspace=.42)
    for suffix in ('png', 'pdf'):
        fig.savefig(output / ('comparison.'+suffix), dpi=180)
    plt.close(fig)


def draw_held(data, output):
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), gridspec_kw={'height_ratios': [1.7, 1]})
    columns = [('rms_growth_percent', 'Complete-cohort RMS growth', 'Growth from fixed step-0 cohort (%)', 'rms_time_fs'),
               ('delta_mean_lab_energy_eV', 'Fixed-cohort lab energy change', 'Mean energy change (eV/electron)', 'time_box_fs'),
               ('sigma_lab_eV', 'Fixed-cohort lab energy spread', 'Energy spread (eV)', 'time_box_fs')]
    for col, (key, title, label, tkey) in enumerate(columns):
        for model in MODELS:
            m = data[model]
            mask = m['rms_steps' if tkey == 'rms_time_fs' else 'steps'] >= 7
            axes[0, col].plot(m[tkey][mask], m[key][mask], color=COLORS[model],
                              ls='--' if model == 'legacy' else '-', lw=1.7, label=LABELS[model])
        m, old = data['cowan'], data['legacy']
        mask = m['rms_steps' if tkey == 'rms_time_fs' else 'steps'] >= 7
        axes[1, col].plot(m[tkey][mask], (m[key]-old[key])[mask], color='#58505f', lw=1.5)
        axes[0, col].set_title(title, loc='left', fontsize=11)
        axes[0, col].set_ylabel(label)
        axes[1, col].set_ylabel('New − old\n' + ('percentage points' if col == 0 else 'eV'))
        axes[1, col].set_xlabel('Boost time (fs)')
        axes[1, col].ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useMathText=True)
    axes[0, 2].legend(frameon=False, fontsize=9)
    fig.suptitle('Held release supplement · fixed 64-particle cohort · full release by step 7', x=.07, ha='left', fontsize=14)
    fig.text(.07, .014, 'Step-5 trajectory sample contains only 49 particles and is excluded. Full-cohort trajectory samples start at step 9.\n'
             'Ledger counts all 64 particles even while held; closure_valid=0. This underresolved, low-charge fixture checks integration only.', fontsize=9)
    fig.subplots_adjust(left=.08, right=.98, bottom=.14, top=.86, hspace=.5, wspace=.4)
    for suffix in ('png', 'pdf'):
        fig.savefig(output / ('held_comparison.'+suffix), dpi=180)
    plt.close(fig)


def serializable(value):
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    result = json.loads((root / 'RESULT.json').read_text())
    require(result['integration_status'] == 'PASS', 'Integration campaign did not pass')
    require(result['ranks'] == 1, 'This figure is explicitly the one-MPI comparison')
    output.mkdir(parents=True, exist_ok=True)
    cases = {}
    for case in CASES + ('held-supplement',):
        count = 64 if case == 'held-supplement' else 4096
        old, new = [load_model(root / case / model, count) for model in MODELS]
        cases[case] = summarize(old, new)
    draw_main(cases, output)
    draw_held(cases['held-supplement'], output)
    review = {'scope': 'Single-MPI full-program initial-field comparison',
              'source_root': str(root), 'integration_result_sha256': digest(root/'RESULT.json'),
              'analysis_script': str(Path(__file__).resolve()),
              'analysis_script_sha256': digest(Path(__file__).resolve()),
              'cases': cases,
              'definitions': {
                  'transverse_rms': 'sqrt(weighted var(x)+weighted var(y)) at common boost time',
                  'relative_lab_spread': 'sigma(gamma_lab)/mean(gamma_lab); total-energy denominator',
                  'energy_difference': 'New minus legacy; absolute changes referenced to the same fixed step-0 particle cohort',
                  'zero_reference': 'Relative differences are undefined at zero reference; use absolute eV difference',
                  'residual': 'Existing ordinary field-energy ledger diagnostic, not the exact Cowan discrete energy invariant'},
              'limitations': ['Finite-domain short Gaussian tests, not full-scale production or pure radiation.',
                  'No zero-field or independent physical reference is run here: new-old agreement alone does not isolate startup excess energy.',
                  'Held fixture is underresolved transversely and has only 64 represented electrons; use as a functional release/restart check.',
                  'Held ledger closure_valid=0; do not interpret its balance residual as a conservation test.',
                  'RMS uses complete cohorts only; no comparison mixes early partially released trajectory samples.']}
    (output / 'PHYSICS_REVIEW.json').write_text(json.dumps(serializable(review), indent=2, allow_nan=False)+'\n')
    for case, values in cases.items():
        print(case, json.dumps(values['comparisons'], sort_keys=True))
    print(output / 'comparison.png')


if __name__ == '__main__':
    main()
