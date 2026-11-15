"""Read-only comparison of closed, independent w355 six-segment runs."""
import csv
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if not (HERE / 'growth_compare.py').exists():
    sys.path.insert(0, str(HERE.parent / 'puffin_six_same_ledger'))
    sys.path.insert(0, str(HERE.parent / 'puffin_eleven_four_slip'))

import h5py
import numpy as np
import yaml
from common import C, E, read, save, need, sha, frozen, dump_at
from growth_compare import particle_stats, field_power, table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD = np.longdouble
REST_EV = LD('9.1093837139e-31') * LD(C)**2 / LD(E)


def identity(path):
    s = Path(path).stat()
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def close_check(root):
    puffin = root / 'six'
    pic = root / 'six-aprl/energy-binomial3'
    plan = read(root / 'six-aprl/PLAN.json')
    need(read(puffin / 'TRACKING_CLOSED.json')['complete'], 'Unclosed Puffin')
    need(read(puffin / 'COMPLETE.json')['complete'], 'Incomplete Puffin')
    frozen(puffin / 'FROZEN_PHYSICS.json')
    frozen(root / 'six-aprl/FROZEN.json')
    launch = read(pic / 'aprl.log.launch.json')
    need(launch['exit_code'] == 0, 'PIC physical process failed')
    log = (pic / 'aprl.log').read_text()
    for marker in ['Configured stop reached:', 'phase=solve-complete',
                   'phase=mpi-finalize-complete', 'phase=detector-close-complete']:
        need(marker in log, 'PIC missing final marker: ' + marker)
    need(not re.search(r'Simulation failure|Segmentation fault|OpenFabrics device', log),
         'PIC runtime error')
    steps = int(re.search(r'Moving window completed: steps=(\d+)', log)[1])
    need(steps == plan['steps'] + 1, 'Unexpected production step discrepancy')
    driver = (root / 'six-aprl/driver-energy-binomial3.log').read_text()
    need(driver.rstrip().endswith('AssertionError: Moving window completed: steps=' + str(plan['steps'])),
         'Unexpected driver failure')
    closed = {}
    for path in sorted((pic / 'detectors').glob('*.h5')):
        with h5py.File(path, 'r') as f:
            for name, g in f.items():
                need(int(g['complete'][()]) == 1, 'Unclosed detector: ' + str(path))
                committed = 'committed_samples' if 'fields' in g else 'committed_records'
                data = 'fields' if 'fields' in g else 'records'
                n = int(g[committed][()])
                need(n == len(g[data]), 'Uncommitted detector: ' + str(path))
                closed[str(path)] = dict(group=name, records=n)
    with h5py.File(pic / 'ledger/energy-ledger.h5', 'r') as f:
        g = f['energy_ledger']
        need(int(g['complete'][()]) == 1, 'Unclosed energy ledger')
        need(len(g['records']) == int(g['committed_records'][()]) == steps + 1,
             'Energy ledger record mismatch')
        need(int(g['records'][-1]['step']) == steps, 'Ledger final step mismatch')
    need((pic / 'worldtube/COMPLETE').read_text().split() == ['1', str(steps)],
         'Unclosed worldtube')
    for rank in range(16):
        for prefix in ['lab-field', 'paired-momentum', 'whole-cowan']:
            need((pic / 'lab-audit' / f'{prefix}-rank-{rank:05d}.complete').read_text().split()
                 == ['1', str(steps)], 'Unclosed audit rank')
    return dict(physical_tracking_complete=True, pic_process_exit_code=0,
                original_pic_batch_state='FAILED 1:0 in post-run expected-step assertion',
                expected_pic_steps=plan['steps'], actual_pic_steps=steps,
                all_detector_groups_closed=closed, all_48_lab_audit_markers_closed=True,
                worldtube_closed=True, ledger_records=steps + 1,
                frozen_inputs_cards_programs_verified=True, physical_acceptance=False)


def gamma_stats(path, kind, info):
    with h5py.File(path, 'r') as f:
        if kind == 'PIC':
            a = f['particle_substep_plane/records'].fields(
                ['gamma_lab', 'charge_C', 'source_id', 'crossing_direction'])[:]
            q = abs(a['charge_C']).astype(LD)
            gamma = a['gamma_lab'].astype(LD)
            need(np.all(a['crossing_direction'] == 1), 'Reverse particle crossings')
            need(len(np.unique(a['source_id'])) == len(a), 'Duplicate source crossings')
        else:
            a = f['electrons'][:]
            q = a[:, 6].astype(LD) * LD(info['scale']['npkbar']) * LD(info['scale']['qe_C'])
            gamma = a[:, 5].astype(LD) * LD(info['scale']['gamma_ref'])
    need(np.isfinite(gamma).all() and np.isfinite(q).all() and np.all(q > 0),
         'Nonfinite gamma/charge')
    mean = np.sum(q * gamma, dtype=LD) / q.sum(dtype=LD)
    return mean


def main(root, out):
    root, out = Path(root).resolve(), Path(out).resolve()
    need(root.name == 'fel_independent_w355_20261114_r1', 'Wrong source campaign')
    out.mkdir(parents=True, exist_ok=False)
    closure = close_check(root)
    save(out / 'CLOSURE_REVIEW.json', closure)
    puffin = root / 'six'
    pic = root / 'six-aprl/energy-binomial3'
    plan, info = read(puffin / 'PLAN.json'), read(puffin / 'INPUT_VALIDATED.json')
    ap = read(root / 'six-aprl/PLAN.json')
    need(info['input_sha256'] == ap['input_sha256'], 'Different input events')
    need(info['independent_fixed_N_source'] and info['zero_transverse_duplicates']
         and info['no_clones'], 'Unqualified source')
    card = yaml.safe_load((pic / 'run.yaml').read_text())
    need(plan['magnetic_lattice'] == card['sources']['magnetic_elements'], 'Different magnets')
    need(read(puffin / 'NATIVE_VALIDATED.json')['all_seven_columns_bitwise'], 'Import mismatch')
    half = LD(str(plan['duration_s'])) / 2
    low, high = LD('.1') * LD(str(math.sqrt(.2))) * half, LD('.9') * LD(str(math.sqrt(.2))) * half
    rows, initial, protected = [], [], {}
    initial_puffin = dump_at(puffin / 'main', 'electrons', 0)
    source_gamma = gamma_stats(initial_puffin, 'Puffin', info)
    for name, kind in [('Puffin', 'Puffin'), ('APRL EC+F3', 'PIC')]:
        pin = initial_puffin if kind == 'Puffin' else pic / 'detectors/und-entry-particles-substeps.h5'
        protected[str(pin)] = identity(pin)
        g0 = gamma_stats(pin, kind, info)
        initial.append(dict(method=name, segment=0, **particle_stats(pin, kind, info, plan, low, high),
                            mean_gamma_lab=float(g0), mean_loss_eV=float((source_gamma - g0) * REST_EV)))
        for item in plan['outputs']:
            k, step = item['segment'], item['step']
            p = dump_at(puffin / 'main', 'electrons', step) if kind == 'Puffin' else pic / 'detectors' / f'p{k}-particles-substeps.h5'
            f = dump_at(puffin / 'main', 'aperp', step) if kind == 'Puffin' else pic / 'detectors' / f'p{k}.h5'
            for path in [p, f]: protected[str(path)] = identity(path)
            ps = particle_stats(p, kind, info, plan, low, high)
            need(abs(ps['actual_z_m'] - item['actual_lab_z_m']) < 1e-11, 'Exit position mismatch')
            fp = field_power(f, kind, info, plan, LD(ps['median_retarded_time_s']), low, high)
            gamma = gamma_stats(p, kind, info)
            loss = (source_gamma - gamma) * REST_EV
            full = ps['records'] == plan['sample_records'] and abs(ps['charge_fraction_of_original'] - 1) < 2e-12
            row = dict(method=name, segment=k, **ps, **fp, mean_gamma_lab=float(gamma),
                       mean_loss_eV=float(loss), mean_loss_from_own_magnetic_entry_eV=float((g0-gamma)*REST_EV),
                       full_cohort_available=full,
                       full_cohort_particle_loss_J=float(loss * LD(plan['source_charge_C'])) if full else None)
            rows.append(row)
            print('SEGMENT', name, k, 'power_W', fp['Hann_pm5pct_power_Eproxy_W'],
                  'local_BF', ps['central_10lambda_rms_b1'], 'loss_eV', float(loss),
                  'coverage', ps['charge_fraction_of_original'], flush=True)
    table(out / 'growth.csv', rows)
    table(out / 'initial_bunching.csv', initial)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout='constrained')
    fig2, ea = plt.subplots(figsize=(7.5, 4.6), layout='constrained')
    for name, color in [('Puffin', '#2679ae'), ('APRL EC+F3', '#d37b28')]:
        rr = [r for r in rows if r['method'] == name]
        bb = [r for r in initial if r['method'] == name] + rr
        axes[0].semilogy([r['actual_z_m'] for r in rr], [r['Hann_pm5pct_power_Eproxy_W'] for r in rr],
                         'o-', color=color, label=name)
        axes[1].semilogy([r['actual_z_m'] for r in bb], [r['central_10lambda_rms_b1'] for r in bb],
                         'o-', color=color, label=name)
        ea.plot([r['actual_z_m'] for r in bb], [r['mean_loss_eV']/1000 for r in bb],
                'o-', color=color, label=name)
    axes[0].set(xlabel='Lab z (m)', ylabel='Central resonant power (W)', title='Hann ±5% band; aperture ±520 µm')
    axes[1].set(xlabel='Lab z (m)', ylabel='Charge-weighted local RMS |b1|', title='Same central window; 10λ bins')
    for ax in axes:
        ax.grid(alpha=.22, which='both'); ax.legend()
    fig.suptitle('Independent random source, w ≈ 355 | six-segment comparison\nSame events, magnets and fixed 0.497 ps gate; no amplitude fitting', fontsize=12)
    fig.savefig(out / 'power_bf_growth.png', dpi=180)
    fig.savefig(out / 'power_bf_growth.pdf'); plt.close(fig)
    ea.set(xlabel='Lab z (m)', ylabel='Mean laboratory energy loss (keV/electron)',
           title='Relative to the same physical input mean energy')
    ea.grid(alpha=.22); ea.legend(); ea.axhline(0, color='gray', lw=.6)
    fig2.savefig(out / 'particle_energy_loss.png', dpi=180)
    fig2.savefig(out / 'particle_energy_loss.pdf'); plt.close(fig2)
    for path, before in protected.items():
        need(identity(path) == before, 'Original closed output changed')
    save(out / 'PROTECTED_OUTPUTS.json', protected)
    save(out / 'RESULT.json', dict(complete=True, physical_acceptance=False, rows=rows, initial=initial,
         closure=closure, source_gamma_mean=str(source_gamma), source_weight=info['actual_weight'],
         input_sha256=info['input_sha256'], fixed_window_relative_to_charge_median_s=[float(low), float(high)],
         fixed_window_duration_s=float(high-low), methods=['Puffin', 'APRL EC+F3'],
         independent_source=True, source_records=plan['sample_records'], wavelength_m=plan['wavelength_m'],
         max_Parseval_relative=max(r['Parseval_relative_error'] for r in rows),
         protected_closed_outputs=len(protected), no_new_physics=True,
         limitations=['PIC E-field band power includes any surviving near-field contribution.',
           'Puffin is forward/paraxial with native filter; its input 1cm drift omits PIC self-field.',
           'PIC uses its stored detector interpolation; Puffin uses native transverse nodes.',
           'Local RMS BF, central global BF and full-beam BF are distinct observables.',
           'Lines join six exit measurements, not a continuous saturation determination.',
           'Particle mean loss is relative to the same imported physical source, not a field energy ledger.']))


if __name__ == '__main__':
    main(*sys.argv[1:])
