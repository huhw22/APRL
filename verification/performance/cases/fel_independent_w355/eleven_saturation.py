"""Closed-only eleven-segment power history and resonant exit statistics."""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if not (HERE / 'growth_compare.py').exists():
    sys.path.insert(0, str(HERE.parent / 'puffin_six_same_ledger'))
    sys.path.insert(0, str(HERE.parent / 'puffin_eleven_four_slip'))

import h5py
import numpy as np
from scipy.signal import find_peaks
from common import C, read, save, need, frozen, dump_at
from growth_compare import particle_stats, field_power, table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD = np.longdouble


def identity(path):
    s = Path(path).stat()
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def native_history(root, out, plan, info, low, high, protected):
    chosen = {}
    for path in (root / 'main').glob('main_integrated_*.h5'):
        with h5py.File(path, 'r') as f:
            step = int(f['runInfo'].attrs['iCsteps'])
        suffix = int(path.stem.rsplit('_', 1)[1])
        if step not in chosen or suffix < chosen[step][0]: chosen[step] = (suffix, path)
    rows = []
    sc = info['scale']
    origin = LD(info['arrival_origin_s']) - LD('.01')/LD(C)
    for step, (_, path) in sorted(chosen.items()):
        protected[str(path)] = identity(path)
        with h5py.File(path, 'r') as f:
            attrs = f['runInfo'].attrs
            power = f['powerSI'][:]
            q = abs(f['Slice Charge SI'][:]).astype(LD)
            bf = f['bunchingFundamentalSI'][:]
            dt = LD(float(attrs['sLengthOfElmZ2'])) * LD(float(attrs['Lc'])) / LD(sc['code_c_m_per_s'])
            ds = LD(sc['lambda_reference_m']) / LD(sc['code_c_m_per_s'])
            z = .01 + float(attrs['zbarTotal']) * float(attrs['Lg'])
        need(np.isfinite(power).all() and np.all(power >= 0) and np.isfinite(q).all()
             and np.isfinite(bf).all() and q.sum() > 0, 'Invalid native profile')
        cumulative = np.cumsum(q, dtype=LD)
        j = int(np.searchsorted(cumulative, cumulative[-1]/2))
        previous = cumulative[j-1] if j else LD(0)
        # Native slices do not retain individual arrival times. Uniform-in-bin
        # interpolation defines a continuous charge median at 1-lambda resolution.
        median = origin + (LD(j) + (cumulative[-1]/2-previous)/q[j]) * ds
        time = origin + np.arange(len(power), dtype=LD) * dt
        rel = time - median
        use = (rel >= low) & (rel < high)
        need(rel[0] <= low and rel[-1] >= high and use.any(), 'Native central window outside record')
        weight = np.sqrt(8/3)*np.sin(np.pi*np.asarray((rel[use]-low)/(high-low), float))**2
        btime = origin + (np.arange(len(q), dtype=LD)+LD('.5'))*ds
        central = (btime-median >= low) & (btime-median < high) & (q > 0)
        rows.append(dict(step=step, lab_z_m=z, charge_median_retarded_time_s=str(median),
            central_Hann_native_allgrid_allband_power_W=float(np.sum(power[use]*weight**2, dtype=LD)*dt/(high-low)),
            central_Rect_native_allgrid_allband_power_W=float(np.sum(power[use], dtype=LD)*dt/(high-low)),
            full_record_retained_field_energy_J=float(np.sum(power, dtype=LD)*dt),
            central_native_1lambda_rms_b1=float(np.sqrt(np.sum(q[central]*bf[central]**2, dtype=LD)/q[central].sum())),
            slice_charge_C=float(q.sum()), source=str(path)))
    need(rows[0]['step'] == 0 and rows[-1]['step'] == 79200, 'Native history is incomplete')
    table(out / 'native_history.csv', rows)
    regular = [r for r in rows if r['step'] % 80 == 0]
    need([r['step'] for r in regular] == list(range(0, 79201, 80)), 'Missing per-period native output')
    power = np.array([r['central_Hann_native_allgrid_allband_power_W'] for r in regular])
    ids, props = find_peaks(power, prominence=.1*power.max(), distance=math.ceil(.5/.028))
    candidates = [dict(lab_z_m=regular[i]['lab_z_m'], step=regular[i]['step'],
                       power_W=float(power[i]), prominence_W=float(props['prominences'][j]))
                  for j, i in enumerate(ids)]
    save(out / 'NATIVE_PEAK_CANDIDATES.json', dict(
        definition='Local maxima of per-period full-grid retained-frequency power in the common Hann gate; prominence >=10% of history maximum, separation >=0.5m.',
        candidates=candidates, quantity_is_not_aperture_band_power=True,
        native_charge_median_resolution_s=float(ds), relative_gate_s=[float(low),float(high)]))
    print('NATIVE_PEAK_CANDIDATES', candidates, flush=True)
    fig, ax = plt.subplots(figsize=(10, 4.6), layout='constrained')
    ax.plot([r['lab_z_m'] for r in rows], [r['central_Hann_native_allgrid_allband_power_W']/1e6 for r in rows])
    for peak in candidates:
        ax.plot(peak['lab_z_m'], peak['power_W']/1e6, 'o')
        ax.annotate(f"{peak['lab_z_m']:.2f}m", (peak['lab_z_m'], peak['power_W']/1e6), xytext=(0, 8), textcoords='offset points')
    ax.set(xlabel='Lab z (m)', ylabel='Native central Hann power (MW)', title='Full transverse grid, all retained frequencies; per-period output')
    ax.grid(alpha=.2); fig.savefig(out/'native_preview.png', dpi=160); plt.close(fig)
    return rows, candidates


def main(root, out):
    root, out = Path(root).resolve(), Path(out).resolve()
    need(root.name == 'eleven-single-node-r1', 'Use the current single-node eleven-segment case')
    need(read(root/'TRACKING_CLOSED.json')['complete'] and read(root/'COMPLETE.json')['complete'], 'Unclosed Puffin')
    need(read(root/'main/puffin.log.launch.json')['exit_code'] == 0, 'Failed physical process')
    need('Puffin simulation completed' in (root/'main/puffin.log').read_text(), 'Missing normal completion')
    frozen(root/'FROZEN_PHYSICS.json')
    plan, info = read(root/'PLAN.json'), read(root/'INPUT_VALIDATED.json')
    need(plan['segments'] == 11 and info['independent_fixed_N_source'] and info['no_clones'], 'Wrong source')
    out.mkdir(parents=True, exist_ok=False)
    half = LD(str(plan['duration_s']))/2
    low, high = LD('.1')*LD(str(math.sqrt(.2)))*half, LD('.9')*LD(str(math.sqrt(.2)))*half
    protected = {}
    dense, candidates = native_history(root, out, plan, info, low, high, protected)
    rows = []
    for item in plan['outputs']:
        particle = dump_at(root/'main', 'electrons', item['step'])
        field = dump_at(root/'main', 'aperp', item['step'])
        for path in [particle, field]: protected[str(path)] = identity(path)
        ps = particle_stats(particle, 'Puffin', info, plan, low, high)
        need(ps['records'] == plan['sample_records'] and abs(ps['charge_fraction_of_original']-1) < 2e-12,
             'Particle population or charge changed')
        # dump_at requires the exact integration step. Iterated floating zbar
        # differs from the analytically computed lab plane by tens of pm after
        # 57k+ steps; 1nm remains far below the 350um integration displacement.
        z_roundoff = ps['actual_z_m']-item['actual_lab_z_m']
        need(abs(z_roundoff) < 1e-9, 'Wrong physical exit')
        fp = field_power(field, 'Puffin', info, plan, LD(ps['median_retarded_time_s']), low, high)
        row = dict(segment=item['segment'], plane_z_roundoff_difference_m=z_roundoff, **ps, **fp)
        rows.append(row)
        table(out/'gain.csv', rows)
        print('EXIT', item['segment'], 'z', ps['actual_z_m'], 'band_power_W', fp['Hann_pm5pct_power_Eproxy_W'],
              'BF', ps['central_10lambda_rms_b1'], flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), layout='constrained')
    z = [r['lab_z_m'] for r in dense]
    power = [r['central_Hann_native_allgrid_allband_power_W']/1e6 for r in dense]
    for ax in axes:
        ax.plot(z, power, color='#7b94aa', lw=1.1, label='Native all-grid / all retained frequencies')
        ax.plot([r['actual_z_m'] for r in rows], [r['Hann_pm5pct_power_Eproxy_W']/1e6 for r in rows],
                'o-', color='#c8792e', label='Hann ±5% band / aperture ±520 µm')
        for item in plan['outputs']: ax.axvline(item['actual_lab_z_m'], color='gray', lw=.5, alpha=.22)
        ax.set(xlabel='Lab z (m)', ylabel='Central power (MW)')
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set_yscale('log'); axes[0].set_ylim(bottom=1e-5)
    axes[0].set_title('Eleven-segment growth')
    axes[1].set_xlim(10, 28); axes[1].set_ylim(bottom=0)
    axes[1].set_title('Nonlinear region; native peaks are a separate observable')
    for peak in candidates:
        if peak['lab_z_m'] >= 10:
            axes[1].annotate(f"{peak['lab_z_m']:.2f}m", (peak['lab_z_m'],peak['power_W']/1e6),
                             xytext=(0,8), textcoords='offset points', ha='center', fontsize=8)
    fig.suptitle('Puffin | independent random source, w ≈ 355 | four-slip bunch\nSame central-window selection rule; fixed gate 0.993 ps, no amplitude fitting', fontsize=12)
    fig.savefig(out/'power_growth.png', dpi=180); fig.savefig(out/'power_growth.pdf'); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.4), layout='constrained')
    ax.plot([r['actual_z_m'] for r in rows], [r['central_10lambda_rms_b1'] for r in rows], 'o-')
    ax.set(xlabel='Lab z (m)', ylabel='Central charge-weighted local RMS |b1| (10λ)', title='Eleven exit planes; same central window')
    ax.grid(alpha=.2); fig.savefig(out/'bunching_growth.png', dpi=180); fig.savefig(out/'bunching_growth.pdf'); plt.close(fig)
    for path, before in protected.items(): need(identity(path) == before, 'Closed original output changed')
    frozen(root/'FROZEN_PHYSICS.json')
    save(out/'PROTECTED_OUTPUTS.json', protected)
    save(out/'RESULT.json', dict(complete=True, physical_acceptance=False, segments=rows,
         native_history=dense, native_peak_candidates=candidates,
         fixed_window_relative_to_charge_median_s=[float(low),float(high)], fixed_window_duration_s=float(high-low),
         source_weight=info['actual_weight'], source_records=plan['sample_records'], input_sha256=info['input_sha256'],
         closed_original_outputs=len(protected), max_Parseval_relative=max(r['Parseval_relative_error'] for r in rows),
         no_physics_rerun=True, physical_job=19799,
         limits=['Native per-period power is full-grid and all retained frequencies, not the ±5% aperture quantity.',
                 'Native charge median is interpolated within one optical-wavelength slice; exit medians use actual particles.',
                 'Native local BF is 1λ, while plotted exit BF is the same10λ RMS definition as the six-segment review.',
                 'No saved intermediate 3D field allows exact dense ±5% aperture-band peak localization.',
                 'First strong peak and later oscillations are distinct from a unique steady-state saturation point.',
                 'One weighted source is not a one-electron SASE noise acceptance; four slips do not guarantee a flat pulse.']))


if __name__ == '__main__':
    main(*sys.argv[1:])
