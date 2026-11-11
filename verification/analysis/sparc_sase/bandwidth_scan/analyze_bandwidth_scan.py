#!/usr/bin/env python3
"""Read completed detector planes once; scan frequency width at fixed time windows."""
import argparse
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
from scipy.constants import c, epsilon_0, h, e, mu_0

import analyze_stable_half as base
import analyze_central_power as central

WIDTHS = (.01, .02, .05, .10, .20, .50)
FLOORS = (.6, .7, .8, .9)
MAIN_FLOOR = .8


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist())+'\n')


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def stat_identity(path):
    st = path.stat()
    return dict(path=str(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns, inode=st.st_ino)


def window_masks(time_s, interval):
    u = (time_s-interval['t_start_s'])/interval['half_duration_s']
    masks = {}
    for floor in FLOORS:
        edge = math.sqrt(1-floor)
        low, high = .1*edge, .9*edge
        masks[f'current_{floor:.1f}'] = (u >= low) & (u < high)
    for guard in (0, .05, .10, .15, .20):
        masks[f'half_guard_{guard:.2f}'] = (u >= guard) & (u < 1-guard)
    for name, mask in masks.items():
        base.need(mask.sum() > 100, 'empty interval: '+name)
    return u, masks


def spectral_power(transform, dt, n, area):
    one = np.full(transform.shape[-1], 2.)
    one[0] = 1.
    if n % 2 == 0:
        one[-1] = 1.
    power_e = np.sum(abs(transform[:, :2])**2, axis=(0, 1))
    cross = (transform[:, 0]*transform[:, 3].conj() -
             transform[:, 1]*transform[:, 2].conj()).real
    power_s = np.sum(cross, axis=0)
    return (power_e*one*dt/n*area*epsilon_0*c,
            power_s*one*dt/n*area/mu_0)


def envelope_power(transform, mask, n, area, workers=8):
    # Positive frequencies doubled -> analytic signal.  Half its modulus
    # squared is the optical-cycle mean, with no amplitude/window rescaling.
    coefficients = np.zeros((*transform.shape[:-1], n), dtype=np.complex128)
    coefficients[:, :, np.flatnonzero(mask)] = 2*transform[:, :, mask]
    analytic = fft.ifft(coefficients, axis=-1, workers=workers, overwrite_x=True)
    ep = .5*epsilon_0*c*area*np.sum(
        abs(analytic[:, 0])**2+abs(analytic[:, 1])**2, axis=0)
    sp = area/(2*mu_0)*np.sum(
        (analytic[:, 0]*analytic[:, 3].conj()-
         analytic[:, 1]*analytic[:, 2].conj()).real, axis=0)
    return ep, sp


def self_test():
    checks = []
    rng = np.random.default_rng(9725)
    for n in (511, 512):
        values = rng.normal(size=(3, 4, n))
        transformed = fft.rfft(values, axis=-1)
        dt, area = .01, .7
        es, ss = spectral_power(transformed, dt, n, area)
        expected_e = epsilon_0*c*area*dt*np.sum(values[:, :2]**2)
        expected_s = area/mu_0*dt*np.sum(
            values[:, 0]*values[:, 3]-values[:, 1]*values[:, 2])
        checks.append(abs(es.sum()/expected_e-1) < 1e-12)
        checks.append(abs(ss.sum()-expected_s) < 1e-12*abs(expected_s))
        mask = np.zeros(transformed.shape[-1], bool)
        mask[20:100] = True
        ep, sp = envelope_power(transformed, mask, n, area, workers=1)
        checks.append(abs(ep.sum()*dt/es[mask].sum()-1) < 1e-12)
        checks.append(abs(sp.sum()*dt-ss[mask].sum()) <
                      1e-12*max(abs(ss[mask].sum()), 1.))
    # A second tone is outside a narrow band and inside a wider one.  Its
    # independent factor-four energy must be recovered after widening.
    n = 1024
    values = np.zeros((1, 4, n))
    t = np.arange(n)/n
    values[0, 0] = np.cos(2*np.pi*100*t)+2*np.cos(2*np.pi*114*t)
    values[0, 3] = values[0, 0]/c
    transformed = fft.rfft(values, axis=-1)
    a = np.arange(transformed.shape[-1])
    narrow = abs(a/100-1) <= .05
    wide = abs(a/100-1) <= .20
    p1, _ = envelope_power(transformed, narrow, n, 1, workers=1)
    p2, _ = envelope_power(transformed, wide, n, 1, workers=1)
    checks.append(abs(p2.mean()/p1.mean()-5) < 1e-12)
    base.need(all(checks), 'FFT normalization/bandwidth self-test failed')
    return dict(passed=len(checks), odd_even_parseval=True,
                envelope_spectral_energy=True, two_tone_width_recovery=True)


def analyze_plane(path, interval, f0, output, module):
    ts, dt, area, cube, native, missing = base.detector_cube(path, 6)
    n, pixels, _ = cube.shape
    frequency = fft.rfftfreq(n, dt)
    masks = [abs(frequency/f0-1) <= width for width in WIDTHS]
    u, tmasks = window_masks(ts, interval)
    core = tmasks['current_0.8']
    m = int(core.sum())
    core_freq = fft.rfftfreq(m, dt)
    window_full = np.hanning(n)
    window_full /= np.sqrt(np.mean(window_full**2))
    window_core = np.hanning(m)
    window_core /= np.sqrt(np.mean(window_core**2))
    full_e, full_s, hann_e, hann_s = [np.zeros(len(frequency)) for _ in range(4)]
    core_e, core_s = [np.zeros(len(core_freq)) for _ in range(2)]
    envelope_e, envelope_s = [np.zeros((len(WIDTHS), n)) for _ in range(2)]
    direct_e = 0.
    direct_s = 0.
    for start in range(0, pixels, 32):
        values = np.array(cube[:, start:start+32, :][:, :, [0, 1, 3, 4]].transpose(1, 2, 0),
                          order='C', copy=True)
        direct_e += epsilon_0*c*area*dt*np.sum(values[:, :2]**2)
        direct_s += area/mu_0*dt*np.sum(
            values[:, 0]*values[:, 3]-values[:, 1]*values[:, 2])
        transformed = fft.rfft(values, axis=-1, workers=8)
        de, ds = spectral_power(transformed, dt, n, area)
        full_e += de
        full_s += ds
        ht = fft.rfft(values*window_full[None, None, :], axis=-1, workers=8)
        de, ds = spectral_power(ht, dt, n, area)
        hann_e += de
        hann_s += ds
        ct = fft.rfft(values[:, :, core]*window_core[None, None, :], axis=-1, workers=8)
        de, ds = spectral_power(ct, dt, m, area)
        core_e += de
        core_s += ds
        del ht, ct
        for index, mask in enumerate(masks):
            ep, sp = envelope_power(transformed, mask, n, area)
            envelope_e[index] += ep
            envelope_s[index] += sp
        del values, transformed
    del cube
    spectrum_file = output/f'p{module:02d}_spectra.npz'
    envelope_file = output/f'p{module:02d}_envelopes.npz'
    np.savez_compressed(spectrum_file, frequency_Hz=frequency,
        photon_energy_eV=frequency*h/e, rect_E_bin_J=full_e, rect_flux_bin_J=full_s,
        hann_E_bin_J=hann_e, hann_flux_bin_J=hann_s,
        core_frequency_Hz=core_freq, core_hann_E_bin_J=core_e,
        core_hann_flux_bin_J=core_s)
    # Small reduced time series permit new time cuts without another 25 GB read.
    np.savez_compressed(envelope_file, widths=WIDTHS, u=u,
        time_offset_from_midpoint_s=ts-interval['t_start_s'], dt_s=dt,
        power_E_W=envelope_e, power_flux_W=envelope_s)
    rows, window_rows, errors = [], [], []
    for index, (width, mask) in enumerate(zip(WIDTHS, masks)):
        ce_mask = abs(core_freq/f0-1) <= width
        ep, sp = envelope_e[index], envelope_s[index]
        energy_e, energy_s = full_e[mask].sum(), full_s[mask].sum()
        envelope_error = abs(ep.sum()*dt/energy_e-1)
        errors.append(envelope_error)
        row = dict(module=module, relative_half_width=width,
            rect_pulse_E_J=float(energy_e), rect_pulse_signed_flux_J=float(energy_s),
            hann_pulse_E_J=float(hann_e[mask].sum()),
            hann_pulse_signed_flux_J=float(hann_s[mask].sum()),
            central_envelope_mean_W=float(ep[core].mean()),
            central_envelope_median_W=float(np.median(ep[core])),
            central_envelope_p90_W=float(np.quantile(ep[core], .9)),
            central_envelope_peak_W=float(ep[core].max()),
            central_signed_flux_mean_W=float(sp[core].mean()),
            central_hann_band_J=float(core_e[ce_mask].sum()),
            central_hann_mean_W=float(core_e[ce_mask].sum()/(m*dt)),
            central_duration_s=m*dt, full_band_bins=int(mask.sum()),
            central_band_bins=int(ce_mask.sum()),
            envelope_vs_spectral_energy_relative_error=float(envelope_error))
        rows.append(row)
        for name, tmask in tmasks.items():
            window_rows.append(dict(module=module, relative_half_width=width,
                time_window=name, samples=int(tmask.sum()), duration_s=tmask.sum()*dt,
                mean_W=float(ep[tmask].mean()), median_W=float(np.median(ep[tmask])),
                peak_W=float(ep[tmask].max()),
                signed_flux_mean_W=float(sp[tmask].mean())))
    region = (frequency/f0 >= .5) & (frequency/f0 <= 1.5)
    region_indices = np.flatnonzero(region)
    peak = region_indices[np.argmax(hann_e[region])]
    cumulative = np.cumsum(hann_e[region])
    quantiles = np.interp([.025, .5, .975], cumulative/cumulative[-1],
                          frequency[region]/f0)
    checks = dict(module=module, samples=n, native_samples=native,
        missing_samples=missing, dt_s=dt,
        parseval_E_relative_error=float(abs(full_e.sum()/direct_e-1)),
        parseval_flux_relative_error=float(abs(full_s.sum()-direct_s)/max(abs(direct_s), 1e-300)),
        envelope_parseval_max_error=float(max(errors)),
        hann_peak_frequency_over_f0=float(frequency[peak]/f0),
        hann_main_region_95_percent_frequency_over_f0=quantiles.tolist())
    base.need(checks['parseval_E_relative_error'] < 1e-10 and max(errors) < 1e-10,
              'Parseval failed on physical data')
    return rows, window_rows, checks


def render(output, rows, checks, xie_power):
    colors = ['#777777', '#4d83ad', '#24648c', '#cf8a24', '#8e4b75', '#75884b']
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'lines.linewidth': 1.7})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for width, color in zip(WIDTHS, colors):
        current = [row for row in rows if row['relative_half_width'] == width]
        z = [row['z_m'] for row in current]
        label = f'±{100*width:g}%'
        axes[0].semilogy(z, [row['central_envelope_mean_W']/1e6 for row in current],
                         'o-', color=color, label=label, markersize=4)
        axes[1].semilogy(z, [row['hann_pulse_E_J']*1e6 for row in current],
                         'o-', color=color, label=label, markersize=4)
    axes[0].axhline(xie_power/1e6, color='#333333', linestyle='--',
                    label=f'Xie: {xie_power/1e6:.2f} MW')
    axes[0].set_ylabel('central mean band power [MW]')
    axes[0].set_title('(a) Same central time window for every width')
    axes[1].set_ylabel('full-record Hann band energy [µJ]')
    axes[1].set_title('(b) Pulse energy and frequency-width dependence')
    for ax in axes:
        ax.set_xlabel('undulator-exit z [m]')
        ax.grid(which='both', alpha=.2)
        ax.legend(fontsize=9)
    fig.savefig(output/'bandwidth_growth.png', dpi=190)
    fig.savefig(output/'bandwidth_growth.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for module, color in [(5, '#24648c'), (6, '#8e4b75'), (8, '#cf8a24'), (11, '#75884b')]:
        with np.load(output/f'p{module:02d}_envelopes.npz') as data:
            u = data['u']
            use = (u >= -.05) & (u <= 1.02)
            for width, ls in [(.05, '-'), (.20, '--'), (.50, ':')]:
                index = list(WIDTHS).index(width)
                axes[0].plot(u[use], data['power_E_W'][index, use]/1e6,
                    color=color, linestyle=ls, label=f'p{module}, ±{100*width:g}%')
        with np.load(output/f'p{module:02d}_spectra.npz') as data:
            energy = data['photon_energy_eV']
            f0energy = float(checks[module-1]['photon_energy_eV'])
            x = energy/f0energy
            use = (x >= .5) & (x <= 1.5)
            density = data['hann_E_bin_J']/(energy[1]-energy[0])
            axes[1].semilogy(x[use], density[use], color=color, label=f'p{module}')
    core_edge = math.sqrt(1-MAIN_FLOOR)
    axes[0].axvspan(.1*core_edge, .9*core_edge, color='#777777', alpha=.12,
                    label='central statistic interval')
    axes[0].set_xlabel('position in later bunch half (0=centre, 1=tail)')
    axes[0].set_ylabel('band-power envelope [MW]')
    axes[0].set_title('(a) Envelopes at selected exits')
    axes[0].legend(fontsize=8, ncol=2)
    axes[1].axvline(.95, color='#555555', linestyle=':', linewidth=.8)
    axes[1].axvline(1.05, color='#555555', linestyle=':', linewidth=.8)
    axes[1].set_xlabel('photon energy / nominal resonance energy')
    axes[1].set_ylabel('full-record Hann spectrum [J/eV]')
    axes[1].set_title('(b) Nominal ±5% edges and spectral movement')
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(output/'envelopes_and_spectra.png', dpi=190)
    fig.savefig(output/'envelopes_and_spectra.pdf')
    plt.close(fig)


def main(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    base.need(not output.exists(), 'output already exists')
    marker = json.loads((root/'energy-binomial3/PIC_COMPLETE.json').read_text())
    base.need(marker['complete'] and marker['run']['exit_code'] == 0, 'case incomplete')
    plan = json.loads((root/'plan.json').read_text())
    prior = json.loads((root/'analysis-energy-binomial3/RESULT.json').read_text())['rows']
    oldhalf = json.loads((root/'analysis-energy-binomial3/stable-half/RESULT.json').read_text())['rows']
    physics = plan['physics']
    wavelength = float(physics['wavelength_m'])
    f0 = c/wavelength
    output.mkdir(parents=True)
    save(output/'SELF_TEST.json', self_test())
    before, after, rows, wrows, checks = [], [], [], [], []
    regressions = []
    started = time.monotonic()
    for module, (plane, z) in enumerate(plan['geometry']['field_planes'], start=1):
        detector = root/'energy-binomial3/detectors'/f'{plane}.h5'
        before.append(stat_identity(detector))
        interval = base.particle_interval_and_bunching(
            root/'energy-binomial3/detectors'/f'{plane}-particles.h5', module, float(z),
            int(plan['NPART']), int(plan['active_slices']), wavelength, 0.)
        newrows, windows, check = analyze_plane(detector, interval, f0, output, module)
        for row in newrows:
            row['z_m'] = float(z)
        for row in windows:
            row['z_m'] = float(z)
        check.update(z_m=float(z), photon_energy_eV=h*f0/e, interval=interval)
        rows.extend(newrows)
        wrows.extend(windows)
        checks.append(check)
        after.append(stat_identity(detector))
        five = next(row for row in newrows if row['relative_half_width'] == .05)
        old = prior[module-1]['hann_E_proxy_band_J']
        legacy = next(row for row in windows if row['relative_half_width'] == .05 and
                      row['time_window'] == 'half_guard_0.10')
        oldmean = oldhalf[module-1]['main']['mean_E_proxy_W']
        regression = dict(module=module, hann_energy_relative_error=five['hann_pulse_E_J']/old-1,
                          legacy_half_mean_relative_error=legacy['mean_W']/oldmean-1)
        base.need(abs(regression['hann_energy_relative_error']) < 1e-10 and
                  abs(regression['legacy_half_mean_relative_error']) < 1e-10,
                  'previous ±5% regression failed')
        base.need(all(np.diff([row['hann_pulse_E_J'] for row in newrows]) >= -1e-20),
                  'nonmonotone integrated band energy')
        regressions.append(regression)
        save(output/f'p{module:02d}_RESULT.json', dict(rows=newrows, windows=windows,
                                                    check=check, regression=regression))
        print(json.dumps(dict(module=module, elapsed_s=time.monotonic()-started,
            means_MW={f'{100*row["relative_half_width"]:g}%':
                      row['central_envelope_mean_W']/1e6 for row in newrows},
            hann_energy_uJ={f'{100*row["relative_half_width"]:g}%':
                           row['hann_pulse_E_J']*1e6 for row in newrows})), flush=True)
    base.need(before == after, 'source files changed during analysis')
    write_csv(output/'bandwidth_scan.csv', rows)
    write_csv(output/'time_window_sensitivity.csv', wrows)
    l1d = float(physics['period_m'])/(4*math.pi*math.sqrt(3)*.00288)
    beam = physics['peak_current_A']*physics['kinetic_energy_MeV']*1e6
    xie = 1.6*.00288*(l1d/.67)**2*beam
    summary = []
    for width in WIDTHS:
        current = [row for row in rows if row['relative_half_width'] == width]
        peak = max(current, key=lambda row: row['central_envelope_mean_W'])
        summary.append(dict(relative_half_width=width, peak_module=peak['module'],
            peak_z_m=peak['z_m'], peak_central_mean_W=peak['central_envelope_mean_W'],
            peak_mean_over_xie=peak['central_envelope_mean_W']/xie,
            final_central_mean_W=current[-1]['central_envelope_mean_W'],
            final_hann_pulse_energy_J=current[-1]['hann_pulse_E_J']))
    render(output, rows, checks, xie)
    save(output/'RESULT.json', dict(complete=True, case='energy-binomial3',
        widths=WIDTHS, band_coordinate='photon energy or frequency, not wavelength',
        central_time_window='later half, parabolic I/Ipk>=0.8, then 10% guards on that interval',
        current_floors=FLOORS, fixed_time_window_across_widths=True,
        xie_saturation_power_W=xie, summary=summary, rows=rows, plane_checks=checks,
        limitations=['one physical noise seed', 'raw total E/B, native B time stagger',
            'all frequencies within 0.5..1.5 resonance are a broad-band diagnostic',
            'local filtered power can vary nonmonotonically with width due to interference',
            'full-record band energy is additive and must increase with width']))
    save(output/'VALIDATION.json', dict(complete=True, source_metadata_unchanged=before==after,
        parseval_E_max_error=max(row['parseval_E_relative_error'] for row in checks),
        envelope_parseval_max_error=max(row['envelope_parseval_max_error'] for row in checks),
        previous_five_percent_regressions=regressions, source_files=before,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_s=time.monotonic()-started))
    (output/'README.md').write_text('# 中心区功率与频带宽度扫描\n\n'
        '共振频率两侧分别取±1%、±2%、±5%、±10%、±20%、±50%。所有宽度共享时间窗口。'
        '主窗口是后半束团中抛物线电流I/Ipk≥0.8的中心段，再从该段两端各裁10%；'
        '另保留0.6–0.9电流阈值及旧后半束团0–20%裁边。\n\n'
        'bandwidth_scan.csv保留完整记录Hann/矩形频带能量、中心时域功率及'
        '中心段自身Hann频谱三种口径。包络未用全记录Hann乘幅；'
        '谱能量沿用原RMS Hann窗定义。envelopes_and_spectra.png展示频谱移动。\n\n'
        '每面缩减包络和频谱缓存允许随后只改统计窗口，不必重新读取大型H5。'
        '完整频带积分随频带加宽应单调上升；局部时域功率包含频率间干涉，'
        '不要求逐点单调。原生总场与单噪声种子限制保持。\n')
    save(output/'COMPLETE.json', dict(complete=True, planes=11, widths=len(WIDTHS)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root', nargs='?')
    parser.add_argument('output', nargs='?')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test()))
    else:
        main(args.root, args.output)
