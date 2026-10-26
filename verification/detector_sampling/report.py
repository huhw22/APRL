#!/usr/bin/env python3
"""Full-plane spectral checks for the analytic MPI detector fixture (no filtering)."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def load(directory):
    meta = json.loads((directory / 'sweep-metadata.json').read_text())
    shape = tuple(meta['shape'])
    with h5py.File(directory / 'sweep.h5') as file:
        group = file['field_plane']
        assert int(group['complete'][()]) == 1
        fields = group['fields'][:]
        actual = np.concatenate((fields['electric_V_per_m'],
                                 fields['magnetic_T'] * 299792458.0), axis=-1)
        assert actual.shape == shape and np.isfinite(actual).all()
        assert int(group['committed_samples'][()]) == shape[0]
    reference = np.fromfile(directory / 'sweep-reference.bin', dtype=np.float64).reshape(shape)
    legacy = np.fromfile(directory / 'sweep-legacy.bin', dtype=np.float64).reshape(shape)
    return meta, actual, reference, legacy


def metrics(actual, reference, tone):
    # FFT the full, unmodified trace at every pixel. No detrend, clipping or window.
    fa = np.fft.rfft(actual, axis=0) / len(actual)
    fr = np.fft.rfft(reference, axis=0) / len(reference)
    a, r = fa[tone], fr[tone]
    mask = np.abs(r) > np.max(np.abs(r)) * 1e-8
    complex_error = np.max(np.abs(a[mask] - r[mask]) / np.abs(r[mask]))
    amplitude_error = np.max(np.abs(np.abs(a[mask]) / np.abs(r[mask]) - 1))
    phase_error = np.max(np.abs(np.angle(a[mask] * r[mask].conj())))
    # All ordered pixel pairs for Ex, including their complex cross spectrum.
    ax, rx = a[..., 0].ravel(), r[..., 0].ravel()
    ca, cr = ax[:, None] * ax.conj()[None, :], rx[:, None] * rx.conj()[None, :]
    csd_error = np.max(np.abs(ca-cr)) / np.max(np.abs(cr))
    component_scale = np.maximum(np.max(np.abs(reference), axis=(0, 1, 2)), 1.)
    error_fft = np.fft.rfft(actual-reference, axis=0) / len(actual)
    return dict(max_scaled_point_error=float(np.max(np.abs(actual-reference)/component_scale)),
                max_tone_complex_relative_error=float(complex_error),
                max_tone_amplitude_relative_error=float(amplitude_error),
                max_tone_phase_error_rad=float(phase_error),
                ex_all_pixel_pairs_csd_relative_error=float(csd_error),
                ex_max_non_dc_error_fft_amplitude=float(np.max(np.abs(error_fft[1:, ..., 0]))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--control', type=Path, help='Actual pre-fix executable output on same layout')
    args = parser.parse_args()
    meta, actual, reference, legacy = load(args.directory)
    tone = meta['tone_bin']
    corrected, old = metrics(actual, reference, tone), metrics(legacy, reference, tone)
    report = dict(scope=meta['fixture'], fft='raw rectangular, no detrending or filtering',
                  ranks=meta['ranks'], samples=len(actual), full_plane_pixels=int(np.prod(actual.shape[1:3])),
                  tone_bin=tone, corrected=corrected, legacy_local_sampler=old)
    if args.control:
        cm, control, cref, _ = load(args.control)
        assert cm == meta and np.array_equal(cref, reference), 'Control input differs'
        report['actual_pre_fix_executable'] = metrics(control, reference, tone)
        report['pre_fix_matches_legacy_max_absolute'] = float(np.max(np.abs(control-legacy)))
        assert np.max(np.abs(control-legacy)) < 1e-9 * np.max(np.abs(legacy))
        legacy = control
    report['pass'] = bool(corrected['max_scaled_point_error'] < 2e-10 and
                          corrected['max_tone_complex_relative_error'] < 2e-7 and
                          corrected['max_tone_phase_error_rad'] < 2e-7 and
                          corrected['ex_all_pixel_pairs_csd_relative_error'] < 4e-7 and
                          old['max_scaled_point_error'] > 1e-5)
    (args.directory / 'spectral-review.json').write_text(json.dumps(report, indent=2)+'\n')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    center = (actual.shape[1]//2, actual.shape[2]//2)
    ind = (slice(None), *center, 0)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    n = np.arange(len(actual))
    axes[0].plot(n, legacy[ind]-reference[ind], color='#b14b37', lw=1, label='Old sampler minus unsplit reference')
    axes[0].plot(n, actual[ind]-reference[ind], color='#087c8c', lw=1, label='Fixed sampler minus unsplit reference')
    axes[0].set(xlabel='Saved sample', ylabel='Lab Ex sampling error (V/m)',
                title='Boundary interpolation error in a strong-background + weak-tone fixture')
    axes[0].legend(loc='upper right')
    for values, label, color, style in [(legacy, 'Old sampler', '#b14b37', '-'),
                                      (actual, 'Fixed sampler', '#087c8c', '-'),
                                      (reference, 'Unsplit reference', '#222222', '--')]:
        transform = np.fft.rfft(values[..., 0], axis=0)/len(values)
        power = np.mean(np.abs(transform)**2, axis=(1, 2))
        axes[1].semilogy(np.arange(len(power))[1:], np.maximum(power[1:], 1e-30),
                         style, color=color, lw=1, label=label)
    axes[1].axvline(tone, color='#222222', alpha=.4, lw=.8)
    axes[1].set(xlabel='Fourier bin (cycles per saved trace)', ylabel='Mean-pixel |FFT(Ex)/N|² [(V/m)²]',
                title='Raw spectra: known weak tone at bin %d; DC outside displayed range' % tone)
    axes[1].legend(loc='upper right')
    fig.suptitle('MPI detector sampling correction — %d ranks, %d pixels, γ = %.2f' %
                 (meta['ranks'], report['full_plane_pixels'], meta['gamma']))
    fig.savefig(args.directory / 'sampling-before-after.png', dpi=180)
    fig.savefig(args.directory / 'sampling-before-after.pdf')
    plt.close(fig)
    print(json.dumps(report, indent=2))
    if not report['pass']:
        raise SystemExit('Spectral acceptance failed')


if __name__ == '__main__':
    main()
