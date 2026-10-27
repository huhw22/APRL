#!/usr/bin/env python3
"""Independently review the small window probes; never launch computation."""
import argparse
import json
from pathlib import Path

import numpy as np


AXES = ('Ex', 'Ey', 'Ez', 'cBx', 'cBy', 'cBz')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def plane_review(root, case, candidate, exact):
    table = np.genfromtxt(root / case['csv'], delimiter=',', names=True)
    require(all(np.isfinite(table[name]).all() for name in table.dtype.names), 'Nonfinite CSV')
    frames = case['steps'] + 1
    require(len(table) == frames * 30, 'Plane is not the complete 30-pixel sequence')
    require(np.array_equal(table['step'], np.repeat(np.arange(frames), 30)), 'Bad step order')
    def fields(prefix):
        return np.stack([table[prefix + '_' + axis] for axis in AXES], axis=-1).reshape(frames, 30, 6)
    ref, got = fields('reference'), fields(candidate)
    wr = ref - fields('reference_background')
    wg = got - fields(candidate + '_background')
    scale = np.max(np.abs(wr))
    require(scale > 0, 'Weak perturbation absent')
    raw = float(np.max(np.abs(got-ref)) / scale)
    weak = float(np.max(np.abs(wg-wr)) / scale)
    negative = float(np.max(np.abs(fields('zero_boundary')-ref)) / scale)
    require(negative > 100, 'Missing-boundary control does not discriminate')
    if exact:
        require(np.array_equal(ref.view(np.uint64), got.view(np.uint64)), 'Raw plane bits differ')
        require(np.array_equal(wr.view(np.uint64), wg.view(np.uint64)), 'Weak plane bits differ')
    else:
        require(raw < 1e-5 and weak < 1e-5, 'Analytic plane tolerance exceeded')
    # Preserve every frequency and every pixel/polarization. No windowing,
    # demeaning, frequency deletion or intensity-only comparison is used.
    fr, fg = np.fft.rfft(wr, axis=0), np.fft.rfft(wg, axis=0)
    fft_relative = float(np.linalg.norm(fg-fr) / np.linalg.norm(fr))
    mask = np.abs(fr) > .01 * np.max(np.abs(fr))
    phase = float(np.max(np.abs(np.angle(fg[mask] * np.conj(fr[mask])))))
    error_norm, reference_norm = 0., 0.
    max_csd_error = 0.
    for r, g in zip(fr.reshape(len(fr), -1), fg.reshape(len(fg), -1)):
        rr = r[:, None] * np.conj(r[None, :])
        gg = g[:, None] * np.conj(g[None, :])
        error_norm += float(np.sum(np.abs(gg-rr)**2))
        reference_norm += float(np.sum(np.abs(rr)**2))
        max_csd_error = max(max_csd_error, float(np.max(np.abs(gg-rr))))
    csd_relative = float(np.sqrt(error_norm/reference_norm))
    if exact:
        require(fft_relative == 0 and max_csd_error == 0, 'FFT/cross spectra differ')
    else:
        require(fft_relative < 1e-5 and csd_relative < 3e-5 and phase < 1e-4,
                'Independent weak-field spectral check failed')
    return dict(csv=case['csv'], steps=case['steps'], beta=case['beta'],
                scalar_plane_values=ref.size, raw_error_over_weak=raw,
                weak_error_over_weak=weak, negative_error_over_weak=negative,
                weak_fft_relative_l2=fft_relative,
                phase_error_rad_above_1pct_spectral_peak=phase,
                full_pixel_polarization_csd_relative_frobenius=csd_relative,
                cross_spectrum_dimension=180, frequencies=len(fr))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--paired', default='probe-r1')
    parser.add_argument('--analytic', default='analytic-r1')
    args = parser.parse_args()
    root = args.campaign.resolve()
    require(not (root/'SCIENCE_REVIEW.json').exists(), 'Refusing to overwrite review')
    paired, analytic = root/args.paired, root/args.analytic
    for run in (paired, analytic):
        require((run/'EXIT_CODE').read_text().strip() == '0', 'Probe runner failed')
        require('[binding] result=PASS failures=0' in (run/'binding.log').read_text(), 'Binding missing')
    pic = json.loads((paired/'pic-dependency.json').read_text())
    oracle = json.loads((paired/'window/oracle-results.json').read_text())
    wave = json.loads((analytic/'analytic/analytic-results.json').read_text())
    require(pic['status'] == 'PASS' and oracle['pass'] and wave['pass'], 'Probe did not pass')
    require(len(pic['samples']) == 48, 'PIC case count changed')
    q, m, eps0, c = 1.602176634e-19, 9.1093837139e-31, 8.8541878128e-12, 299792458.
    dt, volume = 1e-6/c, 4e-6*4e-6*1e-6
    gamma = np.sqrt(1+.2**2); beta = .2/gamma
    response_errors, continuity, responses = [], [], []
    for s in pic['samples']:
        observed = s['target_ex_V_per_m']
        active = s['particle'] and s['source_offset_cells'] == 2 and s['source_ex_V_per_m'] > 0
        if active:
            f, sub = s['initial_z_fraction'], s['substeps']
            slope = -q*q*dt*dt*beta*(f-.5)/(m*gamma*eps0*volume)*.75*(1-f-beta/2)*(sub+1)/(2*sub)
            error = abs(observed/(slope*s['source_ex_V_per_m'])-1)
            require(error < 2e-5, 'Independent PIC analytic prediction failed')
            require(abs(observed)>1e-12 and observed == s['ampere_current_contribution_V_per_m'], 'PIC path not isolated')
            require(s['gathered_ex_V_per_m'] == 0 and s['gathered_by_T'] != 0, 'Wrong gather path')
            response_errors.append(float(error)); responses.append(abs(observed))
            continuity.append(s['continuity_relative_residual'])
        else:
            require(observed == 0 and s['target_jx_A_per_m2'] == 0, 'PIC negative control nonzero')
        require(s['ex_at_k_minus_1_V_per_m'] == 0, 'Farther output control nonzero')
    require(max(continuity)<1e-12, 'Continuity failed')
    oa = [plane_review(paired/'window', case, 'oracle', True) for case in oracle['cases']]
    aa = [plane_review(analytic/'analytic', case, 'window', False) for case in wave['cases']]
    require(len(oa) == len(aa) == 6, 'Expected six cases per window probe')
    require(len({case['window_bytes_per_grid'] for case in wave['cases']}) == 1, 'Analytic grid memory is not fixed')
    for case in wave['cases']:
        require(case['full_reference_grid_count'] == 0, 'Analytic candidate depends on full reference')
        require(case['cell_updates_per_grid'] == wave['nx']*wave['ny']*48*case['steps'], 'Analytic work count changed')
    for case in oracle['cases']:
        require(case['reference_grid_count'] == 2, 'Oracle reference cost omitted')
        require(case['reference_cell_updates_per_grid'] == oracle['nx']*oracle['ny']*(case['steps']+80)*case['steps'], 'Reference count changed')
        require(case['window_cell_updates_per_grid'] == oracle['nx']*oracle['ny']*48*case['steps'], 'Window count changed')
    report = dict(pass_=True, pic_samples=len(pic['samples']), pic_positive_cases=len(responses),
        max_independent_pic_slope_error=max(response_errors), max_pic_continuity_residual=max(continuity),
        pic_target_response_range_V_per_m=[min(responses),max(responses)],
        oracle_state_scalars_checked=sum(c['field_values_checked'] for c in oracle['cases']),
        oracle_plane_scalars=sum(c['scalar_plane_values'] for c in oa),
        oracle=oa, analytic=aa,
        analytic_grid_bytes_per_grid=wave['cases'][0]['window_bytes_per_grid'],
        limitations=['Oracle reference generation is part of total work and remains non-linear in this experiment.',
                     'Oracle source/initial arrays are manufactured, not a Gauss-compatible charged-beam model.',
                     'Independent analytic boundary supply is established only for known discrete vacuum waves.',
                     'Cross spectra are single-realization field products, not an ensemble mutual-intensity measurement.',
                     'Native E/B half-step timing is preserved; no radiation separation or physical energy-flux accuracy is claimed.',
                     'No full PIC window, held injection, moving CPML, MPI repartition, checkpoint support or production speedup.'])
    save(root/'SCIENCE_REVIEW.json', report)
    (root/'PASS').write_text('First-stage field-window and PIC-dependency probes independently verified; production feasibility remains open.\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('oracle','analytic','limitations')},indent=2))


if __name__ == '__main__':
    main()
