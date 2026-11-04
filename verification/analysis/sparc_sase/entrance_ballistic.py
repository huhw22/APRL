"""Read-only entrance phase-space transport, without a PIC/GENESIS simulation.

The primary cohort is the complete input: no selection by later particle loss.
Fawley phase clones are retained exactly, including in permutation controls.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c, e


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist())+'\n')


def table(path, rows):
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    s = Path(path).stat()
    return dict(path=str(Path(path).resolve()), size=s.st_size,
                mtime_ns=s.st_mtime_ns, sha256=h.hexdigest())


def read_plane(path):
    with h5py.File(path, 'r') as f:
        g = f['particle_plane']
        assert int(g['complete'][()]) == 1
        a = g['records'][:]
        assert len(a) == int(g['committed_records'][()])
    a = a[np.argsort(a['source_id'])]
    assert np.all(np.diff(a['source_id']) > 0)
    assert np.all(a['crossing_direction'] == 1)
    for key in a.dtype.names:
        if a[key].dtype.kind == 'f':
            assert np.isfinite(a[key]).all()
    return a


def phases(t, anchor, wavelength):
    cycles = (np.asarray(t, dtype=np.longdouble)-np.longdouble(anchor))*np.longdouble(c/wavelength)
    angle = np.asarray(np.remainder(cycles, 1)*np.longdouble(2*np.pi), float)
    return np.cos(angle)+1j*np.sin(angle)


def csum(values, groups, count):
    return (np.bincount(groups, weights=values.real, minlength=count)
            +1j*np.bincount(groups, weights=values.imag, minlength=count))


def summarize(weighted, weights, width, nbase=64):
    # weighted already contains charge*complex phase, one value per beamlet.
    groups = np.arange(len(weights))//(nbase*width)
    ng = groups[-1]+1
    charge = np.bincount(groups, weights=weights, minlength=ng)
    valid = charge > 0
    local = np.full(ng, np.nan+1j*np.nan)
    local[valid] = csum(weighted, groups, ng)[valid]/charge[valid]
    global_b = weighted.sum()/weights.sum()
    rms = np.sqrt(np.sum(charge[valid]*abs(local[valid])**2)/charge.sum())
    mean_abs = np.sum(charge[valid]*abs(local[valid]))/charge.sum()
    assert abs(global_b) <= rms+1e-12
    return dict(global_b=float(abs(global_b)), global_real=float(global_b.real),
        global_imag=float(global_b.imag), local_rms=float(rms),
        phase_alignment=float(abs(global_b)/mean_abs) if mean_abs else 0.), local


def lag(gamma, q):
    uz = np.sqrt(gamma*gamma-q)
    return q/(uz*(gamma+uz))


def run(top, out):
    top, out = Path(top).resolve(), Path(out).resolve()
    out.mkdir(exist_ok=False)
    root = top/'compact-production-r1'
    plan = json.loads((root/'plan.json').read_text())
    src = top/'compact-preflight-r1/particles.h5'
    paths = [src, root/'plan.json', root/'boris/run.yaml',
        top/'compact-preflight-r1/NOISE_VALIDATION.json']
    for case in ['boris', 'energy']:
        complete = json.loads((root/case/'PIC_COMPLETE.json').read_text())
        assert complete['run']['exit_code'] == 0
        paths += [root/case/'detectors'/(name+'-particles.h5')
                  for name in ['entry', 'und-entry', 'p1']]
    before = [digest(p) for p in paths]
    with h5py.File(src, 'r') as f:
        a = f['particles/records'][:]
    n, nbins, nslice = len(a), plan['NBINS'], plan['active_slices']
    assert n == plan['input_particles'] == nslice*plan['NPART']
    assert np.array_equal(a['source_id'], np.arange(1, n+1))
    assert nbins == 16 and plan['NPART'] == 1024
    v = a['proper_velocity'].reshape(-1, nbins, 3)
    xy = a['plane_position_m'].reshape(-1, nbins, 2)
    w = a['macro_weight'].reshape(-1, nbins)
    assert np.all(v == v[:, :1]) and np.all(xy == xy[:, :1])
    assert np.all(w == w[:, :1]) and np.all(w > 0)
    # Compression is exact because every clone shares the same external orbit.
    u, pos = v[:, 0].copy(), xy[:, 0].copy()
    weights = w.sum(axis=1)
    gamma = np.sqrt(1+np.sum(u*u, axis=1))
    gamma_bar = np.sum(weights*gamma)/weights.sum()
    delta = (gamma-gamma_bar)/gamma_bar
    sigma = np.sqrt(np.sum(weights*delta**2)/weights.sum())
    anchor = float(np.median(a['arrival_time_offset_s']))
    wavelength = plan['physics']['wavelength_m']
    K, period = plan['physics']['K'], plan['physics']['period_m']
    ku = 2*np.pi/period
    prep, length = .01, 77*period
    z_end = prep+length
    q0 = 1+u[:, 0]**2+u[:, 1]**2
    drift_lag = lag(gamma, q0)
    y_entry = pos[:, 1]+prep*u[:, 1]/u[:, 2]
    q_focus = .5*K*K*ku*ku*y_entry*y_entry
    # End-envelope scalar quadrature only, no particle trajectory integration.
    grid = np.linspace(0, length, 77*512+1)
    dist = np.minimum(grid, length-grid)
    env = np.where(dist < period, .5*(1-np.cos(np.pi*dist/period)), 1.)
    g = env*np.cos(ku*grid)
    from scipy.integrate import cumulative_simpson
    g2_int = cumulative_simpson(g*g, x=grid, initial=0.)
    assert abs(g2_int[-1]-(length/2-5*period/8)) < 1e-12
    ref_gamma = plan['physics']['gamma']

    def shifts(z, model):
        d, l = min(z, prep), max(0., z-prep)
        reference = d*lag(ref_gamma, 1.)+l*lag(ref_gamma, 1.+K*K/2)
        if model == 'drift_exact':
            return z*drift_lag-reference
        if model == 'energy_only':
            return d*lag(gamma, 1.)+l*lag(gamma, 1.+K*K/2)-reference
        if model == 'undulator_average':
            return d*drift_lag+l*lag(gamma, q0+K*K/2)-reference
        if model in ['undulator_ramp', 'undulator_focus']:
            avg_g2 = float(np.interp(l, grid, g2_int))/l if l else 0.
            extra = K*K*avg_g2
            if model == 'undulator_focus':
                # Harmonic transverse action: uy^2+(gamma*kappa*y)^2.
                extra = extra+q_focus*(2*avg_g2)
            return d*drift_lag+l*lag(gamma, q0+extra)-reference
        raise ValueError(model)

    old = top/'bunching_growth_20261105/reviewed/analysis'
    old_result = json.loads((old/'RESULT.json').read_text())
    waves = dict(old_result['wavelengths_nm'])
    with np.load(old/'cohort.npz') as f:
        common_ids = f['source_id']
    common = np.zeros(n, bool); common[common_ids.astype(int)-1] = True
    masks = {'full': np.ones((len(weights), nbins), bool),
             'legacy_common': common.reshape(-1, nbins)}
    cohort_weights = {name: (w*mask).sum(axis=1) for name, mask in masks.items()}
    # The legacy subset can lose complete tail groups; summarize omits zero charge.
    rows, local_arrays, validation = [], {}, {}
    design_phase = phases(a['arrival_time_offset_s'], anchor, wavelength).reshape(-1, nbins)
    initial_weighted = (w*design_phase).sum(axis=1)
    for width in [5, 10, 20]:
        _, loc = summarize(initial_weighted, weights, width)
        local_arrays[f'entry_full_{width}'] = loc
    direct = np.sum(a['macro_weight']*design_phase.ravel())/weights.sum()
    validation['beamlet_compression_complex_error'] = float(abs(direct-initial_weighted.sum()/weights.sum()))
    assert validation['beamlet_compression_complex_error'] < 1e-14
    validation['physical_Ne'] = float(weights.sum())
    validation['physical_global_noise_reference'] = float(1/np.sqrt(weights.sum()))
    validation['wrong_independent_macro_reference'] = float(np.sqrt(np.sum(a['macro_weight']**2))/weights.sum())
    validation['physical_local10_noise_reference'] = float(np.sqrt(len(local_arrays['entry_full_10'])/weights.sum()))

    models = ['drift_exact', 'energy_only', 'undulator_average', 'undulator_ramp', 'undulator_focus']
    zgrid = np.unique(np.r_[np.linspace(0, z_end, 121), prep, z_end])
    for wave, wave_nm in waves.items():
        wave_m = wave_nm*1e-9
        ph = phases(a['arrival_time_offset_s'], anchor, wave_m).reshape(-1, nbins)
        for cohort, mask in masks.items():
            weighted0 = (w*mask*ph).sum(axis=1)
            cw = cohort_weights[cohort]
            for model in models:
                use_z = zgrid if wave == 'design' and cohort == 'full' else [0., prep, z_end]
                for z in use_z:
                    weighted = weighted0*np.exp(2j*np.pi/wave_m*shifts(z, model))
                    summary, loc = summarize(weighted, cw, 10)
                    rows.append(dict(kind='analytic', cohort=cohort, model=model,
                        wave=wave, wavelength_nm=wave_nm, z_m=float(z), **summary))
                    if wave == 'design' and z == z_end:
                        local_arrays[f'{model}_{cohort}_p1'] = loc
                        for width in [5, 20]:
                            ss, ll = summarize(weighted, cw, width)
                            rows.append(dict(kind='analytic_width'+str(width), cohort=cohort,
                                model=model, wave=wave, wavelength_nm=wave_nm, z_m=float(z), **ss))
            ss, _ = summarize(weighted0, cw, 10)
            if cohort == 'legacy_common':
                prev = next(r for r in old_result['rows'] if r['case']=='boris' and
                    r['plane']=='entry-particles' and r['wavelength_label']==wave)
                validation['legacy_entry_'+wave] = float(abs(ss['local_rms']-prev['local_rms_10lambda']))
                assert validation['legacy_entry_'+wave] < 2e-10

    # Read only completed early planes. Their whole ID population is unchanged.
    input_t = a['arrival_time_offset_s']
    time_offset = None
    pic_cache = {}
    for case in ['boris', 'energy']:
        for plane, z in [('entry', 0.), ('und-entry', prep), ('p1', z_end)]:
            p = read_plane(root/case/'detectors'/(plane+'-particles.h5'))
            assert np.array_equal(p['source_id'], a['source_id'])
            assert np.max(abs(p['weight']-a['macro_weight'])) == 0
            if plane == 'entry':
                offset = np.median(p['time_s'].astype(np.longdouble)-input_t.astype(np.longdouble))
                validation[case+'_entry_time_nonconstant_max_s'] = float(np.max(abs(
                    p['time_s'].astype(np.longdouble)-input_t.astype(np.longdouble)-offset)))
                validation[case+'_entry_u_max'] = float(np.max(abs(p['proper_velocity']-a['proper_velocity'])))
                if time_offset is None:
                    time_offset = offset
                assert abs(offset-time_offset) < 1e-22
            # One physical reference phase for all models and PIC at each face.
            d, l = min(z, prep), max(0., z-prep)
            reference = d*lag(ref_gamma, 1.)+l*lag(ref_gamma, 1.+K*K/2)
            t_anchor = np.longdouble(anchor)+time_offset+np.longdouble((z+reference)/c)
            if plane == 'p1':
                pic_cache[case] = (p['time_s'].copy(), t_anchor)
            for wave, wave_nm in waves.items():
                ph = phases(p['time_s'], t_anchor, wave_nm*1e-9).reshape(-1, nbins)
                for cohort, mask in masks.items():
                    weighted = (w*mask*ph).sum(axis=1)
                    summary, loc = summarize(weighted, cohort_weights[cohort], 10)
                    rows.append(dict(kind='PIC', cohort=cohort, model=case, wave=wave,
                        wavelength_nm=wave_nm, z_m=z, **summary))
                    if wave == 'design' and plane == 'p1':
                        local_arrays[f'PIC_{case}_{cohort}_p1'] = loc
                        for width in [5, 20]:
                            ss, _ = summarize(weighted, cohort_weights[cohort], width)
                            rows.append(dict(kind='PIC_width'+str(width), cohort=cohort,
                                model=case, wave=wave, wavelength_nm=wave_nm, z_m=z, **ss))
                    if cohort == 'legacy_common':
                        prev = next(r for r in old_result['rows'] if r['case']==case and
                            r['plane']==plane+'-particles' and r['wavelength_label']==wave)
                        err = max(abs(summary['global_b']-prev['global_b']),
                                  abs(summary['local_rms']-prev['local_rms_10lambda']))
                        assert err < 2e-10
                        validation[f'legacy_{case}_{plane}_{wave}'] = float(err)
            del p

    # Connected moments: remove density times the uniform variance.
    moment_rows = []
    grp = np.arange(len(weights))//640
    charge = np.bincount(grp, weights=weights)
    mean = np.bincount(grp, weights=weights*delta)/charge
    centered = delta-mean[grp]
    var = np.bincount(grp, weights=weights*centered**2)/charge
    for wave, wave_nm in waves.items():
        ph = phases(input_t, anchor, wave_nm*1e-9).reshape(-1, nbins)
        wp = (w*ph).sum(axis=1)
        for name, values in [('density', np.ones_like(delta)), ('energy_connected', centered),
                             ('variance_connected', centered**2-var[grp])]:
            # Global moments have their own global centering.
            gv = np.ones_like(delta) if name=='density' else (
                delta if name=='energy_connected' else delta**2-sigma**2)
            global_m = np.sum(wp*gv)/weights.sum()
            lm = csum(wp*values, grp, len(charge))/charge
            norm = 1. if name=='density' else sigma if name=='energy_connected' else sigma**2
            moment_rows.append(dict(wave=wave, wavelength_nm=wave_nm, moment=name,
                global_real=float(global_m.real), global_imag=float(global_m.imag),
                global_abs=float(abs(global_m)),
                local_rms=float(np.sqrt(np.sum(charge*abs(lm)**2)/charge.sum())),
                local_normalized_rms=float(np.sqrt(np.sum(charge*abs(lm/norm)**2)/charge.sum()))))
            if wave=='design':
                local_arrays[name+'_moment'] = lm
    # Exact energy-only phase shear versus the hierarchy expansion.
    k = 2*np.pi/wavelength
    R = -(prep+length*(1+K*K/2))/gamma_bar**2
    alpha = k*R
    hierarchy_rows = []
    partial = np.ones(len(delta), complex)
    terms = np.ones(len(delta), complex)
    energy_linear_exact = initial_weighted*np.exp(1j*alpha*delta)
    exact_local = summarize(energy_linear_exact, weights, 10)[1]
    for order in range(7):
        if order:
            terms *= 1j*alpha*delta/order
            partial += terms
        ss, local = summarize(initial_weighted*partial, weights, 10)
        termlocal = summarize(initial_weighted*terms, weights, 10)[1]
        hierarchy_rows.append(dict(order=order, **ss,
            term_local_rms=float(np.sqrt(np.sum(charge*abs(termlocal)**2)/charge.sum())),
            complex_error_local_rms=float(np.sqrt(np.sum(charge*abs(local-exact_local)**2)/charge.sum()))))
    assert hierarchy_rows[-1]['complex_error_local_rms'] < 1e-7

    # Permute complete Fawley families within each optical slice; never split
    # the 16 equal-momentum clones. These are sensitivity controls, not new SASE seeds.
    rng = np.random.default_rng(20261105)
    controls = []
    shifts_end = {model: shifts(z_end, model).reshape(nslice, 64)
                  for model in ['energy_only', 'undulator_focus']}
    for repeat in range(128):
        permutation = np.argsort(rng.random((nslice, 64)), axis=1)
        for model, shift in shifts_end.items():
            randomized = np.take_along_axis(shift, permutation, axis=1).ravel()
            ss, _ = summarize(initial_weighted*np.exp(1j*k*randomized), weights, 10)
            controls.append(dict(repeat=repeat, model=model, **ss))
    # Keep a near-resonance spectrum. Local finite-window edges are not a new noise floor.
    spectra = []
    scan_nm = np.unique(np.r_[np.linspace(486., 497., 89), wavelength*1e9,
                              list(waves.values())])
    for wave_nm in scan_nm:
        wave_m = wave_nm*1e-9
        ph = phases(input_t, anchor, wave_m).reshape(-1, nbins)
        wp = (w*ph).sum(axis=1)
        for label, pred in [('entry', wp), ('ballistic', wp*np.exp(2j*np.pi/wave_m*shifts(z_end, 'undulator_focus')))]:
            ss, _ = summarize(pred, weights, 10)
            spectra.append(dict(model=label, wavelength_nm=float(wave_nm), **ss))
        for case, (t, ta) in pic_cache.items():
            ph = phases(t, ta, wave_m).reshape(-1, nbins)
            ss, _ = summarize((w*ph).sum(axis=1), weights, 10)
            spectra.append(dict(model=case, wavelength_nm=float(wave_nm), **ss))
    table(out/'bunching.csv', rows)
    table(out/'connected_moments.csv', moment_rows)
    table(out/'moment_hierarchy.csv', hierarchy_rows)
    table(out/'family_permutations.csv', controls)
    table(out/'bunching_spectrum.csv', spectra)
    np.savez_compressed(out/'local_complex.npz', **local_arrays,
        group_charge_e=charge, group_mean_delta=mean, group_variance_delta=var)
    summary = dict(input_particles=n, input_electrons=float(weights.sum()),
        gamma_mean=float(gamma_bar), sigma_delta=float(sigma), wavelength_m=wavelength,
        first_segment_end_m=z_end, energy_mixing_parameter_77=float(4*np.pi*77*sigma),
        energy_mixing_length_m=float(period/(4*np.pi*sigma)),
        noise_power_memory_loss_77=float(1-np.exp(-(4*np.pi*77*sigma)**2)),
        integral_g2_m=float(g2_int[-1]), model_names=models,
        source='Full frozen input; exact equal-orbit Fawley beamlet compression',
        primary_cohort='All input IDs; no conditioning on later losses',
        no_new_particle_or_field_simulation=True,
        model_limit='Constant particle energy; exact unmagnetized flight and averaged undulator flight. Natural focusing uses quadratic transverse action in the constant-amplitude bulk, with short-end approximation and scalar ramp integrals. No self-fields, radiation feedback, or exact 3D trajectories; model differences are not rigorous error bounds.',
        control_limit='128 deterministic within-slice whole-beamlet permutations; sensitivity distribution, not physical confidence intervals or independent SASE shots.',
        spectrum_limit='Raw particle bunching, no Hann window; off-design local BF includes material-window edge leakage.',
        rows_at_design_p1=[r for r in rows if r['wave']=='design' and r['z_m']==z_end and r['kind'] in ['analytic','PIC']],
        connected_moments=moment_rows, moment_hierarchy=hierarchy_rows)
    save(out/'RESULT.json', summary)
    after = [digest(p) for p in paths]
    assert before == after
    validation['sources_unchanged'] = True
    validation['source_manifest'] = before
    save(out/'VALIDATION.json', validation)
    render(out, rows, moment_rows, spectra, controls, local_arrays, summary)
    document(out, summary, rows, moment_rows, controls)
    save(out/'COMPLETE.json', dict(complete=True, new_simulations=False))
    print(json.dumps(summary['rows_at_design_p1'], indent=2), flush=True)


def render(out, rows, moments, spectra, controls, arrays, summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False,
        'axes.spines.right':False, 'lines.linewidth':1.8})
    fig, ax = plt.subplots(2, 2, figsize=(12.6, 8.8), layout='constrained')
    styles = [('drift_exact','Free drift','#929292',':'),
              ('energy_only','Energy dispersion only','#b18428','--'),
              ('undulator_focus','Undulator, angles + focusing','#222222','-')]
    for model, label, color, line in styles:
        r = [x for x in rows if x['cohort']=='full' and x['wave']=='design' and
             x['kind']=='analytic' and x['model']==model]
        for a, field in [(ax[0,0], 'local_rms'), (ax[0,1], 'global_b')]:
            a.semilogy([x['z_m'] for x in r], [x[field] for x in r], line, color=color, label=label)
    for case, label, color, marker in [('boris','PIC explicit','#24648c','o'),
                                      ('energy','PIC energy coupled','#8e4b75','^')]:
        r = [x for x in rows if x['kind']=='PIC' and x['cohort']=='full' and
             x['wave']=='design' and x['model']==case]
        for a, field in [(ax[0,0], 'local_rms'), (ax[0,1], 'global_b')]:
            a.semilogy([x['z_m'] for x in r], [x[field] for x in r], marker,
                       ls='none', ms=7, color=color, label=label)
    ax[0,0].set_title('(a) Local RMS BF: same input material groups')
    ax[0,1].set_title('(b) Global BF: complete bunch')
    for a in ax[0]:
        a.set_xlabel('Laboratory plane z (m)'); a.set_ylabel('Bunching factor')
        a.set_xlim(-.05, 2.26); a.grid(alpha=.2); a.legend(fontsize=8, loc='upper left')
    ax[0,0].set_ylim(3e-4, 9e-3); ax[0,1].set_ylim(1e-5, 9e-3)
    idx = np.arange(len(arrays['entry_full_10']))
    for key, label, color, line in [('entry_full_10','Entrance','#929292',':'),
        ('undulator_focus_full_p1','Ballistic undulator','#222222','-'),
        ('PIC_boris_full_p1','PIC explicit','#24648c','-'),
        ('PIC_energy_full_p1','PIC energy coupled','#8e4b75','--')]:
        ax[1,0].semilogy(idx, abs(arrays[key]), line, color=color, label=label)
    ax[1,0].set_title('(c) Local BF at first segment exit')
    ax[1,0].set_xlabel('Fixed input group (10 optical slices)')
    ax[1,0].set_ylabel('Local bunching factor'); ax[1,0].set_ylim(1e-5, 2e-2)
    ax[1,0].grid(alpha=.2); ax[1,0].legend(fontsize=8, ncol=2)
    for model, color, label in [('energy_only','#b18428','Energy only'),
                              ('undulator_focus','#222222','Angles + focusing')]:
        values = [x['local_rms']*1e4 for x in controls if x['model']==model]
        ax[1,1].hist(values, bins=18, histtype='step', lw=1.8, color=color, label=label+' permutations')
        observed = next(x['local_rms']*1e4 for x in rows if x['kind']=='analytic' and
            x['cohort']=='full' and x['wave']=='design' and x['model']==model and
            x['z_m']==summary['first_segment_end_m'])
        ax[1,1].axvline(observed, color=color, ls='--', label=label+' actual input')
    ax[1,1].set_title('(d) Preserve 16-clone families; permute within slices')
    ax[1,1].set_xlabel('Predicted local RMS BF ($10^{-4}$)'); ax[1,1].set_ylabel('Count (128 controls)')
    ax[1,1].legend(fontsize=8); ax[1,1].grid(alpha=.2)
    fig.suptitle('Entrance phase-space transport: no self-field or radiation feedback\n'
                 'Design wavelength 491.317 nm; full input; PIC shown only at stored planes', fontsize=13)
    fig.savefig(out/'entrance_ballistic_comparison.png', dpi=220, bbox_inches='tight'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.9), layout='constrained')
    for model, label, color, line in [('entry','Entrance','#929292',':'),
        ('ballistic','Ballistic undulator at p1','#222222','-'),
        ('boris','PIC explicit at p1','#24648c','-'),
        ('energy','PIC energy coupled at p1','#8e4b75','--')]:
        rr = [r for r in spectra if r['model']==model]
        for a, field in zip(axes, ['global_b','local_rms']):
            a.semilogy([r['wavelength_nm'] for r in rr], [r[field] for r in rr],
                       line, color=color, label=label)
    for a, title in zip(axes, ['Global BF','Local RMS BF (10 input slices)']):
        a.set_title(title); a.axvline(summary['wavelength_m']*1e9, color='#b18428', lw=1, ls=':')
        a.set_xlabel('Wavelength (nm)'); a.set_ylabel('Bunching factor'); a.grid(alpha=.2)
        a.legend(fontsize=8)
    fig.suptitle('Raw particle bunching spectrum; off-design local values include finite-window edges', fontsize=12)
    fig.savefig(out/'entrance_ballistic_spectrum.png', dpi=220, bbox_inches='tight'); plt.close(fig)


def document(out, result, rows, moments, controls):
    end = result['first_segment_end_m']
    rr = [r for r in rows if r['wave']=='design' and r['cohort']=='full'
          and r['z_m']==end and r['kind'] in ['analytic','PIC']]
    labels = {'drift_exact':'纯自由漂移', 'energy_only':'波荡器色散，仅实际能量',
        'undulator_average':'平均波荡器，含实际横向动量', 'undulator_ramp':'加入首尾包络积分',
        'undulator_focus':'再加入近轴自然聚焦作用量', 'boris':'已有PIC显式版', 'energy':'已有PIC能量版'}
    text = '# 实际入口相空间：无自场传播能否解释首段BF增长\n\n'
    baseline = next(r['local_rms'] for r in moments if r['wave']=='design' and r['moment']=='density')
    ballistic = next(r for r in rr if r['model']=='undulator_focus')
    observed = [r for r in rr if r['kind']=='PIC']
    text += (f"**入口原有能量、能散及横向相空间结构的自由输运不足以解释PIC首段增长。** "
             f"完整输入的局部RMS BF为{baseline:.8g}；含实际能量、角度、端部和近轴聚焦的平均轨道预测为"
             f"{ballistic['local_rms']:.8g}，仅上升{100*(ballistic['local_rms']/baseline-1):.3f}%。"
             f"两份已有PIC结果则为{observed[0]['local_rms']:.8g}与{observed[1]['local_rms']:.8g}，"
             f"约为入口的{observed[0]['local_rms']/baseline:.2f}和{observed[1]['local_rms']/baseline:.2f}倍。\n\n")
    text += '这个结论排除了这里计算的“既有入口相关项直接经无反馈输运抬高BF”解释，未判定PIC增长的其余来源；传播中产生的自场、辐射反馈或数值作用不在本参考内。\n\n'
    text += '只对已有输入和已闭合粒子面后处理，没有运行新的PIC、GENESIS或粒子轨迹仿真。\n\n'
    text += f"实际相对能散为{result['sigma_delta']:.8g}，77周期混合参数{result['energy_mixing_parameter_77']:.6f}，充分相混尺度约{result['energy_mixing_length_m']:.3f} m。\n\n"
    text += '主口径保留全部1,419,264宏粒子（约5.035亿电子）。入口到首段末两版均完整，不用后期退场筛选早期样本。原跨全部面的共同子集只作复现副口径，见CSV。\n\n'
    text += '![入口自由传播对照](entrance_ballistic_comparison.png)\n\n'
    text += '|2.166 m处，设计波长|全局BF|10切片局部RMS|\n|---|---:|---:|\n'
    for r in rr:
        text += f"|{labels[r['model']]}|{r['global_b']:.9g}|{r['local_rms']:.9g}|\n"
    text += '\n## 实际能量与能散相关项\n\n'
    text += r'采用 $C_1=\langle(\delta-\mu)e^{i\phi}\rangle$、$C_2=\langle[(\delta-\mu)^2-\sigma^2]e^{i\phi}\rangle$；全局和局部各自中心化，扣除密度BF乘均匀方差项。'+ '\n\n'
    text += '|设计波长相关项|全局模|局部RMS|除以相应能散尺度后的局部RMS|\n|---|---:|---:|---:|\n'
    for r in moments:
        if r['wave']=='design':
            text += f"|{r['moment']}|{r['global_abs']:.9g}|{r['local_rms']:.9g}|{r['local_normalized_rms']:.9g}|\n"
    text += '\n这些有限样本矩不是独立的额外种子。最终BF由原始相位和完整逐粒子飞行时间复数求和；高阶项及交叉项均包含。moment_hierarchy.csv另记录能量线性色散的0–6阶展开与精确指数之差。\n\n'
    text += '## 无反馈参考的物理范围\n\n'
    text += r'从原始z=0过面事件出发，前1 cm严格无磁飞行；波荡器段保持每个粒子的能量，计算平均纵向飞行时间。$c\Delta t=\int(\gamma/u_z-1)dz$；纯漂移用实际$u_z$，平均磁段用$u_z^2=\gamma^2-1-u_x^2-u_y^2-K^2/2$。'+'\n\n'
    text += r'自然聚焦修正在常幅平均磁段内保留$u_y^2+K^2k_u^2y^2/2$的近轴作用量不变量；首尾采用短端部近似与标量积分$\int g^2dz=L/2-5\lambda_u/8$，不声称端部内作用量严格保持。未求原三维磁场精确轨迹，未包含自场、空间电荷或辐射反馈。段内细线是平均轨道预测，不是新增PIC采样；几个参考之间的差值不是严格误差上界。'+'\n\n'
    text += '## Fawley噪声与置换\n\n'
    text += '每个光学切片含64个beamlet，每组16条相位克隆共享完全相同的动量和横向位置。计算严格保留这个结构。128次置换只在同一光学切片内交换完整beamlet的飞行时延，不拆散克隆；它们是相关性敏感性对照，不是独立真实电子噪声样本或SASE置信区间。不得用1/√Nmacro解释本输入的物理噪声。\n\n'
    for model in ['energy_only','undulator_focus']:
        val = np.array([x['local_rms'] for x in controls if x['model']==model])
        text += f"- {labels[model]}置换局部RMS：中位{np.median(val):.8g}，128次最小/最大{val.min():.8g}/{val.max():.8g}。\n"
    text += '\n![附近频段](entrance_ballistic_spectrum.png)\n\n'
    text += '这是粒子BF谱，不是场能谱，不使用Hann窗。非设计波长处的局部材料切片有有限窗口边缘项，必须与同波长入口对照。\n\n'
    text += '数据和核验：bunching.csv、connected_moments.csv、moment_hierarchy.csv、family_permutations.csv、bunching_spectrum.csv、local_complex.npz、RESULT.json、VALIDATION.json。全部输入及读取的原始证据已做前后SHA256一致核对。\n'
    (out/'README.md').write_text(text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('top'); parser.add_argument('out')
    args = parser.parse_args()
    run(args.top, args.out)
