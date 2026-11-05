"""Bounded closed-only raw-plane and strict-ledger critical-mode comparison.

A projection on the staggered axial cutoff signature is a compatibility test,
not a unique modal energy measurement. It is interpreted with the independent
old/new field overlap and with the actual-substep BF control.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re

import h5py
import numpy as np

LD = np.longdouble
C = LD(299792458)


def require(value, message):
    if not value:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def identity(path):
    p = Path(path)
    s = p.stat()
    return dict(path=str(p.resolve()), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def fingerprint(ev, hv, fraction, centered=False):
    """Project a (transverse E,cB) sample onto the native Nyquist signature."""
    a, b = np.asarray(ev, dtype=LD), np.asarray(hv, dtype=LD)
    require(a.shape == b.shape and a.size > 0, 'Invalid field pair')
    require(np.isfinite(a).all() and np.isfinite(b).all(), 'Nonfinite selected field')
    if centered:
        a, b = a-a.mean(), b-b.mean()
    f = LD(fraction) % 1
    ee, hh = 1-2*f, -2*min(f, 1-f)
    norm = np.sum(a*a+b*b)
    template_norm = ee*ee+hh*hh
    require(template_norm > 0, 'Degenerate cutoff signature')
    projection = (ee*a+hh*b)/np.sqrt(template_norm)
    orthogonal = (ee*b-hh*a)/np.sqrt(template_norm)
    projected_energy = np.sum(projection*projection)
    orthogonal_energy = np.sum(orthogonal*orthogonal)
    residual = orthogonal_energy/norm if norm else LD(0)
    require(abs(projected_energy+orthogonal_energy-norm) <= LD('1e-14')*max(norm, LD('1e-100')), 'Projection identity differs')
    return dict(grid_fraction=float(f),
        E_boost_rms_V_per_m=float(np.sqrt(np.mean(a*a))),
        cB_boost_rms_V_per_m=float(np.sqrt(np.mean(b*b))),
        compatible_projection_rms_V_per_m=float(np.sqrt(projected_energy/a.size)),
        orthogonal_rms_V_per_m=float(np.sqrt(orthogonal_energy/a.size)),
        signature_residual_amplitude=float(np.sqrt(residual)),
        total_pair_rms_V_per_m=float(np.sqrt(norm/a.size)),
        predicted_H_over_E=float(hh/ee) if abs(ee) > LD('1e-12') else None,
        observed_regression_H_over_E=float(np.sum(a*b)/np.sum(a*a)) if np.any(a) else None)


def field_overlap(sums):
    """Exact native positive-energy overlap with ledger's actual endpoints."""
    old = sums['native_begin_J']
    # E endpoint in TE is after electric CPML and before inlet overwrite.
    # B endpoint in TB is before magnetic CPML. Do not subtract electric CPML.
    new = sums['native_end_J']-sums['remap_J']-sums['inlet_overwrite_J']-sums['magnetic_cpml_J']
    difference = sums['time_electric_J']-sums['time_magnetic_J']
    require(np.all(old > 0) and np.all(new > 0), 'Nonpositive overlap norm')
    cosine = (old+new-difference)/(2*np.sqrt(old*new))
    require(np.isfinite(cosine).all(), 'Nonfinite field overlap')
    require(np.all(abs(cosine) <= 1+LD('1e-8')), 'Field overlap exceeds inner-product bound')
    return old, new, difference, cosine


def read_ledger(directory, expected_steps, expected_ranks):
    paths = sorted((directory/'ledger').glob('discrete-field-rank-*.csv'))
    require([p.name for p in paths] == [f'discrete-field-rank-{k:05d}.csv' for k in range(expected_ranks)], 'Missing/extra strict-ledger ranks')
    keys = ['native_begin_J', 'native_end_J', 'remap_J', 'inlet_overwrite_J', 'magnetic_cpml_J', 'time_electric_J', 'time_magnetic_J']
    sums = {k: np.zeros(expected_steps, dtype=LD) for k in keys}
    for p in paths:
        require(p.with_suffix('.complete').read_text().split() == ['1', str(expected_steps)], 'Unclosed strict ledger: '+str(p))
        with p.open() as f:
            rows = list(csv.DictReader(f))
        require([int(x['step']) for x in rows] == list(range(1, expected_steps+1)), 'Missing/repeated ledger step')
        for k in keys:
            values = np.array([LD(row[k]) for row in rows])
            require(np.isfinite(values).all(), 'Nonfinite ledger')
            sums[k] += values
    return sums


def actual_geometry(root, directory, geometry):
    plan = json.loads((directory/'run-plan.json').read_text())
    import yaml
    cfg = yaml.safe_load((directory/'run.yaml').read_text())
    g, dz, dt = LD(geometry['boost_gamma']), LD(geometry['dz_m']), LD(geometry['dt_box_s'])
    beta = np.sqrt(1-1/(g*g))
    # Use the exact double beta written by the running program, below.
    energy = directory/'detectors/entry-transport-energy.h5'
    with h5py.File(energy, 'r') as f:
        group = f['energy_plane']
        require(int(group['complete'][()]) == 1, 'Unclosed energy plane')
        t0 = LD(group.attrs['time_lab_origin_s'])
        actual_gamma = LD(group.attrs['boost_gamma'])
        beta = LD(group.attrs['boost_beta'])
    require(abs(actual_gamma-g) < LD('1e-12')*g, 'Actual boost differs')
    require(abs(t0-LD(geometry['predicted_frame_time_origin_lab_s'])) < dt/g*LD('1e-5'), 'Frame origin differs from prepared geometry')
    anchor = LD(cfg['beam']['reference']['activation_plane_z'])
    mesh_center = LD(geometry['mesh_center_z_m'])
    reference_origin = mesh_center-LD(geometry['reference_nz'])*dz/2
    log = (directory/'aprl.log').read_text()
    hit = re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+)', log)
    require(hit is not None, 'Missing moving-window geometry')
    require(int(hit[1]) == geometry['reference_nz'] and int(hit[2]) == geometry['active_nz'], 'Actual field/reference mesh differs')
    offset = int(hit[3])
    with (directory/'reference/window-rank-0.csv').open() as f:
        first = next(csv.DictReader(f))
    require(int(first['step']) == 1, 'Missing first remap')
    origin0 = LD(first['origin_box_z_m'])+dz
    require(abs(origin0-(reference_origin+offset*dz)) < dz*LD('1e-9'), 'Actual moving origin differs')
    translation = LD(geometry['predicted_box_translation_m'])
    # The current program does not emit boxZOrigin at full precision. This
    # value is reconstructed from the exact frozen input by prepare.py.
    f0 = ((translation-origin0)/dz-anchor/(g*dz)) % 1
    predicted_f = LD(geometry['predicted_entry_fraction'])
    require(abs(((f0-predicted_f+LD('.5')) % 1)-LD('.5')) < LD('1e-8'), 'Prepared plane/grid fraction differs')
    return dict(g=g, beta=beta, dz=dz, dt=dt, t0=t0, anchor=anchor,
                translation=translation, origin0=origin0, plan=plan,
                coordinate_note='Frame time/gamma/beta and moving-grid origin checked against completed output; box translation reconstructed from frozen input by prepare.py, not directly recorded at full precision by the program.')


def require_case_complete(directory):
    p = Path(directory)/'PIC_COMPLETE.json'
    require(p.exists(), 'Missing final PIC completion marker: '+str(p))
    complete = json.loads(p.read_text())
    require(complete.get('analysis_complete') is True and complete.get('stage') == 'full', 'Full case has not completed final integrity checks: '+str(directory))


def analyze_mode(root, mode, geometry, targets):
    directory = root/mode
    require_case_complete(directory)
    result = json.loads((directory/'RESULT.json').read_text())
    require(result.get('program_integrity_pass') is True and result['stage'] == 'full', 'Require a completed accepted full control before opening raw fields: '+mode)
    steps, ranks = result['steps'], result['ranks']
    geo = actual_geometry(root, directory, geometry)
    plan = geo['plan']
    require(abs(steps-plan['expected_steps']) <= 2 and ranks == plan['ranks'], 'Actual layout/steps differs')
    sums = read_ledger(directory, steps, ranks)
    old, new, difference, cosine = field_overlap(sums)
    overlap = []
    for i in range(steps):
        overlap.append(dict(mode=mode, step=i+1, native_old_J=float(old[i]),
            matched_native_new_J=float(new[i]), actual_delta_square_J=float(difference[i]),
            energy_weighted_old_new_cosine=float(cosine[i])))
    rows, sources = [], []
    for plane in ['entry', 'und-entry', 'p1']:
        path = directory/'detectors'/(plane+'.h5')
        before = identity(path)
        with h5py.File(path, 'r') as f:
            group = f['field_plane']
            require(int(group['complete'][()]) == 1, 'Unclosed raw field plane')
            times = group['time_s'][:].astype(LD)
            require(len(times) == len(group['fields']) == int(group['committed_samples'][()]), 'Raw field sample commit count differs')
            z = LD(group.attrs['plane_z_m'])
            n_float = (times-geo['t0']-geo['beta']*(z-geo['anchor'])/C)*geo['g']/geo['dt']
            n = np.rint(n_float).astype(np.int64)
            require(len(n) and np.max(abs(n_float-n)) < LD('1e-3'), 'Raw time does not correspond to an integer Maxwell step')
            require(len(n) == len(np.unique(n)), 'Duplicate Maxwell-step field frames')
            chosen = sorted(set(int(q) for q in targets if q in n))
            for wanted in chosen:
                index = int(np.flatnonzero(n == wanted)[0])
                v = group['fields'][index]
                electric, magnetic = v['electric_V_per_m'].astype(LD), v['magnetic_T'].astype(LD)*C
                ex = geo['g']*(electric[...,0]-geo['beta']*magnetic[...,1])
                ey = geo['g']*(electric[...,1]+geo['beta']*magnetic[...,0])
                hx = geo['g']*(magnetic[...,0]+geo['beta']*electric[...,1])
                hy = geo['g']*(magnetic[...,1]-geo['beta']*electric[...,0])
                xi = (geo['translation']-geo['origin0'])/geo['dz']+(z-geo['anchor'])/(geo['g']*geo['dz'])+(1-geo['beta'])*wanted
                nx, ny = ex.shape
                require(nx >= 40 and ny >= 40, 'Field mesh too small for selected apertures')
                apertures = [('physical', (slice(10,-10), slice(10,-10))),
                             ('beam_core_0p8mm', (slice(nx//2-10,nx//2+10),slice(ny//2-10,ny//2+10)))]
                for region, sl in apertures:
                    for component, ev, hv in [('x', ex[sl], hy[sl]), ('y', ey[sl], -hx[sl])]:
                        for centered in [False, True]:
                            rows.append(dict(mode=mode, plane=plane, z_lab_m=float(z), maxwell_step=wanted,
                                region=region, spatial_mean_removed=centered, component=component,
                                **fingerprint(ev, hv, xi%1, centered)))
        require(before == identity(path), 'Raw file changed during closed analysis')
        sources.append(dict(**before, selected_frames=chosen))
    return dict(rows=rows, overlap=overlap, sources=sources,
        min_field_overlap=float(cosine.min()), field_overlap_final=float(cosine[-1]),
        coordinate_note=geo['coordinate_note'])


def write_csv(path, rows):
    require(bool(rows), 'No rows to write')
    with Path(path).open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare(root):
    root = Path(root).resolve()
    geometry = json.loads((root/'GEOMETRY.json').read_text())
    for mode in ['off', 'on']:
        require_case_complete(root/mode)
        p = root/mode/'RESULT.json'
        require(p.exists() and json.loads(p.read_text()).get('program_integrity_pass') is True, 'Both full control results must be closed before mode comparison')
    steps = geometry['steps']
    targets = sorted(set([0,10,20,50,100,200,400,800,1000,1500,2000,2500,3000,3500,steps]
                         + list(np.rint(np.linspace(0,steps,21)).astype(int))))
    targets = [int(n) for n in targets]
    cases = {mode: analyze_mode(root, mode, geometry, targets) for mode in ['off','on']}
    rows = [row for data in cases.values() for row in data['rows']]
    overlaps = [row for data in cases.values() for row in data['overlap']]
    keys = lambda row: tuple(row[x] for x in ['plane','maxwell_step','region','spatial_mean_removed','component'])
    maps = {mode: {keys(row):row for row in data['rows']} for mode,data in cases.items()}
    require(set(maps['off']) == set(maps['on']), 'Paired raw sample support differs')
    comparisons=[]
    for key in sorted(maps['off']):
        off,on = maps['off'][key],maps['on'][key]
        comparisons.append(dict(plane=key[0],maxwell_step=key[1],region=key[2],spatial_mean_removed=key[3],component=key[4],
            projection_rms_off_V_per_m=off['compatible_projection_rms_V_per_m'],
            projection_rms_on_V_per_m=on['compatible_projection_rms_V_per_m'],
            on_over_off_projection_rms=on['compatible_projection_rms_V_per_m']/off['compatible_projection_rms_V_per_m'] if off['compatible_projection_rms_V_per_m'] else None,
            residual_amplitude_off=off['signature_residual_amplitude'],
            residual_amplitude_on=on['signature_residual_amplitude']))
    out = root/'mode-review'
    out.mkdir(exist_ok=False)
    write_csv(out/'fingerprint.csv', rows)
    write_csv(out/'field_overlap.csv', overlaps)
    write_csv(out/'paired_fingerprint.csv', comparisons)
    theta = np.pi/10
    response = 1-np.sin(theta/2)**6
    plot(out, cases)
    save(out/'RESULT.json', dict(complete=True, physical_acceptance=False, frame='boost',
        new_simulation=False, selected_maxwell_steps=targets, cases=cases, comparisons=comparisons,
        signature='native axial critical Hy=-Ex,Hx=Ey; after staggered interpolation e=1−2f,h=−2min(f,1−f), so e*H−h*E=0 for H=(cBy,−cBx)',
        overlap_formula='cos=(Uold+Unew−(TE−TB))/(2 sqrt(Uold Unew)); Unew=end−remap−inlet−magneticCPML; electricCPML is already included in TE',
        target_shape_response=response, target_source_times_gather_transfer=response**2,
        transfer_note="This is a spatial coupling transfer, not a predicted final FEL power error.",
        cautions=['A two-component signature projection is compatibility, not a unique measured mode energy. Physical fields can project onto it.',
                  'Actual reduction of the signature and of field-step sign reversal must be checked; a filter switch alone does not prove suppression.',
                  'No causal conclusion for the original full bunch unless the unfiltered crop reproduces the early fast BF rise.',
                  'BF must be checked with actual converged-substep endpoint crossing interpolation as well as legacy full-step interpolation.',
                  'The cutoff filter leaves the Cowan axial vacuum dispersion relation unchanged, but changes spatial particle coupling.']))
    return out


def plot(out, cases):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-critical-mode-review-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.6),layout='constrained')
    for mode,color in [('off','#24648c'),('on','#b18428')]:
        a=cases[mode]['overlap']
        axes[0,0].plot([r['step'] for r in a],[r['energy_weighted_old_new_cosine'] for r in a],color=color,label=mode)
        for plane,style in [('entry',':'),('und-entry','--'),('p1','-')]:
            r=[q for q in cases[mode]['rows'] if q['plane']==plane and q['region']=='beam_core_0p8mm' and q['component']=='x' and not q['spatial_mean_removed']]
            n=[q['maxwell_step'] for q in r]
            axes[0,1].semilogy(n,[max(q['compatible_projection_rms_V_per_m'],1e-30) for q in r],style,color=color,label=mode+' '+plane)
            axes[1,0].plot(n,[q['signature_residual_amplitude'] for q in r],style,color=color,label=mode+' '+plane)
            axes[1,1].semilogy(n,[max(q['E_boost_rms_V_per_m'],1e-30) for q in r],style,color=color,label=mode+' '+plane)
    axes[0,0].set_title('Actual field-step old/new overlap')
    axes[0,0].set_ylabel('Native-energy-weighted cosine')
    axes[0,0].set_ylim(-1.05,1.05)
    axes[0,1].set_title('Compatible cutoff-signature projection')
    axes[0,1].set_ylabel('Boost pair RMS (V/m)')
    axes[1,0].set_title('Cutoff-signature mismatch amplitude')
    axes[1,0].set_ylabel('Orthogonal / total pair RMS')
    axes[1,0].set_ylim(0,1.05)
    axes[1,1].set_title('Raw boost Ex in beam-core aperture')
    axes[1,1].set_ylabel('RMS Ex (V/m)')
    for ax in axes.flat:
        ax.set_xlabel('Maxwell step')
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    fig.suptitle('Critical-mode control: measured fields and independent field-step overlap\nBeam-core 0.8 mm aperture; x pair; no near-field subtraction',fontsize=12)
    fig.savefig(out/'mode_comparison.png',dpi=200)
    fig.savefig(out/'mode_comparison.pdf')
    plt.close(fig)


def self_test(output):
    output=Path(output)
    rng=np.random.default_rng(20261105)
    pattern=rng.normal(size=(20,20))
    results=[]
    for f in [.02,.25,.5,.8,.95]:
        ee,hh=1-2*f,-2*min(f,1-f)
        a=fingerprint(ee*pattern,hh*pattern,f)
        require(a['signature_residual_amplitude'] < 1e-14, 'Exact cutoff signature not recognized')
        b=fingerprint((ee-.2*hh)*pattern,(hh+.2*ee)*pattern,f)
        require(abs(b['signature_residual_amplitude']-.2/np.sqrt(1.04)) < 1e-14, 'Orthogonal perturbation was not detected')
        results.append(dict(fraction=f,exact_residual=a['signature_residual_amplitude'],perturbed_residual=b['signature_residual_amplitude']))
    x=np.array([1,2,3],dtype=LD)
    for cosine in [-1,-.6,0,.4,1]:
        sums=dict(native_begin_J=x,native_end_J=x,remap_J=np.zeros(3,dtype=LD),inlet_overwrite_J=np.zeros(3,dtype=LD),magnetic_cpml_J=np.zeros(3,dtype=LD),time_electric_J=2*x*(1-cosine),time_magnetic_J=np.zeros(3,dtype=LD))
        *_,got=field_overlap(sums)
        require(np.max(abs(got-cosine))<1e-15,'Overlap sign/magnitude differs')
    # Deliberately nonzero CPML and remap avoid a test that only mirrors a
    # zero-exchange special case. TE already contains electric CPML.
    sums=dict(native_begin_J=x,native_end_J=2*x,remap_J=.2*x,inlet_overwrite_J=.3*x,magnetic_cpml_J=.5*x,time_electric_J=4*x,time_magnetic_J=np.zeros(3,dtype=LD))
    require(np.max(abs(field_overlap(sums)[-1]+1))<1e-15,'Matched-endpoint exchange signs differ')
    save(output,dict(passed=True,synthetic_fingerprint=results,overlap_cosines_checked=[-1,-.6,0,.4,1],nonzero_exchange_endpoint_test=True,scope='Analysis algebra only; not a PIC physical acceptance'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path,nargs='?')
    p.add_argument('--self-test',type=Path)
    a=p.parse_args()
    if a.self_test:
        self_test(a.self_test)
    else:
        require(a.root is not None,'Provide a completed campaign root')
        print(compare(a.root))
