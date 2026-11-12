"""Separate same-input bunching, full-pulse power, and reference-code controls."""
import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c

from prepare import save


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def grouped(current,moment,width=10):
    starts=np.arange(0,len(current),width)
    q=np.add.reduceat(current,starts)
    b=np.add.reduceat(current*moment,starts)
    use=q>0
    return float(np.sqrt(np.sum(abs(b[use])**2/q[use])/q[use].sum()))


def pic_bunching(path,wavelength):
    with h5py.File(path,'r') as f:
        g=f['particle_plane'];assert int(g['complete'][()])==1
        a=g['records'][:];z=float(g.attrs['plane_z_m'])
    phase=np.exp(2j*np.pi*c*(a['time_s'].astype(np.longdouble)-np.longdouble(z)/c)/wavelength).astype(np.complex128)
    bins=(a['source_id'].astype(np.int64)-1)//(10*1024)
    q=np.bincount(bins,weights=a['weight'])
    b=np.bincount(bins,weights=a['weight']*phase.real)+1j*np.bincount(bins,weights=a['weight']*phase.imag)
    use=q>0
    return dict(count=len(a),represented_electrons=float(a['weight'].sum()),
        global_b1=float(abs(b.sum())/q.sum()),
        source_group_local_rms_b1=float(np.sqrt(np.sum(abs(b[use])**2/q[use])/q.sum())))


def gen_series(d,p):
    with h5py.File(d/'reference.out.h5','r') as f:
        z=f['Lattice/zplot'][:];current=f['Beam/current'][:]
        moment=f['Beam/bunching'][:]*np.exp(1j*f['Beam/bunchingphase'][:])
        power=f['Field/power'][:]
        ds=float(f['Global/lambdaref'][0])*float(f['Global/sample'][0])
        global_energy=f['Field/Global/energy'][:]
        xs=f['Beam/xsize'][:];ys=f['Beam/ysize'][:]
        energy=f['Beam/energy'][:];spread=f['Beam/energyspread'][:]
        assert np.all(current==current[:1]),'Current unexpectedly changed: sorting or remapping?'
        assert np.allclose(power.sum(axis=1)*ds/c,global_energy,rtol=2e-12,atol=1e-24)
        lo=512;hi=lo+p['active_slices'];slice_indices=np.arange(lo,hi)[::-1]
        I=current[:,slice_indices];b=moment[:,slice_indices]
        totals=I.sum(axis=1)
        global_b=abs((I*b).sum(axis=1))/totals
        rms=np.array([grouped(x,y) for x,y in zip(I,b)])
        peak=power.max(axis=1)
        bx=np.sqrt((I*xs[:,slice_indices]**2).sum(axis=1)/totals)
        by=np.sqrt((I*ys[:,slice_indices]**2).sum(axis=1)/totals)
        eg=(I*energy[:,slice_indices]).sum(axis=1)/totals
        sg=np.sqrt((I*spread[:,slice_indices]**2).sum(axis=1)/totals)
    return dict(z_m=z,pulse_energy_J=global_energy,peak_power_W=peak,
                global_b1=global_b,source_group_local_rms_b1=rms,
                current_weighted_sigma_x_m=bx,current_weighted_sigma_y_m=by,
                current_weighted_gamma=eg,current_weighted_slice_sigma_gamma=sg)


def dump_spectrum(path,wavelength,halfwidth=.00056):
    with h5py.File(path,'r') as f:
        n=int(f['slicecount'][0]);nx=int(f['gridpoints'][0]);dx=float(f['gridsize'][0])
        dt=float(f['slicespacing'][0])/c
        x=(np.arange(nx)-(nx-1)/2)*dx
        # Fractional quadrature-cell intersections retain the exact PIC aperture.
        w=np.maximum(0,np.minimum(x+dx/2,halfwidth)-np.maximum(x-dx/2,-halfwidth))/dx
        ix=np.flatnonzero(w>0);area_weight=w[ix,None]*w[None,ix]
        values=np.empty((n,len(ix),len(ix)),np.complex128)
        for k in range(n):
            g=f[f'slice{k+1:06d}']
            field=g['field-real'][:]+1j*g['field-imag'][:]
            values[k]=field.reshape(nx,nx)[np.ix_(ix,ix)]
        # GENESIS dumped field units are sqrt(W per transverse quadrature cell).
        direct=abs(values)**2
        edges=float((direct[:8].sum()+direct[-8:].sum())/max(direct.sum(),1e-300))
        transform=np.fft.fft(values,axis=0)*dt
        frequency=np.fft.fftfreq(n,dt)
        # The envelope sign changes which sideband is positive, not this symmetric band.
        band=abs(frequency)<=.05*c/wavelength
        spectral=float(np.sum(abs(transform[band])**2*area_weight)/(n*dt))
        total=float(np.sum(direct*area_weight)*dt)
        parseval=float(np.sum(abs(transform)**2*area_weight)/(n*dt))
        assert abs(parseval-total)<2e-12*max(total,1e-300)
    return dict(rect_5pct_aperture_J=spectral,rect_all_envelope_aperture_J=total,
                full_envelope_edge8_fraction=edges,grid_dx_m=dx,gridpoints=nx,
                aperture_halfwidth_m=halfwidth,native_slices=n)


def analyze(root,p):
    root=Path(root);pic=Path(p['pic_campaign']);out=root/'analysis';out.mkdir(exist_ok=False)
    lam=p['original_plan']['physics']['wavelength_m']
    with (pic/'analysis-energy-binomial3/growth.csv').open() as f:
        old=list(csv.DictReader(f))[:6]
    pic_rows=[]
    entry=pic_bunching(pic/'energy-binomial3/detectors/entry-particles.h5',lam)
    pic_rows.append(dict(module=0,z_m=0.,rect_5pct_flux_J=0.,**entry))
    for j in range(1,7):
        bf=pic_bunching(pic/f'energy-binomial3/detectors/p{j}-particles.h5',lam)
        with np.load(pic/f'analysis-energy-binomial3/bandwidth-scan-r1/p{j:02d}_spectra.npz') as f:
            mask=abs(f['frequency_Hz']/(c/lam)-1)<=.05
            W=float(f['rect_flux_bin_J'][mask].sum())
        pic_rows.append(dict(module=j,z_m=float(old[j-1]['z_m']),rect_5pct_flux_J=W,**bf))
    table(out/'pic_exits.csv',pic_rows)
    series={};exits=[];spectra=[]
    for case in p['cases']:
        name=case['name'];d=Path(case['directory'])
        assert json.loads((d/'COMPLETE.json').read_text())['complete']
        s=gen_series(d,p);series[name]=s
        np.savez_compressed(out/(name+'-growth.npz'),**s)
        table(out/(name+'-growth.csv'),[dict(zip(s,[float(s[k][j]) for k in s])) for j in range(len(s['z_m']))])
        for j,z in [('entry',0.),*[(i,float(old[i-1]['z_m'])) for i in range(1,7)]]:
            k=int(np.argmin(abs(s['z_m']-z)))
            exits.append(dict(case=name,module=j,requested_z_m=z,
                              **{key:float(value[k]) for key,value in s.items()}))
        files=sorted(d.glob('reference.*.fld.h5'),key=lambda f:int(f.name.split('.')[1]))
        assert len(files)==2,(name,files)
        for module,f in zip([1,6],files):
            spectra.append(dict(case=name,module=module,path=str(f),**dump_spectrum(f,lam)))
    table(out/'genesis_exits.csv',exits);table(out/'genesis_spectra.csv',spectra)
    input_bf=pic_rows[0]['source_group_local_rms_b1']
    for name,s in series.items():
        assert abs(s['source_group_local_rms_b1'][0]/input_bf-1)<2e-7,(name,input_bf,s['source_group_local_rms_b1'][0])
    ratios=[]
    for r in spectra:
        a=pic_rows[r['module']]
        ratios.append(dict(case=r['case'],module=r['module'],
            PIC_rect_5pct_aperture_J=a['rect_5pct_flux_J'],
            GENESIS_rect_5pct_aperture_J=r['rect_5pct_aperture_J'],
            PIC_over_GENESIS=a['rect_5pct_flux_J']/r['rect_5pct_aperture_J']))
    table(out/'same_band_aperture_comparison.csv',ratios)
    controls=[]
    for a,b in [('longrange-installed','longrange-fixed'),('sc-off-halfstep','sc-off-fixed'),('longrange-fixed','sc-off-fixed')]:
        ra=next(x for x in exits if x['case']==a and x['module']==6)
        rb=next(x for x in exits if x['case']==b and x['module']==6)
        controls.append(dict(numerator=a,denominator=b,
            pulse_energy_ratio=ra['pulse_energy_J']/rb['pulse_energy_J'],
            local_bf_ratio=ra['source_group_local_rms_b1']/rb['source_group_local_rms_b1']))
    save(out/'RESULT.json',dict(complete=True,pic_rows=pic_rows,genesis_exits=exits,
        same_band_and_aperture=ratios,controls=controls,
        initial_complex_bf_magnitude_matched=True,physical_acceptance=False,
        initial_noise_seeds=1,power_growth_metric='Unwindowed full-pulse ±5% radiation band at matched ±560um aperture for p1 and p6; continuous GENESIS curve has full transverse aperture and full envelope bandwidth.',
        source_group_metric='Fixed source-ID groups of ten original optical cells, not rebinned by shifted arrival times. PIC transverse losses reported per plane.'))
    render(out,pic_rows,series,spectra)


def render(out,pic,series,spectra):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    fig,axes=plt.subplots(1,3,figsize=(16,4.9),layout='constrained')
    colors={'sc-off-fixed':'#337ba6','longrange-fixed':'#b8792e','longrange-installed':'#a7537d','sc-off-halfstep':'#657a48'}
    for name,s in series.items():
        style='--' if name in ['longrange-installed','sc-off-halfstep'] else '-'
        axes[0].semilogy(s['z_m'],s['source_group_local_rms_b1'],style,color=colors[name],label=name)
        axes[1].semilogy(s['z_m'],s['pulse_energy_J']*1e6,style,color=colors[name])
        axes[2].plot(s['z_m'],s['current_weighted_sigma_x_m']*1e6,style,color=colors[name],label=name+' x')
    axes[0].semilogy([x['z_m'] for x in pic],[x['source_group_local_rms_b1'] for x in pic],'o-',color='#333333',label='APRL F3: same input')
    axes[1].semilogy([x['z_m'] for x in pic[1:]],[x['rect_5pct_flux_J']*1e6 for x in pic[1:]],'o-',color='#333333')
    for name in ['sc-off-fixed','longrange-fixed']:
        rows=[r for r in spectra if r['case']==name]
        axes[1].semilogy([pic[r['module']]['z_m'] for r in rows],[r['rect_5pct_aperture_J']*1e6 for r in rows],'s',mfc='white',color=colors[name],markersize=7)
    for ax in axes:ax.set_xlabel('Laboratory z (m)');ax.grid(alpha=.2);ax.set_xlim(0,15.1)
    axes[0].set_ylabel('Local RMS BF, fixed original 10-cell groups');axes[0].set_title('(a) Same-input bunching');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Full-pulse energy (µJ)');axes[1].set_title('(b) Pulse growth; squares: matched band/aperture')
    axes[2].set_ylabel('Current-weighted horizontal rms size (µm)');axes[2].set_title('(c) GENESIS beam-envelope controls');axes[2].legend(fontsize=8)
    fig.suptitle('SPARC six-segment gain audit: unchanged 1,734,656 input macroparticles')
    fig.savefig(out/'same_input_comparison.png',dpi=180);plt.close(fig)
