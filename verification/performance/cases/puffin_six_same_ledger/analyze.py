"""Closed-only lab particle losses and retained Puffin radiation energies."""
import csv
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import C,E,read,save,need,dump_at

LD=np.longdouble
M=LD('9.1093837139e-31')  # Same electron mass as the frozen APRL production code.


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def particles(path,plan,info):
    with h5py.File(path) as f:
        a=f['electrons'][:];attrs=f['runInfo'].attrs
        z=.01+float(attrs['zbarTotal'])*float(attrs['Lg'])
    need(len(a)==plan['sample_records'] and np.isfinite(a).all(),'Particle population/finite check failed')
    sc=info['scale'];q=a[:,6].astype(LD)*LD(sc['npkbar'])*LD(sc['qe_C']);charge=q.sum(dtype=LD)
    need(abs(charge/LD(info['total_charge_C'])-1)<LD('2e-12'),'Charge changed')
    g=a[:,5].astype(LD)*LD(plan['gamma_ref']);n=q/LD(E)
    delta_g=(g-LD(info['physical_gamma_initial_mean']))
    delta=LD(M*C*C)*np.sum(n*delta_g,dtype=LD)
    mean_eV=delta/charge
    tau=LD(info['arrival_origin_s'])+a[:,2].astype(LD)*LD(sc['cooperation_length_m'])/LD(sc['code_c_m_per_s'])
    phase=np.exp(2j*np.pi*np.asarray(tau,dtype=float)*C/plan['wavelength_m'])
    return dict(actual_lab_z_m=z,records=len(a),charge_C=float(charge),
        particle_delta_K_lab_J=float(delta),mean_delta_K_eV_per_electron=float(mean_eV),
        mean_gamma_lab=float(np.sum(q*g,dtype=LD)/charge),
        global_b1=float(abs(np.sum(q*phase))/charge),
        rest_energy_definition='APRL physical SI constants; physical charge-weighted mean gamma change')


def fields(path,plan,info,out,tag):
    sc=info['scale'];norm=LD(sc['eps0'])*LD(sc['code_c_m_per_s'])*LD(sc['field_scale_V_per_m'])**2
    with h5py.File(path) as f:
        g=f['aperp'];attrs=f['runInfo'].attrs;need(g.shape[0]==2,'Unknown Puffin field layout')
        nz,ny,nx=g.shape[1:];need(nx==ny==81,'Unexpected transverse grid')
        dx=float(attrs['sLengthOfElmX'])*sc['transverse_scale_m'];dy=float(attrs['sLengthOfElmY'])*sc['transverse_scale_m']
        dt=float(attrs['sLengthOfElmZ2'])*float(attrs['Lc'])/sc['code_c_m_per_s']
        x=(np.arange(nx)-(nx-1)/2)*dx;y=(np.arange(ny)-(ny-1)/2)*dy
        wx=np.ones(nx);wy=np.ones(ny);wx[[0,-1]]=.5;wy[[0,-1]]=.5
        area=wy[:,None]*wx[None,:]*dx*dy;full=LD(0);edge=LD(0);front=LD(0);back=LD(0)
        # Native fPower_3D integrates nodal intensity with m_trapz2D.
        for lo in range(0,nz,256):
            a=g[:,lo:min(lo+256,nz),:,:];need(np.isfinite(a).all(),'Nonfinite full-grid field')
            t=np.sum(a*a*area[None,None,:,:],axis=(0,2,3),dtype=LD)
            full+=t.sum(dtype=LD)
            if lo==0:front+=t[:40].sum(dtype=LD)
            if lo+256>=nz:back+=t[-40:].sum(dtype=LD)
            edge+=np.sum(a[:,:,:2,:]**2,dtype=LD)+np.sum(a[:,:,-2:,:]**2,dtype=LD)+np.sum(a[:,:,:,:2]**2,dtype=LD)+np.sum(a[:,:,:,-2:]**2,dtype=LD)
        full_J=float(full*norm*LD(dt));rows=[];cache={}
        for half in plan['apertures_half_width_m']:
            ix=np.flatnonzero(abs(x)<=half+1e-12);iy=np.flatnonzero(abs(y)<=half+1e-12)
            b=np.empty((nz,len(iy),len(ix),2))
            for k in range(2):b[:,:,:,k]=g[k,:,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
            ax=np.maximum(0,np.minimum(x[ix]+dx/2,half)-np.maximum(x[ix]-dx/2,-half))
            ay=np.maximum(0,np.minimum(y[iy]+dy/2,half)-np.maximum(y[iy]-dy/2,-half))
            aa=ay[:,None]*ax[None,:];duration=nz*dt;freq=np.fft.rfftfreq(nz,dt);df=1/duration
            direct=float(norm*np.sum(b*b*aa[None,:,:,None],dtype=LD)*LD(dt))
            f0=C/plan['wavelength_m'];band=(freq>=f0*.95)&(freq<=f0*1.05)
            for taper in ['Rect','Hann']:
                window=np.ones(nz) if taper=='Rect' else np.hanning(nz)
                window/=np.sqrt(np.mean(window*window))
                fft=np.fft.rfft(b*window[:,None,None,None],axis=0)*dt
                factors=np.full(len(freq),2.);factors[0]=1.
                if nz%2==0:factors[-1]=1.
                band_density=2*float(norm)*np.sum(abs(fft[band])**2*aa[None,:,:,None],axis=(1,2,3))
                band_J=float(np.sum(band_density)*df)
                if taper=='Rect':
                    total=float(norm)*float(np.sum(abs(fft)**2*aa[None,:,:,None]*factors[:,None,None,None],dtype=LD))*df
                    error=abs(total-direct)/max(abs(direct),1e-300);need(error<3e-12,'Parseval failed')
                else:error=0.
                rows.append(dict(aperture_half_width_m=half,taper=taper,record_duration_s=duration,
                    all_retained_field_energy_J=direct,resonant_pm5pct_energy_J=band_J,
                    record_mean_band_power_W=band_J/duration,Parseval_relative_error=error))
                key=f'{round(half*1e6)}um_{taper}';cache[key+'_frequency_Hz']=freq[band];cache[key+'_density_J_per_Hz']=band_density
                del fft
    np.savez_compressed(out/(tag+'-spectra.npz'),**cache)
    return rows,dict(all_grid_retained_forward_field_energy_J=full_J,dt_s=dt,field_nodes_z2=nz,
        first40_time_nodes_energy_fraction=float(front/full) if full else 0.,
        last40_time_nodes_energy_fraction=float(back/full) if full else 0.,
        full_grid_spatial_quadrature='Native nodal trapezoid, two physical transverse polarizations',
        transverse_edge_unweighted_square_sum=float(edge))


def integrated(root,info,out):
    chosen={}
    for p in (root/'main').glob('main_integrated_*.h5'):
        with h5py.File(p) as f:step=int(f['runInfo'].attrs['iCsteps'])
        number=int(p.stem.rsplit('_',1)[1])
        if step not in chosen or number<chosen[step][0]:chosen[step]=(number,p)
    rows=[]
    for step,(_,p) in sorted(chosen.items()):
        with h5py.File(p) as f:
            a=f['runInfo'].attrs;dt=float(a['sLengthOfElmZ2'])*float(a['Lc'])/info['scale']['code_c_m_per_s']
            power=f['powerSI'][:];z=.01+float(a['zbarTotal'])*float(a['Lg'])
            need(np.isfinite(power).all(),'Nonfinite integrated field power')
            rows.append(dict(step=step,lab_z_m=z,all_grid_native_field_energy_J=float(np.sum(power,dtype=LD)*LD(dt)),
                peak_power_W=float(power.max()),record_duration_s=len(power)*dt))
    need(rows and rows[0]['step']==0,'Missing initial integrated field')
    table(out/'native_energy_history.csv',rows);return rows


def main(root):
    root=Path(root);need(read(root/'TRACKING_CLOSED.json')['complete'],'Refuse running H5')
    plan=read(root/'PLAN.json');info=read(root/'INPUT_VALIDATED.json');out=root/'analysis';out.mkdir()
    dense=integrated(root,info,out);rows=[];bands=[];checks={}
    # Initial full-field check also confirms the absence of an imported seed.
    initial_bands,initial=fields(dump_at(root/'main','aperp',0),plan,info,out,'initial')
    initial_particle=particles(dump_at(root/'main','electrons',0),plan,info)
    need(abs(initial_particle['mean_delta_K_eV_per_electron'])<1e-5,'Physical initial particle energy mismatch')
    for item in plan['outputs']:
        step=item['step'];tag=f"segment-{item['segment']:02d}"
        p=particles(dump_at(root/'main','electrons',step),plan,info)
        band,field=fields(dump_at(root/'main','aperp',step),plan,info,out,tag)
        for b in band:b.update(segment=item['segment'],actual_lab_z_m=p['actual_lab_z_m'])
        bands.extend(band);field_gain=field['all_grid_retained_forward_field_energy_J']-initial['all_grid_retained_forward_field_energy_J']
        row=dict(segment=item['segment'],nominal_lab_z_m=item['nominal_lab_z_m'],**p,
            retained_forward_field_gain_all_grid_J=field_gain,
            particle_plus_retained_field_change_J=p['particle_delta_K_lab_J']+field_gain)
        for half in [520,560]:
            for taper in ['Rect','Hann']:
                b=next(b for b in band if abs(b['aperture_half_width_m']-half*1e-6)<1e-12 and b['taper']==taper)
                row[f'field_energy_{half}um_J']=b['all_retained_field_energy_J']
                row[f'{taper}_resonant_energy_{half}um_J']=b['resonant_pm5pct_energy_J']
        rows.append(row);checks[tag]=field
        native=next((a for a in dense if a['step']==step),None)
        if native is not None:
            error=abs(native['all_grid_native_field_energy_J']-field['all_grid_retained_forward_field_energy_J'])/max(abs(field['all_grid_retained_forward_field_energy_J']),1e-300)
            need(error<2e-11,'Independent full-grid integral disagrees with native power')
            checks[tag]['native_power_integration_relative_error']=error
        print('SEGMENT_ENERGY',row,flush=True)
    table(out/'energy_gain.csv',rows);table(out/'band_statistics.csv',bands)
    fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout='constrained');z=[r['actual_lab_z_m'] for r in rows]
    axes[0].plot(z,[-r['mean_delta_K_eV_per_electron']/1000 for r in rows],'o-',label='Puffin mean loss')
    axes[0].set(xlabel='Lab z (m)',ylabel='Physical electron mean loss (keV)')
    axes[1].plot([r['lab_z_m'] for r in dense],[r['all_grid_native_field_energy_J']*1e6 for r in dense],label='All grid, all retained modes')
    for key,label in [('field_energy_520um_J','Aperture ±520µm'),('Hann_resonant_energy_520um_J','Hann ±5%, aperture ±520µm')]:
        axes[1].plot(z,[r[key]*1e6 for r in rows],'o-',label=label)
    axes[1].set(xlabel='Lab z (m)',ylabel='Forward field energy (µJ)')
    axes[2].plot(z,[r['particle_plus_retained_field_change_J']*1e6 for r in rows],'o-',label='Particle + retained forward field')
    axes[2].set(xlabel='Lab z (m)',ylabel='Uncorrected forward-model sum (µJ)')
    comparison=[];reference=root/'reference';ref=read(reference/'RESULT.json')['cases']
    for name,color in [('energy-binomial3','#888888'),('energy-momentum','#aa8866')]:
        r=ref[name];loss=-r['lab_particle_change_full_cohort_J']/plan['source_charge_C']/1000
        axes[0].scatter([plan['outputs'][-1]['nominal_lab_z_m']],[loss],label='APRL '+name,color=color,marker='s')
        with (reference/'lab-v2'/name/'lab_energy_absolute.csv').open() as f:last=list(csv.DictReader(f))[-1]
        comparison.append(dict(method='APRL '+name,particle_delta_K_lab_J=r['lab_particle_change_full_cohort_J'],
            mean_delta_eV=r['lab_particle_change_full_cohort_J']/plan['source_charge_C'],
            exit_aperture520_field_J=float(last['field_exit_signed_J']),scope='PIC full field including near field; not Puffin retained forward modes'))
    last=rows[-1];comparison.append(dict(method='Puffin historical same source',particle_delta_K_lab_J=last['particle_delta_K_lab_J'],
        mean_delta_eV=last['mean_delta_K_eV_per_electron'],exit_aperture520_field_J=last['field_energy_520um_J'],
        scope='Forward paraxial, no space charge or velocity near field; native low-frequency filter retained'))
    table(out/'final_comparison.csv',comparison)
    for a in axes:a.grid(alpha=.2);a.legend(fontsize=7);a.axhline(0,color='grey',lw=.6)
    fig.suptitle('Six-segment Puffin | exactly the historical APRL input and macro charges\nFull pulse; no amplitude fit; field model and unmeasured side/filter exchanges limit interpretation of the sum')
    fig.savefig(out/'energy_comparison.png',dpi=170);fig.savefig(out/'energy_comparison.pdf');plt.close(fig)
    save(out/'RESULT.json',dict(analysis_complete=True,physical_acceptance=False,segments=rows,field_checks=checks,
        initial_forward_field=initial,initial_particle=initial_particle,comparison=comparison,
        limits=['Historical correlated input; not a new independent-random benchmark.',
            'Puffin is forward/paraxial, without PIC longitudinal self-field or velocity near field.',
            'The native filter removes low frequencies; retained-field sum is not a full physical PIC worldtube balance.',
            'No cumulative transverse field escape or filter-loss channel is recorded; their unknown contributions are not filled from the residual.',
            'Original events drift ballistically by1cm to the Puffin entrance; PIC can have self-fields there.',
            'Persistent particle IDs are unavailable after Puffin redistribution; weighted full-cohort means are used.',
            'Output planes follow nominal exits by at most0.35mm; exact coordinates are stored.']))


if __name__=='__main__':main(sys.argv[1])
