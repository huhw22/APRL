"""Only after normal closure: gain, local BF and direct node-field spectra."""
import csv
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import C,E,read,save,need,dump_at


def table(path,rows):
    with Path(path).open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def particle_statistics(path,plan,info):
    with h5py.File(path) as f:
        a=f['electrons'][:];actual_z=.01+float(f['runInfo'].attrs['zbarTotal'])*float(f['runInfo'].attrs['Lg'])
    need(len(a)==plan['sample_records'] and np.isfinite(a).all(),'Particle population/nonfinite changed')
    sc=info['scale'];q=a[:,6]*sc['npkbar']*sc['qe_C']
    need(abs(q.sum()/info['total_charge_C']-1)<2e-12,'Particle charge changed')
    tau=info['arrival_origin_s']+a[:,2]*sc['cooperation_length_m']/sc['code_c_m_per_s']
    f0=C/plan['wavelength_m'];phase=np.exp(2j*np.pi*f0*tau)
    lo,hi=plan['steady_retarded_window_s'];use=(tau>=lo)&(tau<hi)
    need(use.any(),'Empty central particle window')
    group=np.floor((tau[use]-lo)/(10/f0)).astype(int);weight=q[use]
    counts=np.bincount(group,weights=weight)
    re=np.bincount(group,weights=weight*phase[use].real);im=np.bincount(group,weights=weight*phase[use].imag)
    occupied=counts>0
    local2=(re[occupied]**2+im[occupied]**2)/counts[occupied]**2
    return dict(actual_lab_z_m=actual_z,records=len(a),charge_C=float(q.sum()),
        global_b1=float(abs(np.sum(q*phase))/q.sum()),
        central_Eulerian_b1=float(abs(np.sum(weight*phase[use]))/weight.sum()),
        central_Eulerian_10lambda_rms_b1=float(np.sqrt(np.average(local2,weights=counts[occupied]))),
        central_records=int(use.sum()),central_charge_C=float(weight.sum()),
        rms_x_m=float(np.std(a[:,0])*sc['transverse_scale_m']),
        rms_y_m=float(np.std(a[:,1])*sc['transverse_scale_m']))


def field_statistics(path,plan,info,out,tag):
    sc=info['scale'];half=plan['aperture_half_width_m'];f0=C/plan['wavelength_m']
    with h5py.File(path) as f:
        attrs=f['runInfo'].attrs;g=f['aperp'];nz=g.shape[1];n=g.shape[2]
        dx=float(attrs['sLengthOfElmX'])*sc['transverse_scale_m']
        dy=float(attrs['sLengthOfElmY'])*sc['transverse_scale_m']
        dt=float(attrs['sLengthOfElmZ2'])*float(attrs['Lc'])/sc['code_c_m_per_s']
        x=(np.arange(n)-(n-1)/2)*dx;y=(np.arange(n)-(n-1)/2)*dy
        ix=np.flatnonzero(abs(x)<=half+1e-12);iy=np.flatnonzero(abs(y)<=half+1e-12)
        a=np.empty((nz,len(iy),len(ix),2))
        for j in range(2): a[:,:,:,j]=g[j,:,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
    need(np.isfinite(a).all(),'Nonfinite closed transverse field')
    wx=np.maximum(0,np.minimum(x[ix]+dx/2,half)-np.maximum(x[ix]-dx/2,-half))/dx
    wy=np.maximum(0,np.minimum(y[iy]+dy/2,half)-np.maximum(y[iy]-dy/2,-half))/dy
    area=wy[:,None]*wx[None,:]*dx*dy
    norm=sc['eps0']*sc['code_c_m_per_s']*sc['field_scale_V_per_m']**2
    tau=info['arrival_origin_s']+np.arange(nz)*dt
    masks={'full_record':np.ones(nz,dtype=bool),
           'central_window':(tau>=plan['steady_retarded_window_s'][0])&(tau<plan['steady_retarded_window_s'][1])}
    rows=[];saved={}
    for window,use in masks.items():
        b=a[use];N=len(b);duration=N*dt;frequency=np.fft.rfftfreq(N,dt)
        select=(frequency>=f0*.95)&(frequency<=f0*1.05);df=1/duration
        need(select.sum()>=20,'Insufficient resonant frequency resolution')
        for taper in ['Rect','Hann']:
            w=np.ones(N) if taper=='Rect' else np.hanning(N)
            # Fixed energy-window convention, not a fit to another solver.
            w=w/np.sqrt(np.mean(w*w))
            spectrum=np.fft.rfft(b*w[:,None,None,None],axis=0)*dt
            spectral=2*norm*np.sum(abs(spectrum[select])**2*area[None,:,:,None],axis=(1,2,3))
            energy=float(spectral.sum()*df)
            if taper=='Rect':
                factors=np.full(len(frequency),2.);factors[0]=1.
                if N%2==0:factors[-1]=1.
                fft_total=float(norm*np.sum(abs(spectrum)**2*area[None,:,:,None]*factors[:,None,None,None])*df)
                direct=float(norm*np.sum(b*b*area[None,:,:,None])*dt)
                error=abs(fft_total-direct)/max(abs(direct),1e-300)
                need(error<3e-12,'Real-field Parseval failed')
            else:error=0.
            rows.append(dict(window=window,taper=taper,band_energy_J=energy,
                record_mean_band_power_W=energy/duration,duration_s=duration,frequency_bins=int(select.sum()),
                carrier_spectral_density_J_per_Hz=float(np.interp(f0,frequency[select],spectral)),
                real_field_Parseval_relative=error))
            key=window+'_'+taper;saved[key+'_frequency_Hz']=frequency[select];saved[key+'_density_J_per_Hz']=spectral
            del spectrum
    np.savez_compressed(out/(tag+'-spectra.npz'),**saved)
    return rows


def history(root,plan,info,out):
    chosen={}
    for path in (root/'main').glob('main_integrated_*.h5'):
        with h5py.File(path) as f: step=int(f['runInfo'].attrs['iCsteps'])
        suffix=int(path.stem.rsplit('_',1)[1])
        if step not in chosen or suffix<chosen[step][0]:chosen[step]=(suffix,path)
    rows=[];power_profiles=[];b_profiles=[];sc=info['scale']
    for step,(_,path) in sorted(chosen.items()):
        with h5py.File(path) as f:
            attrs=f['runInfo'].attrs;z=.01+float(attrs['zbarTotal'])*float(attrs['Lg'])
            power=f['powerSI'][:];dt=float(attrs['sLengthOfElmZ2'])*float(attrs['Lc'])/sc['code_c_m_per_s']
            t=info['arrival_origin_s']+np.arange(len(power))*dt
            use=(t>=plan['steady_retarded_window_s'][0])&(t<plan['steady_retarded_window_s'][1])
            bf=f['bunchingFundamentalSI'][:];q=abs(f['Slice Charge SI'][:])
            tb=info['arrival_origin_s']+(np.arange(len(bf))+.5)*sc['lambda_reference_m']/sc['code_c_m_per_s']
            core=(tb>=plan['steady_retarded_window_s'][0])&(tb<plan['steady_retarded_window_s'][1])&(q>0)
            rows.append(dict(lab_z_m=z,step=step,central_mean_native_allband_power_W=float(np.mean(power[use])),
                central_peak_native_allband_power_W=float(power[use].max()),
                central_native_1lambda_rms_b1=float(np.sqrt(np.average(bf[core]**2,weights=q[core]))) if core.any() else 0.))
            if step==0 or step in {x['step'] for x in plan['outputs']}:
                power_profiles.append((z,t.copy(),power));b_profiles.append((z,tb.copy(),bf))
    table(out/'native_gain_history.csv',rows)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for z,t,power in power_profiles:ax[0].plot(t*1e12,power,label=f'{z:.2f}m')
    ax[0].set(xlabel='Entrance retarded coordinate (ps)',ylabel='Native all-grid power (W)',title='Pulse profiles, native power')
    ax[0].axvspan(*np.array(plan['steady_retarded_window_s'])*1e12,color='grey',alpha=.12)
    ax[0].legend(fontsize=6,ncol=2)
    for z,t,b in b_profiles:ax[1].plot(t*1e12,b,label=f'{z:.2f}m')
    ax[1].set(xlabel='Entrance retarded coordinate (ps)',ylabel='Native 1-wavelength |b1|',title='Local bunching profiles')
    fig.savefig(out/'pulse_profiles.png',dpi=160);plt.close(fig)
    return rows


def main(root):
    root=Path(root);need(read(root/'TRACKING_CLOSED.json')['complete'],'Refuse running H5')
    plan=read(root/'PLAN.json');info=read(root/'INPUT_VALIDATED.json');out=root/'analysis';out.mkdir()
    dense=history(root,plan,info,out);rows=[];all_fields=[]
    for item in plan['outputs']:
        tag=f"segment-{item['segment']:02d}";step=item['step']
        particle=particle_statistics(dump_at(root/'main','electrons',step),plan,info)
        fields=field_statistics(dump_at(root/'main','aperp',step),plan,info,out,tag)
        for row in fields:row.update(segment=item['segment'],lab_z_m=particle['actual_lab_z_m'])
        all_fields.extend(fields)
        central=next(x for x in fields if x['window']=='central_window' and x['taper']=='Hann')
        full=next(x for x in fields if x['window']=='full_record' and x['taper']=='Hann')
        rows.append(dict(segment=item['segment'],**particle,central_Hann_band_energy_J=central['band_energy_J'],
            central_Hann_band_power_W=central['record_mean_band_power_W'],full_Hann_band_energy_J=full['band_energy_J']))
        print('ANALYZED',tag,rows[-1],flush=True)
    table(out/'gain.csv',rows);table(out/'band_statistics.csv',all_fields)
    fig,ax=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
    z=[x['actual_lab_z_m'] for x in rows]
    ax[0].semilogy([x['lab_z_m'] for x in dense],[max(x['central_mean_native_allband_power_W'],1e-30) for x in dense],label='Native all-grid/all-band')
    ax[0].semilogy(z,[x['central_Hann_band_power_W'] for x in rows],'o-',label='Hann ±5%, ±560µm')
    ax[0].set(xlabel='Lab z (m)',ylabel='Central-window mean power (W)');ax[0].legend(fontsize=8)
    ax[1].semilogy(z,[x['central_Hann_band_energy_J']*1e9 for x in rows],'o-',label='Central window')
    ax[1].semilogy(z,[x['full_Hann_band_energy_J']*1e9 for x in rows],'s-',label='Full record')
    ax[1].set(xlabel='Lab z (m)',ylabel='Hann ±5% energy (nJ)');ax[1].legend(fontsize=8)
    ax[2].plot(z,[x['central_Eulerian_10lambda_rms_b1'] for x in rows],'o-',label='Central 10λ RMS')
    ax[2].plot(z,[x['global_b1'] for x in rows],'s-',label='Full-pulse global')
    ax[2].set(xlabel='Lab z (m)',ylabel='|b1|');ax[2].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.suptitle('Puffin | independent w≈64 | eleven segments | bunch length = 4 nominal total slips')
    fig.supxlabel('Central statistics use a fixed retarded-time window; all-grid native power and aperture/band power are different observables.',fontsize=9)
    fig.savefig(out/'gain.png',dpi=170);plt.close(fig)
    save(out/'RESULT.json',dict(complete=True,physical_acceptance=False,segments=rows,
        steady_window_s=plan['steady_retarded_window_s'],scope='Puffin only; no amplitude fitting; node intensity, no added transverse averaging.',
        limitations=['One weighted sample; no physical one-electron SASE noise acceptance.',
            'Central Eulerian bunching is not a fixed-ID cohort.','Actual output planes are the next integration step after nominal exits.',
            'Puffin remains forward/paraxial without space charge.','A four-slip parabola does not guarantee a flat temporal plateau.']))


if __name__=='__main__':main(sys.argv[1])
