"""Closed six-segment power/BF comparison with the same Eulerian definition."""
import csv
import json
import math
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy import fft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import C,read,save,need,dump_at,frozen,sha

LD=np.longdouble
EPS0=8.8541878128e-12


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def particle_stats(path,kind,info,plan,low,high):
    with h5py.File(path,'r') as f:
        if kind=='PIC':
            g=f['particle_substep_plane'];n=int(g['committed_records'][()]);need(int(g['complete'][()])==1 and n==len(g['records']),'Unclosed PIC particle plane')
            a=g['records'].fields(['time_s','charge_C'])[:];z=float(g.attrs['plane_z_m'])
            tau=a['time_s'].astype(LD)-LD(z)/LD(C);q=abs(a['charge_C']).astype(LD)
        else:
            a=f['electrons'][:];attrs=f['runInfo'].attrs;z=.01+float(attrs['zbarTotal'])*float(attrs['Lg'])
            sc=info['scale'];tau=LD(info['arrival_origin_s'])+a[:,2].astype(LD)*LD(sc['cooperation_length_m'])/LD(sc['code_c_m_per_s'])-LD('.01')/LD(C)
            q=a[:,6].astype(LD)*LD(sc['npkbar'])*LD(sc['qe_C']);n=len(a)
    need(np.isfinite(tau).all() and np.isfinite(q).all() and np.all(q>0),'Invalid particle times or charges')
    order=np.argsort(tau);s=np.cumsum(q[order],dtype=LD);median=tau[order[np.searchsorted(s,s[-1]/2)]]
    rel=tau-median;use=(rel>=low)&(rel<high);need(use.any(),'Empty central particle window')
    f0=LD(C)/LD(plan['wavelength_m']);cycle=np.remainder(rel*f0,LD(1));phase=np.exp(2j*np.pi*np.asarray(cycle,dtype=float))
    global_b=float(abs(np.sum(q*phase))/q.sum())
    central_b=float(abs(np.sum(q[use]*phase[use]))/q[use].sum())
    group=np.floor((rel[use]-low)*f0/10).astype(int)
    qq=np.bincount(group,weights=np.asarray(q[use],float));re=np.bincount(group,weights=np.asarray(q[use]*phase[use].real,float));im=np.bincount(group,weights=np.asarray(q[use]*phase[use].imag,float))
    good=qq>0;local=float(np.sqrt(np.sum((re[good]**2+im[good]**2)/qq[good])/np.sum(qq[good])))
    need(0<=global_b<=1 and central_b<=local+1e-12 and local<=1,'Invalid BF or RMS inequality')
    return dict(actual_z_m=z,records=n,charge_C=float(q.sum()),charge_fraction_of_original=float(q.sum()/LD(plan['source_charge_C'])),
        median_retarded_time_s=str(median),central_records=int(use.sum()),central_charge_C=float(q[use].sum()),
        global_b1=global_b,central_b1=central_b,central_10lambda_rms_b1=local,central_bins_10lambda=int(good.sum()))


def field_power(path,kind,info,plan,median,low,high):
    half=.00052;width=float(high-low)
    with h5py.File(path,'r') as f:
        if kind=='PIC':
            g=f['field_plane'];need(int(g['complete'][()])==1,'Unclosed PIC field')
            attrs=g.attrs;n=int(g['committed_samples'][()]);need(n==len(g['fields'])==len(g['time_s']),'PIC field length mismatch')
            z=LD(attrs['plane_z_m']);t=g['time_s'][:].astype(LD)-z/LD(C)
            dx=float(attrs['dx_m']);dy=float(attrs['dy_m']);ny,nx=g['fields'].shape[1:]
            x=float(attrs['x_first_m'])+np.arange(nx)*dx;y=float(attrs['y_first_m'])+np.arange(ny)*dy
        else:
            g=f['aperp'];attrs=f['runInfo'].attrs;sc=info['scale'];n,ny,nx=g.shape[1:]
            dx=float(attrs['sLengthOfElmX'])*sc['transverse_scale_m'];dy=float(attrs['sLengthOfElmY'])*sc['transverse_scale_m']
            dt=float(attrs['sLengthOfElmZ2'])*float(attrs['Lc'])/sc['code_c_m_per_s']
            t=LD(info['arrival_origin_s'])+np.arange(n,dtype=LD)*LD(dt)-LD('.01')/LD(C)
            x=(np.arange(nx)-(nx-1)/2)*dx;y=(np.arange(ny)-(ny-1)/2)*dy
        rel=t-median;indices=np.flatnonzero((rel>=low)&(rel<high));need(len(indices)>500,'Insufficient time coverage')
        need(rel[0]<=low and rel[-1]>=high,'Requested window extends outside measured field')
        first,last=int(indices[0]),int(indices[-1])+1
        ax=np.maximum(0,np.minimum(x+dx/2,half)-np.maximum(x-dx/2,-half))
        ay=np.maximum(0,np.minimum(y+dy/2,half)-np.maximum(y-dy/2,-half))
        ix=np.flatnonzero(ax>dx*1e-8);iy=np.flatnonzero(ay>dy*1e-8);area=ay[iy,None]*ax[None,ix]
        if kind=='PIC':
            b=g['fields'][first:last,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]['electric_V_per_m'][...,:2].copy()
        else:
            b=np.empty((last-first,len(iy),len(ix),2))
            for k in range(2):b[:,:,:,k]=g[k,first:last,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]*info['scale']['field_scale_V_per_m']
    rel=rel[first:last];d=np.diff(rel);dt=float(np.median(d));ticks=np.rint((rel-rel[0])/LD(dt)).astype(int)
    need(np.all(np.diff(ticks)>0),'Duplicated field times')
    missing=np.setdiff1d(np.arange(ticks[-1]+1),ticks)
    need(len(missing)<=1,'Too many gaps for the existing interpolation convention')
    if len(missing):
        full=np.empty((ticks[-1]+1,*b.shape[1:]));full[ticks]=b
        for k in missing:
            nearest=np.argsort(abs(ticks-k))[:8];nodes=ticks[nearest].astype(float);w=np.ones(8)
            for j in range(8):w[j]=np.prod((k-np.delete(nodes,j))/(nodes[j]-np.delete(nodes,j)))
            full[k]=np.tensordot(w,b[nearest],axes=(0,0))
        b=full;rel=rel[0]+np.arange(len(b),dtype=LD)*LD(dt)
    need(np.isfinite(b).all(),'Nonfinite window field')
    # The same continuous gate and RMS-Hann amplitude for every solver.
    # Translation follows the charge median, never the observed field peak.
    window=np.sqrt(8/3)*np.sin(np.pi*np.asarray((rel-low)/(high-low),float))**2
    weighted=b*window[:,None,None,None];N=len(weighted)
    spec=fft.rfft(weighted,axis=0,workers=4)*dt;freq=fft.rfftfreq(N,dt);df=1/(N*dt)
    f0=C/plan['wavelength_m'];use=(freq>=.95*f0)&(freq<=1.05*f0);need(use.sum()>=20,'Unresolved ±5% band')
    factors=np.full(len(freq),2.);factors[0]=1.
    if N%2==0:factors[-1]=1.
    density=EPS0*C*np.sum(abs(spec)**2*area[None,:,:,None],axis=(1,2,3))*factors
    total=float(density.sum()*df);direct=float(EPS0*C*np.sum(weighted*weighted*area[None,:,:,None],dtype=LD)*dt)
    error=abs(total-direct)/max(abs(direct),1e-300);need(error<3e-12,'Hann power Parseval mismatch')
    energy=float(density[use].sum()*df)
    return dict(Hann_pm5pct_power_Eproxy_W=energy/width,Hann_pm5pct_energy_J=energy,
        power_denominator_s=width,FFT_samples=N,native_dt_s=dt,actual_sample_duration_s=N*dt,
        missing_samples_reconstructed=int(len(missing)),band_bins=int(use.sum()),Parseval_relative_error=error,
        aperture_area_m2=float(area.sum()),Hann_weight_convention='sqrt(8/3)*sin²(pi*(t-low)/fixed_width)',
        field_window_relative_low_s=float(low),field_window_relative_high_s=float(high))


def main(root,out):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=False);plan=read(root/'PLAN.json');info=read(root/'INPUT_VALIDATED.json')
    need(read(root/'TRACKING_CLOSED.json')['complete'],'Refuse unclosed Puffin');frozen(root/'FROZEN.json')
    aprl=Path(plan['APRL_parent']);ap=read(aprl/'PLAN.json');need(plan['input_sha256']==ap['input_sha256'],'Inputs do not match')
    half=LD(str(plan['duration_s']))/2;s=LD(str(math.sqrt(.2)));low=LD('.1')*s*half;high=LD('.9')*s*half
    rows=[];initial=[];protected={}
    for name,kind in [('Puffin','Puffin'),('energy-binomial3','PIC'),('energy-momentum','PIC')]:
        if kind=='PIC':need(read(aprl/name/'PIC_COMPLETE.json')['complete'],'Refuse incomplete PIC')
        for item in plan['outputs']:
            k=item['segment'];step=item['step']
            p=dump_at(root/'main','electrons',step) if kind=='Puffin' else aprl/name/'detectors'/f'p{k}-particles-substeps.h5'
            f=dump_at(root/'main','aperp',step) if kind=='Puffin' else aprl/name/'detectors'/f'p{k}.h5'
            for path in [p,f]:
                st=path.stat();protected[str(path)]=dict(size=st.st_size,mtime_ns=st.st_mtime_ns,inode=st.st_ino)
            ps=particle_stats(p,kind,info,plan,low,high);fp=field_power(f,kind,info,plan,LD(ps['median_retarded_time_s']),low,high)
            row=dict(method=name,segment=k,**ps,**fp);rows.append(row);print('GROWTH',row,flush=True)
        p=dump_at(root/'main','electrons',0) if kind=='Puffin' else aprl/name/'detectors/und-entry-particles-substeps.h5'
        initial.append(dict(method=name,segment=0,**particle_stats(p,kind,info,plan,low,high)))
    table(out/'growth.csv',rows);table(out/'initial_bunching.csv',initial)
    colors=['#2879aa','#d38a35','#6a9270'];labels=['Puffin','APRL: energy + three-point','APRL: energy + momentum'];methods=['Puffin','energy-binomial3','energy-momentum']
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    fig2,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    for name,color,label in zip(methods,colors,labels):
        rr=[x for x in rows if x['method']==name];bf=[x for x in initial if x['method']==name]+rr
        axes[0].semilogy([x['actual_z_m'] for x in rr],[x['Hann_pm5pct_power_Eproxy_W'] for x in rr],'o-',color=color,label=label)
        axes[1].semilogy([x['actual_z_m'] for x in bf],[max(x['central_10lambda_rms_b1'],1e-12) for x in bf],'o-',color=color,label=label)
        ax.semilogy([x['actual_z_m'] for x in bf],[max(x['global_b1'],1e-12) for x in bf],'o-',color=color,label=label)
    axes[0].set(xlabel='Lab z (m)',ylabel='Hann ±5% transverse-field power proxy (W)',title='Same central high-current window; aperture ±520µm')
    axes[1].set(xlabel='Lab z (m)',ylabel='Charge-weighted local RMS |b1| (10λ bins)',title='Same Eulerian window relative to beam charge median')
    for a in axes:a.grid(alpha=.2,which='both');a.legend(fontsize=8)
    fig.suptitle('Six-segment Puffin / APRL comparison | same historical source and macro charges\nLater-half high-current interval, fixed width, no field-peak selection or amplitude fitting')
    fig.savefig(out/'power_bf_growth.png',dpi=180);fig.savefig(out/'power_bf_growth.pdf');plt.close(fig)
    ax.set(xlabel='Lab z (m)',ylabel='Global full-available-beam |b1|',title='Global BF: separate observable from local RMS')
    ax.grid(alpha=.2,which='both');ax.legend();fig2.savefig(out/'global_bf_growth.png',dpi=180);fig2.savefig(out/'global_bf_growth.pdf');plt.close(fig2)
    for path,st in protected.items():
        s=Path(path).stat();need(st==dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino),'Original output changed')
    save(out/'RESULT.json',dict(complete=True,physical_acceptance=False,rows=rows,initial_bunching=initial,
        fixed_window_relative_to_charge_median_s=[float(low),float(high)],fixed_window_duration_s=float(high-low),
        field_definition='epsilon0*c*integral |E_transverse|²; same continuous RMS-Hann gate, physical ±5% band, ±520µm aperture; divide by that fixed gate duration.',
        BF_definition='Central Eulerian gate, charge-weighted RMS over10λ groups. Whole-available-beam global BF is reported separately; no persistent-ID claim.',
        window_rule='Later half of nominal parabolic profile, I/Ipk>=0.8, remove10% at each end of the selected interval; translate with charge median, not radiation amplitude.',
        protected_closed_outputs=len(protected),max_Parseval_relative=max(x['Parseval_relative_error'] for x in rows),
        limits=['Historical correlated input, not independent-random SASE acceptance.',
            'PIC uses its saved Yee-to-cell detector interpolation; Puffin uses native nodes with aperture quadrature, no added four-point averaging.',
            'PIC narrow-band E-field proxy can retain near-field content; not a pure-radiation separation.',
            'Eulerian central statistics are not fixed source-ID cohorts.',
            'Full-beam global BF uses available PIC crossings; charge coverage and incomplete p6 tails are reported.',
            'Power means divide by the common0.497ps gate, not the entire simulation or padding duration.']))


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
