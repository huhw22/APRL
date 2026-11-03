"""Read only closed PIC files. Do not infer radiation from full raw field energy."""
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c,e,h,m_e,epsilon_0,mu_0
from prepare import save,require

def read_particles(path,expected):
    with h5py.File(path,'r') as f:
        g=f['particle_plane']
        require(int(g['complete'][()])==1,'Unclosed particle file')
        n=int(g['committed_records'][()])
        require(n==len(g['records'])==expected,'Missing particle crossings')
        a=g['records'][:]
    for name in a.dtype.names:
        if a[name].dtype.kind in 'fc':require(np.isfinite(a[name]).all(),'Nonfinite particles')
    a=a[np.argsort(a['source_id'])]
    require(np.array_equal(a['source_id'],np.arange(1,expected+1)),'Mismatched particle IDs')
    return a

def interpolate_one(nodes,values,target):
    x=np.asarray(nodes,dtype=float)-float(target)
    weights=np.array([np.prod(-np.delete(x,j)/(x[j]-np.delete(x,j)))
                      for j in range(len(x))])
    return np.tensordot(weights,values,axes=(0,0))

def analyze(root,case,output_name='analysis'):
    root=Path(root);plan=json.loads((root/'campaign.json').read_text())
    d=Path(case['directory']);out=d/output_name;out.mkdir(exist_ok=False)
    pars=plan['parameters'];geom=plan['geometry']
    lam=pars['resonant_wavelength_m'];omega=2*np.pi*c/lam;photon=h*c/(e*lam)
    n=pars['macroparticles'];entry=read_particles(d/'detectors/entry.h5',n)
    u=entry['proper_velocity'].astype(np.longdouble);g0=np.sqrt(1+np.sum(u*u,axis=1))
    weights=abs(entry['charge_C']);weights/=weights.sum()
    particle=[]
    for name,z in [('entry',0.),*geom['field_planes']]:
        a=entry if name=='entry' else read_particles(d/'detectors'/(name+'-particles.h5'),n)
        u=a['proper_velocity'].astype(np.longdouble);gamma=np.sqrt(1+np.sum(u*u,axis=1))
        dg=gamma-g0;mean=float(np.sum(weights*dg));mean_g=float(np.sum(weights*gamma))
        phase=omega*(a['time_s']+geom['reference_time_offset_s'])
        bunch=np.sum(weights*np.exp(1j*phase))
        modulation=np.sum(weights*(gamma-mean_g)*np.exp(1j*phase))
        particle.append(dict(plane=name,z_m=z,count=n,bunching=float(abs(bunch)),
          bunching_complex=[float(bunch.real),float(bunch.imag)],
          mean_kinetic_change_eV=mean*m_e*c*c/e,
          rms_kinetic_change_eV=float(np.sqrt(np.sum(weights*(dg-mean)**2)))*m_e*c*c/e,
          harmonic_gamma_modulation=[float(modulation.real),float(modulation.imag)],
          sigma_x_m=float(np.sqrt(np.sum(weights*(a['position_m'][:,0]-np.sum(weights*a['position_m'][:,0]))**2))),
          sigma_y_m=float(np.sqrt(np.sum(weights*(a['position_m'][:,1]-np.sum(weights*a['position_m'][:,1]))**2)))))
    save(out/'particles.json',particle)
    # All planes use a common retarded-time interval fixed from the input and
    # the theoretical maximum slippage, never from the measured peak location.
    ref=geom['reference_time_offset_s']
    start=-ref
    finish=-ref+pars['pulse_cycles']*lam/c+pars['periods']*lam/c+.35/(2*pars['gamma']**2*c)+1e-15
    summaries=[];cache={}
    for name,z in geom['field_planes']:
        path=d/'detectors'/(name+'.h5')
        with h5py.File(path,'r') as f:
            g=f['field_plane'];attrs=dict(g.attrs)
            require(int(g['complete'][()])==1,'Unclosed field file')
            total=int(g['committed_samples'][()])
            require(total==len(g['fields'])==len(g['time_s']),'Uncommitted field samples')
            frame=attrs['frame'];frame=frame.decode() if isinstance(frame,bytes) else frame
            require(frame=='laboratory','Wrong field frame')
            t=g['time_s'][:]-z/c
            nominal=geom['expected_record_s']/geom['field_steps']
            ticks=np.rint((t-t[0])/nominal).astype(np.int64)
            require(np.all(np.diff(ticks)>0),'Repeated or reversed field times')
            dt=float((t[-1]-t[0])/(ticks[-1]-ticks[0]))
            jitter=float(np.max(abs(t-(t[0]+ticks*dt)))/dt)
            require(jitter<.005,'Times do not lie on the native step lattice')
            i0=int(np.ceil((start-t[0])/dt));nc=int(np.floor((finish-start)/dt))
            target=np.arange(i0,i0+nc,dtype=np.int64)
            require(i0>=ticks[0] and target[-1]<=ticks[-1],'Common time interval is unavailable')
            times=t[0]+target*dt
            nx,ny=int(attrs['nx']),int(attrs['ny']);dx,dy=float(attrs['dx_m']),float(attrs['dy_m'])
            x=float(attrs['x_first_m'])+np.arange(nx)*dx
            y=float(attrs['y_first_m'])+np.arange(ny)*dy
            radius=np.hypot(y[:,None],x[None,:])
            roi=(radius>=90e-6)&(radius<130e-6)
            require(roi.sum()>100,'Insufficient off-beam radiation pixels')
            lo=max(0,int(np.searchsorted(ticks,target[0]))-5)
            hi=min(total,int(np.searchsorted(ticks,target[-1]))+6)
            raw=np.empty((hi-lo,int(roi.sum()),4))
            for a in range(lo,hi,16):
                b=min(hi,a+16);data=g['fields'][a:b]
                E=data['electric_V_per_m'];B=data['magnetic_T']
                require(np.isfinite(E).all() and np.isfinite(B).all(),'Nonfinite plane')
                raw[a-lo:b-lo,:,:2]=E[:,roi,:2];raw[a-lo:b-lo,:,2:]=B[:,roi,:2]
            lookup={int(k):i for i,k in enumerate(ticks[lo:hi])}
            missing=[int(k) for k in target if int(k) not in lookup]
            require(len(missing)<=max(1,nc//1000),'Too many missing native field samples')
            cube=np.empty((nc,int(roi.sum()),4))
            for j,k in enumerate(target):
                if int(k) in lookup:
                    cube[j]=raw[lookup[int(k)]]
                else:
                    near=np.argsort(abs(ticks[lo:hi]-k))[:8]
                    require(len(near)==8,'Insufficient interpolation neighbors')
                    cube[j]=interpolate_one(ticks[lo:hi][near],raw[near],k)
        energy=h*np.fft.rfftfreq(nc,dt)/e;de=energy[1]-energy[0]
        stats={};cv=[]
        for label,win in [('rectangular',np.ones(nc)),('Hann',np.hanning(nc))]:
            norm=float(np.mean(win**2))
            F=np.fft.rfft(cube*win[:,None,None],axis=0)*dt/np.sqrt(norm)
            E2=2*epsilon_0*c*e/h*np.sum(abs(F[:,:,:2])**2,axis=(1,2))*dx*dy
            flux=2*e/(h*mu_0)*np.sum(np.real(F[:,:,0]*np.conj(F[:,:,3])-F[:,:,1]*np.conj(F[:,:,2])),axis=1)*dx*dy
            E2[0]*=.5;flux[0]*=.5
            if nc%2==0:E2[-1]*=.5;flux[-1]*=.5
            direct=epsilon_0*c*np.sum(cube[:,:,:2]**2*win[:,None,None]**2)*dt*dx*dy/norm
            require(abs(E2.sum()*de/direct-1)<2e-11,'E-only Parseval failed')
            selected=abs(energy-photon)<=3.
            peak=selected&(energy>0)
            stats[label]=dict(E_only_band_J=float(E2[selected].sum()*de),
                native_signed_Poynting_band_J=float(flux[selected].sum()*de),
                band_centroid_eV=float(np.sum(energy[selected]*E2[selected])/np.sum(E2[selected])),
                peak_bin_eV=float(energy[np.flatnonzero(peak)[np.argmax(E2[peak])]]),
                spectral_bin_eV=de,window_energy_normalization=norm)
            if missing and label=='Hann':
                probes=set(np.linspace(5,len(raw)-6,16,dtype=int))
                for k in missing:
                    mid=int(np.searchsorted(ticks[lo:hi],k))
                    probes.update(range(max(4,mid-5),min(len(raw)-4,mid+5)))
                denominator=float(np.linalg.norm(F[selected,:,:2]))
                for j in sorted(probes):
                    k=ticks[lo+j]
                    if not target[0]<=k<=target[-1]:continue
                    candidates=np.delete(np.arange(len(raw)),j)
                    near=candidates[np.argsort(abs(ticks[lo:hi][candidates]-k))[:8]]
                    estimate=interpolate_one(ticks[lo:hi][near],raw[near],k)
                    error=estimate-raw[j]
                    # Empirical leave-one-out bound for one filled sample,
                    # including its actual window weight and the target band.
                    relative=np.sqrt(int(selected.sum()))*dt*win[int(k-target[0])]/np.sqrt(norm)*np.linalg.norm(error[:,:2])/denominator
                    cv.append(float(2*relative+relative*relative))
            cache[f'{name}_{label}_energy_eV']=energy
            cache[f'{name}_{label}_E_only_J_per_eV']=E2
            cache[f'{name}_{label}_native_Poynting_J_per_eV']=flux
        # Retain exact aperture/timestamps for independent GENESIS windowing.
        cache[name+'_time_s']=times;cache[name+'_x_m']=x;cache[name+'_y_m']=y;cache[name+'_roi']=roi
        summaries.append(dict(plane=name,z_m=z,n_samples=nc,roi_pixels=int(roi.sum()),
                              time_grid_jitter_steps=jitter,missing_native_steps=missing,
                              missing_sample_reconstruction='8-neighbor Lagrange; empirical held-out checks',
                              max_empirical_Hann_band_energy_error_bound=max(cv) if cv else 0.,
                              pixel_area_m2=dx*dy,retarded_first_s=float(times[0]),results=stats))
        print('Closed FEL plane analyzed',case['name'],name,flush=True)
    np.savez(out/'spectra.npz',**cache)
    save(out/'radiation.json',summaries)
    try:
        from report_discrete import load,calculate
        rows,sums,ranks,steps=load(d/'ledger',d/'ledger/energy-ledger.h5')
        index,q,budget=calculate(rows,sums,steps)
        save(out/'energy_ledger.json',budget)
    except Exception as exc:
        save(out/'LEDGER_ANALYSIS_FAILED.json',dict(error=str(exc),physical_acceptance=False))
        raise
    result=dict(case=case['name'],complete=True,physical_acceptance=False,
       raw_lab_total_fields=True,near_field_subtraction=False,
       primary_observable='Off-beam 90-130 um annulus, transverse E-only radiative energy proxy; native signed E x B separately reported',
       warning='E-only selection assumes a forward paraxial radiation band; it is not a separation of the complete raw field')
    save(out/'RESULT.json',result)
    return result
