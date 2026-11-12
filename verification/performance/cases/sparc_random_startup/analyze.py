"""Raw first-segment metrics, with a fixed central input-ID cohort."""
import csv
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c,mu_0

from common import load, require, save


def table(path,rows):
    with Path(path).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def moments(ids,weights,phase,prefix,first,last):
    group=np.searchsorted(prefix,ids-1,side='right')-1
    use=(group>=first)&(group<last)
    group=group[use]-first;weights=weights[use];phase=phase[use]
    q=np.bincount(group,weights=weights,minlength=last-first)
    b=np.bincount(group,weights=weights*phase.real,minlength=last-first)+1j*np.bincount(group,weights=weights*phase.imag,minlength=last-first)
    active=q>0;total=q.sum()
    rms2=float(np.sum(abs(b[active])**2/q[active])/total)
    glob=complex(b.sum()/total)
    return dict(records=int(use.sum()),represented_electrons=float(total),
        global_b1=abs(glob),local_rms_b1=np.sqrt(rms2),
        global_b1_real=glob.real,global_b1_imag=glob.imag,
        common_phase_bf_fraction=abs(glob)**2/max(rms2,1e-300),
        noise_b_rms=float(np.sqrt(np.sum(weights**2))/total))


def pic_spectrum(path,wavelength):
    with h5py.File(path,'r') as f:
        g=f['field_plane'];require(int(g['complete'][()])==1,'Incomplete field')
        t=g['time_s'][:];raw=g['fields'];dt=float(np.median(np.diff(t)))
        ticks=np.rint((t-t[0])/dt).astype(int)
        require(np.all(np.diff(ticks)>0),'Repeated time samples')
        require(np.max(abs((t-t[0])/dt-ticks))<1e-4,'Irregular field sampling')
        N=int(ticks[-1]+1);missing=np.setdiff1d(np.arange(N),ticks)
        require(len(missing)<=1,'Too many missing field samples')
        # Original40x40, six CPML cells per side: exact ±560um aperture.
        area=float(g.attrs['dx_m'])*float(g.attrs['dy_m'])
        values=np.empty((N,28*28,6))
        for start in range(0,len(t),32):
            block=raw[start:start+32];stop=start+len(block)
            values[ticks[start:stop],:,:3]=block['electric_V_per_m'][:,6:34,6:34].reshape(len(block),-1,3)
            values[ticks[start:stop],:,3:]=block['magnetic_T'][:,6:34,6:34].reshape(len(block),-1,3)
        for k in missing:
            nodes=ticks[np.argsort(abs(ticks-k))[:8]]
            coefficients=np.array([np.prod([(k-j)/(i-j) for j in nodes if j!=i]) for i in nodes])
            values[k]=np.tensordot(coefficients,values[nodes],axes=1)
    frequency=np.fft.rfftfreq(N,dt);band=abs(frequency/(c/wavelength)-1)<=.05
    factor=np.full(len(frequency),2.);factor[0]=1.
    if N%2==0:factor[-1]=1.
    result=dict(samples=N,missing_ticks=missing.tolist(),duration_s=N*dt)
    for label,window in [('rect',np.ones(N)),('hann',np.hanning(N))]:
        window=window/np.sqrt(np.mean(window**2))
        spectrum=np.zeros(len(frequency));direct=0.
        for start in range(0,784,32):
            v=values[:,start:start+32,:]*window[:,None,None]
            direct+=float(np.sum(v[:,:,0]*v[:,:,4]-v[:,:,1]*v[:,:,3]))*dt*area/mu_0
            z=np.fft.rfft(v,axis=0)*dt
            spectrum+=(z[:,:,0]*z[:,:,4].conj()-z[:,:,1]*z[:,:,3].conj()).real.sum(axis=1)
        spectrum*=factor*area/(mu_0*N*dt)
        error=abs(spectrum.sum()-direct)/max(abs(spectrum).sum(),abs(direct),1e-300)
        require(error<2e-10,'Poynting Parseval failed')
        result[label+'_band_J']=float(spectrum[band].sum())
        result[label+'_full_frequency_J']=float(spectrum.sum())
        result[label+'_parseval_relative']=error
    return result


def genesis_spectrum(path,wavelength):
    with h5py.File(path,'r') as f:
        n=int(f['slicecount'][0]);ng=int(f['gridpoints'][0]);dx=float(f['gridsize'][0])
        dt=float(f['slicespacing'][0])/c
        x=(np.arange(ng)-(ng-1)/2)*dx
        weights=np.maximum(0,np.minimum(x+dx/2,.00056)-np.maximum(x-dx/2,-.00056))/dx
        ids=np.flatnonzero(weights>0);area=weights[ids,None]*weights[None,ids]
        values=np.empty((n,len(ids),len(ids)),complex)
        for j in range(n):
            group=f[f'slice{j+1:06d}']
            field=group['field-real'][:]+1j*group['field-imag'][:]
            values[j]=field.reshape(ng,ng)[np.ix_(ids,ids)]
    transform=np.fft.fft(values,axis=0)*dt
    band=abs(np.fft.fftfreq(n,dt))<=.05*c/wavelength
    native=float(np.sum(abs(values)**2*area)*dt)
    spectral=np.sum(abs(transform)**2*area,axis=(1,2))/(n*dt)
    require(abs(spectral.sum()-native)<2e-11*max(native,1e-300),'GENESIS Parseval failed')
    return dict(rect_band_J=float(spectral[band].sum()),rect_all_envelope_J=native)


def analyze(root,entry):
    root=Path(root);plan=load(root/'PLAN.json');d=Path(entry['directory'])
    output=d/'analysis';output.mkdir(exist_ok=False)
    info=load(d/'INPUT_VALIDATED.json');cohort=load(d/'COHORT.json')
    prefix=np.array(cohort['bin_prefix'],dtype=np.uint64);rows=[];spectra=[]
    for name,z in plan['particle_planes']:
        with h5py.File(d/'detectors'/(name+'-particles.h5'),'r') as f:
            g=f['particle_plane'];require(int(g['complete'][()])==1,'Particle plane incomplete')
            a=g['records'][:]
        require(len(np.unique(a['source_id']))==len(a),'Duplicate particle crossings')
        phase=np.exp(2j*np.pi*np.longdouble(c)*(a['time_s'].astype(np.longdouble)-np.longdouble(z)/c)/plan['wavelength_m']).astype(complex)
        for label,first,last in [('all',0,231),('core',77,154)]:
            result=moments(a['source_id'],a['weight'],phase,prefix,first,last)
            rows.append(dict(program='PIC',scope=label,plane=name,z_m=z,**result))
        if name.startswith('u') and name!='und-entry':
            spectra.append(dict(program='PIC',plane=name,z_m=z,**pic_spectrum(d/'detectors'/(name+'.h5'),plan['wavelength_m'])))
    table(output/'pic_bunching.csv',rows);save(output/'pic_spectra.json',spectra)
    genrows=[];genenergy=[]
    for mode in ('off','longrange'):
        gd=d/('genesis-'+mode);mapping=load(d/'NATIVE_MAPPING.json')
        with h5py.File(gd/'reference.out.h5','r') as f:
            z=f['Lattice/zplot'][:];current=f['Beam/current'][:]
            bf=f['Beam/bunching'][:]*np.exp(1j*f['Beam/bunchingphase'][:])
            power=f['Field/Global/energy'][:];slices=np.arange(current.shape[1])+mapping['integer_origin']
            # Eulerian native slices retained here; these are not fixed IDs.
            for k,zk in enumerate(z):
                w=current[k];total=w.sum();global_b=complex((w*bf[k]).sum()/total)
                rms2=float((w*abs(bf[k])**2).sum()/total)
                genrows.append(dict(program='GENESIS-'+mode,z_m=float(zk),
                    global_b1=abs(global_b),local_rms_b1=np.sqrt(rms2),
                    common_phase_bf_fraction=abs(global_b)**2/max(rms2,1e-300),
                    full_envelope_energy_J=float(power[k])))
        dumps=sorted(gd.glob('reference.*.fld.h5'),key=lambda p:int(p.name.split('.')[1]))
        if len(dumps)==4: dumps.append(gd/'exit.fld.h5')
        require(len(dumps)==5,'Missing intermediate GENESIS field dumps')
        for dump,(name,zk) in zip(dumps,plan['particle_planes'][2:]):
            genenergy.append(dict(program='GENESIS-'+mode,plane=name,z_m=zk,
                                 **genesis_spectrum(dump,plan['wavelength_m'])))
    table(output/'genesis_bunching.csv',genrows);table(output/'genesis_band_energy.csv',genenergy)
    save(output/'RESULT.json',dict(complete=True,physical_acceptance=False,
        actual_weight=info['actual_weight'],physical_seed=entry['physical_seed'],
        scope='PIC fixed original-bin/ID cohort; GENESIS Eulerian optical slices, all source. Raw band energies use matched ±560um aperture; no near-field subtraction or amplitude fitting.',
        cautions=['Weighted input has enhanced shot noise; no one-electron SASE acceptance.',
                  'Finite hard crop, one segment and only two seeds; no saturation claim.',
                  'GENESIS slice slippage labels differ from fixed source IDs; do not equate core/local curves without subsequent matched-cohort review.']))
