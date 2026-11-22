"""Validate closed native outputs and summarize the finite virtual-screen aperture."""
from pathlib import Path
import csv
import hashlib
import json
import sys

import h5py
import numpy as np
from scipy.constants import c, mu_0

CASE=Path(sys.argv[1]).resolve()
ROOT=CASE.parent


def main():
    assert (CASE/'EXIT_CODE').read_text().strip()=='0'
    text=(CASE/'aprl.log').read_text()
    assert 'Configured stop reached:' in text
    crossings={};bunching={}
    frequencies=np.array([8.19e12,4.095e12,2.73e12])
    for name in ('input','undulator-entry','middle','exit'):
        with h5py.File(CASE/'detectors'/f'{name}.h5') as f:
            g=f['particle_plane'];assert int(g['complete'][()])==1
            a=g['records'][:]
            assert len(a)==int(g['committed_records'][()])==80000
            assert len(np.unique(a['source_id']))==80000
            for field in ('time_s','position_m','proper_velocity','charge_C'):
                assert np.isfinite(a[field]).all()
            charge=float(np.sum(a['charge_C'],dtype=np.longdouble))
            assert abs(charge/.7e-9+1)<1e-12
            tt=a['time_s'].astype(np.longdouble);tt=tt-np.mean(tt)
            b=[float(abs(np.mean(np.exp(2j*np.pi*f*np.asarray(tt,dtype=float))))) for f in frequencies]
            crossings[name]=dict(records=len(a),charge_C=charge,
                 rms_xy_m=np.std(a['position_m'][:,:2],axis=0).tolist(),
                 xy_min_m=np.min(a['position_m'][:,:2],axis=0).tolist(),
                 xy_max_m=np.max(a['position_m'][:,:2],axis=0).tolist())
            bunching[name]=b
    with h5py.File(CASE/'radiation.h5') as f:
        assert int(f.attrs['complete'])==1
        e=f['E_time'][:];b=f['B_time'][:]
        assert e.shape==b.shape and e.shape[1:]==(19,25,3)
        assert np.isfinite(e).all() and np.isfinite(b).all() and np.max(abs(e))>0
        dt=float(f.attrs['time_step_s'])
        signed_fluence=np.sum(np.cross(e,b)[...,2],axis=0)*dt/mu_0
        electric_exposure=np.sum(e*e,axis=(0,3))*dt
        edge_norm=float(np.sum(e[:8]**2)+np.sum(e[-8:]**2))
        full_norm=float(np.sum(e**2))
        field_window=dict(time_bins=len(e),time_step_s=dt,
              time_bounds_s=[float(f.attrs['time_start_s']),float(f.attrs['time_start_s'])+len(e)*dt],
              first_last_8_bins_electric_square_fraction=edge_norm/full_norm,
              accelerated_segments=int(f.attrs['accelerated_segments']),
              committed_substeps=int(f.attrs['committed_substeps']),
              finite_aperture_only=True)
        np.savez_compressed(CASE/'screen-summary.npz',x_m=f['x_m'][:],y_m=f['y_m'][:],
             signed_acceleration_fluence_J_m2=signed_fluence,electric_exposure_V2s_m2=electric_exposure,
             time_s=f['time_s'][:],central_E_V_m=e[:,9,12,:],central_B_T=b[:,9,12,:])
    with h5py.File(CASE/'spectrum.h5') as f:
        assert int(f.attrs['complete'])==1
        freq=f['frequency_hz'][:]
        E=f['E_spectrum'][...];E=E[...,0]+1j*E[...,1]
        B=f['B_spectrum'][...];B=B[...,0]+1j*B[...,1]
        assert np.isfinite(E).all() and np.isfinite(B).all()
        density=2*np.cross(E,B.conj())[...,2].real/mu_0
        # Finite-aperture flux from the acceleration-only fields, not total
        # radiation energy, near-field energy or a full solid-angle integral.
        x=f['x_m'][:];y=f['y_m'][:]
        flux=np.trapz(np.trapz(density,x=x,axis=2),x=y,axis=1)
        central=np.sum(abs(E[:,9,12,:])**2,axis=1)
        with (CASE/'spectrum.csv').open('w') as out:
            writer=csv.writer(out);writer.writerow(['frequency_THz','central_E_spectral_norm2','aperture_signed_flux_J_Hz'])
            writer.writerows(zip(freq/1e12,central,flux))
        np.savez_compressed(CASE/'spectral-maps.npz',frequency_Hz=freq,x_m=x,y_m=y,
                electric_spectral_norm2=np.sum(abs(E)**2,axis=3),signed_flux_density_J_Hz_m2=density)
    result=dict(complete=True,parameters=json.loads((CASE/'PARAMETERS.json').read_text()),
         crossings=crossings,bunching_frequencies_Hz=frequencies.tolist(),bunching=bunching,
         screen=field_window,source_revision=json.loads((ROOT/'PLAN.json').read_text())['source_revision'])
    (CASE/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    plot()
    print(json.dumps(result,ensure_ascii=False),flush=True)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=np.load(CASE/'screen-summary.npz');s=np.load(CASE/'spectral-maps.npz')
    par=json.loads((CASE/'PARAMETERS.json').read_text());target=par['target_frequency_Hz']
    f=s['frequency_Hz'];k=int(np.argmin(abs(f-target)))
    fig,axes=plt.subplots(1,3,figsize=(15,4.2),constrained_layout=True)
    central=s['electric_spectral_norm2'][:,9,12]
    axes[0].plot(f/1e12,central);axes[0].axvline(target/1e12,color='gray',ls='--')
    axes[0].set(xlabel='Frequency [THz]',ylabel='Central spectral |E|^2 [(V s/m)^2]',xlim=(0,40))
    im=axes[1].pcolormesh(s['x_m']*1e3,s['y_m']*1e3,s['electric_spectral_norm2'][k],shading='nearest')
    axes[1].set(xlabel='x [mm]',ylabel='y [mm]',title=f'{f[k]/1e12:g} THz spectral |E|^2')
    fig.colorbar(im,ax=axes[1])
    im=axes[2].pcolormesh(r['x_m']*1e3,r['y_m']*1e3,r['electric_exposure_V2s_m2'],shading='nearest')
    axes[2].set(xlabel='x [mm]',ylabel='y [mm]',title='Time-integrated |E|^2 [V^2 s/m^2]')
    fig.colorbar(im,ax=axes[2])
    fig.suptitle(f'APRL 2.1: 20 x 5 cm, K={par["K"]:.5f}; acceleration-only finite aperture')
    for ext in ('png','pdf'):fig.savefig(CASE/f'radiation-summary.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
