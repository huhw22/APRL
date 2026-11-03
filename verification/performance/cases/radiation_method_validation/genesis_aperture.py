"""Absolute GENESIS dump normalization and the same finite PIC aperture/window."""
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c,e,h
from scipy.interpolate import RectBivariateSpline
from prepare import require,save

def scalar(f,key):
    return float(np.asarray(f[key]).ravel()[0])

def compare(root,case,analysis_name='analysis-regularized-r1'):
    root=Path(root);plan=json.loads((root/'campaign.json').read_text())
    d=Path(case['directory']);out=d/analysis_name
    cache=np.load(out/'spectra.npz')
    pic=json.loads((out/'radiation.json').read_text())
    baseline=root/'boris'/analysis_name
    if case['name']!='boris' and (baseline/'GENESIS_APERTURE.json').exists():
        reference_cache=np.load(baseline/'spectra.npz')
        same=all(np.array_equal(cache[p['plane']+suffix],reference_cache[p['plane']+suffix])
                 for p in pic for suffix in ['_time_s','_x_m','_y_m','_roi'])
        if same:
            rows=json.loads((baseline/'GENESIS_APERTURE.json').read_text())
            measured={p['plane']:p for p in pic}
            for row in rows:
                value=measured[row['plane']]['results'][row['window']]['E_only_band_J']
                row.update(PIC_E_only_band_J=value,PIC_over_reference=value/row['band_J'],
                           reused_identical_geometry_time_reference=str(baseline))
            save(out/'GENESIS_APERTURE.json',rows)
            return rows
    entrance=np.load(root/'entrance_reference.npz')
    tau0=float(np.mean(entrance['tau_s']))
    gamma=plan['parameters']['gamma']
    ref=plan['geometry']['reference_time_offset_s']
    lam=plan['parameters']['resonant_wavelength_m']
    rows=[]
    for reference in plan['genesis_references']:
        if reference['name'] not in ['strong','weak','space-refined']:continue
        base=Path(reference['directory'])
        launch=json.loads((base/'genesis.log.launch.json').read_text())
        require(launch['exit_code']==0,'Reference process did not close')
        with h5py.File(base/'reference.out.h5','r') as f:
            zplot=f['Lattice/zplot'][:]
            standard=f['Field/power'][:]
            lattice_z=f['Lattice/z'][:];slip=f['Lattice/slippage'][:]
        for plane in pic:
            name=plane['plane']
            if name=='far':
                fieldfile=base/'far.fld.h5';steps=len(lattice_z)
                zsample=zplot[-1]
            else:
                periods=int(name[1:]);steps=periods*reference['substeps_per_period']
                fieldfile=base/f'reference.{steps}.fld.h5'
                zsample=periods*plan['parameters']['period_m']
            if not fieldfile.exists():continue
            times=cache[name+'_time_s']
            x=cache[name+'_x_m'];y=cache[name+'_y_m'];roi=cache[name+'_roi']
            xx,yy=np.meshgrid(x,y)
            with h5py.File(fieldfile,'r') as f:
                ng=int(scalar(f,'gridpoints'));ns=int(scalar(f,'slicecount'))
                dx=scalar(f,'gridsize');ds=scalar(f,'slicespacing')
                s0=scalar(f,'refposition')
                axis=(np.arange(ng)-(ng-1)/2)*dx
                optical=np.empty((ns,int(roi.sum())),complex)
                total=np.empty(ns)
                for i in range(ns):
                    g=f[f'slice{i+1:06d}']
                    a=g['field-real'][:]+1j*g['field-imag'][:]
                    a=a.reshape(ng,ng)
                    # GENESIS dump is sqrt(W) per grid cell. This is checked
                    # against its independent standard power output below.
                    total[i]=float(np.sum(abs(a)**2))
                    real=RectBivariateSpline(axis,axis,a.real,kx=3,ky=3)
                    imag=RectBivariateSpline(axis,axis,a.imag,kx=3,ky=3)
                    optical[i]=(real.ev(yy[roi],xx[roi])+1j*imag.ev(yy[roi],xx[roi]))/dx
            iz=int(np.argmin(abs(zplot-zsample)))
            # The step-refined diagnostic may omit the last free-drift sample;
            # main/weak/spatial dumps used here do sample the final endpoint.
            require(abs(zplot[iz]-zsample)<1e-9,'Power/dump times differ')
            mismatch=float(np.max(abs(total-standard[iz]))/max(np.max(total),1e-300))
            require(mismatch<1e-9,'GENESIS dump normalization differs from reported power')
            # Imported current is a cell average. Half a slice aligns its
            # charge centroid with the known PIC crossing-event input.
            tgen=(s0+(np.arange(ns)+.5)*ds+np.sum(slip[:steps])*lam)/c+tau0-ref
            normscale=(plan['parameters']['peak_current_A']/reference['current_A'])**2
            pdt=float(np.median(np.diff(times)))
            area=plane['pixel_area_m2']
            for label in ['rectangular','Hann']:
                if label=='Hann':
                    T=times[-1]-times[0]
                    win=.5-.5*np.cos(2*np.pi*(tgen-times[0])/T)
                    win[(tgen<times[0])|(tgen>times[-1])]=0.
                    norm=np.mean(np.hanning(len(times))**2)
                else:
                    win=((tgen>=times[0])&(tgen<=times[-1])).astype(float);norm=1.
                # Envelope has both baseband signs. A symmetric +/-3 eV band
                # is invariant under carrier-complex-conjugation convention.
                F=np.fft.fft(optical*win[:,None],axis=0)*(ds/c)/np.sqrt(norm)
                baseband=h*np.fft.fftfreq(ns,ds/c)/e
                selected=abs(baseband)<=3.
                integral=float(np.sum(abs(F[selected])**2)*area/(ns*ds/c))*normscale
                full=float(np.sum(abs(optical)**2*win[:,None]**2)*area*ds/c/norm)*normscale
                rows.append(dict(reference=reference['name'],plane=name,z_m=plane['z_m'],window=label,
                  current_squared_rescale=normscale,band_J=integral,all_envelope_J=full,
                  dump_vs_standard_power_max_relative=mismatch,
                  PIC_E_only_band_J=plane['results'][label]['E_only_band_J'],
                  PIC_over_reference=plane['results'][label]['E_only_band_J']/integral,
                  reference_slice_s=ds/c,carrier_phase_convention_not_used_for_symmetric_band=True))
            print('GENESIS aperture',case['name'],reference['name'],name,flush=True)
    save(out/'GENESIS_APERTURE.json',rows)
    return rows
