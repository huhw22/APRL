"""Closed-output, absolute finite-aperture comparison. Never fit PIC amplitude."""
import json,csv,sys,time
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import c,e,h,hbar,m_e,mu_0
from theory import gaussian_point_reference,energy_density_spectrum
from prepare import require,save

def closed(path,group):
    with h5py.File(path,'r') as f:
        g=f[group];require(int(g['complete'][()])==1,'Incomplete '+str(path))
        a=g['records'][:];require(int(g['committed_records'][()])==len(a),'Uncommitted records')
    for k in a.dtype.names:
        if a[k].dtype.kind in 'fc':require(np.isfinite(a[k]).all(),'Nonfinite '+str(path))
    return a

def analyze(root):
    root=Path(root);plan=json.loads((root/'PLAN.json').read_text());p=plan['parameters'];d=Path(plan['case_directory'])
    out=root/'analysis';out.mkdir(exist_ok=False)
    from verify_reference import verify
    save(out/'REFERENCE_CHECKS.json',verify())
    with h5py.File(d/'detectors/screen.h5','r') as f:
        g=f['field_plane'];a=dict(g.attrs);require(int(g['complete'][()])==1,'Incomplete field plane')
        frame=a['frame'].decode() if isinstance(a['frame'],bytes) else a['frame'];require(frame=='laboratory','Wrong frame')
        n=int(g['committed_samples'][()]);require(n==len(g['fields']),'Uncommitted field tail')
        times=g['time_s'][:n];dt=float((times[-1]-times[0])/(n-1));require(n>32 and dt>0,'Insufficient time samples')
        jitter=float(np.max(abs(times-(times[0]+np.arange(n)*dt)))/dt);require(jitter<1e-4,'Nonuniform timestamps')
        nx=int(a['nx']);ny=int(a['ny']);dx=float(a['dx_m']);dy=float(a['dy_m'])
        x=float(a['x_first_m'])+np.arange(nx)*dx;y=float(a['y_first_m'])+np.arange(ny)*dy
        ix=np.flatnonzero(abs(x)<50e-6);iy=np.flatnonzero(abs(y)<50e-6)
        xs=x[ix];ys=y[iy];cube=np.empty((n,len(iy),len(ix),6))
        for start in range(0,n,16):
            z=min(n,start+16);records=g['fields'][start:z]
            E=records['electric_V_per_m'];B=records['magnetic_T']
            require(np.isfinite(E).all() and np.isfinite(B).all(),'Nonfinite full-plane field')
            cube[start:z,:,:,:3]=E[:,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
            cube[start:z,:,:,3:]=B[:,iy[0]:iy[-1]+1,ix[0]:ix[-1]+1]
    require(np.isfinite(times).all() and np.all(np.diff(times)>0),'Invalid time axis')
    energy=h*np.fft.rfftfreq(n,dt)/e;de=energy[1]-energy[0]
    magtime=-plan['reference_time_offset_s']+plan['magnet_entrance_z_m']/(p['beta']*c)
    first=float(times[0]-magtime-p['screen']/c);span=n*dt
    F=np.conj(np.fft.rfft(cube,axis=0))*dt*np.exp(1j*energy*e/hbar*first)[:,None,None,None]
    density=energy_density_spectrum(F[...,:3],F[...,3:]);density[0]*=.5
    if n%2==0:density[-1]*=.5
    time_energy=float(np.sum(cube[...,0]*cube[...,4]-cube[...,1]*cube[...,3])*dt*dx*dy/mu_0)
    fft_energy=float(density.sum()*de*dx*dy)
    parseval=abs(fft_energy-time_energy)/max(abs(time_energy),1e-300);require(parseval<1e-11,'Measured field Parseval failed')
    roi=(abs(ys[:,None])<10e-6)&(abs(xs[None,:])<10e-6)
    require(abs(float(roi.sum())*dx*dy-4e-10)<1e-20,
            'Registered 20 um square aperture area changed')
    band_limits=plan['registered_band_eV']
    band=(energy>=band_limits[0])&(energy<=band_limits[1]);ind=np.flatnonzero(band);require(len(ind)>=4,'Insufficient band bins')
    upper=200. if plan['sigma_z_m']==0 else 18.
    broader=(energy>=4)&(energy<=upper);frequencies=energy[broader]
    shapes=(len(frequencies),len(ys),len(xs),3)
    acc=np.zeros(shapes,complex);accB=np.zeros_like(acc);vel=np.zeros_like(acc);velB=np.zeros_like(acc)
    # Only the registered ROI requires absolute theory. A wider PIC image is
    # exported to expose aperture clipping; it is not a whole-space comparison.
    ntime=4*n
    for j,i in zip(*np.nonzero(roi)):
        print('Reference pixel',j,i,flush=True)
        vals=gaussian_point_reference(xs[i],ys[j],frequencies,p,
          sigma_xy=plan['sigma_x_m'],sigma_z=plan['sigma_z_m'],
          order=8 if plan['sigma_x_m'] else 1,time_first=first,time_span=span,ntime=ntime)
        for arr,v in zip((acc,accB,vel,velB),vals):arr[:,j,i]=v*plan['represented_electrons']
    rad=energy_density_spectrum(acc,accB);total=energy_density_spectrum(acc+vel,accB+velB)
    velonly=energy_density_spectrum(vel,velB);cross=total-rad-velonly
    measured=density[broader];bb=band[broader]
    aperture=lambda v:np.sum(v[:,roi],axis=1)*dx*dy
    sr,st,sv,sc,sm=map(aperture,(rad,total,velonly,cross,measured))
    integ=lambda v:float(np.sum(v[bb])*de)
    Wr,Wt,Wv,Wc,Wm=map(integ,(sr,st,sv,sc,sm))
    # Numerical reference uncertainty is measured separately from PIC error.
    j,i=tuple(np.argwhere(roi)[len(np.argwhere(roi))//2]);en=frequencies[bb]
    a12,b12,v12,vb12=gaussian_point_reference(xs[i],ys[j],en,p,
      sigma_xy=plan['sigma_x_m'],sigma_z=plan['sigma_z_m'],
      order=12 if plan['sigma_x_m'] else 1,time_first=first,time_span=span,ntime=8*n)
    base=total[bb,j,i]/plan['represented_electrons']**2
    conv=energy_density_spectrum(a12+v12,b12+vb12)
    conv_error=float(np.linalg.norm(conv-base)/np.linalg.norm(conv))
    np.savez(out/'field_spectra.npz',energy_eV=energy,raw_spectral_density_J_per_eV_m2=density,
      x_m=xs,y_m=ys,roi=roi,reference_energy_eV=frequencies,
      acceleration_E=acc,acceleration_B=accB,velocity_E=vel,velocity_B=velB,
      measured_Fourier_fields=F,time_s=times)
    with (out/'aperture_spectra.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['photon_eV','PIC_J_per_eV','acceleration_J_per_eV','total_LW_J_per_eV','velocity_J_per_eV','interference_J_per_eV'])
        writer.writerows(zip(frequencies,sm,sr,st,sv,sc))
    particles=[];entry=None
    for name in ['entry','bend-exit','exit']:
        z=closed(d/'detectors'/(name+'.h5'),'particle_plane');z=z[np.argsort(z['source_id'])]
        require(len(z)==plan['macroparticles'] and np.array_equal(z['source_id'],np.arange(1,len(z)+1)),'Missing particle crossings')
        require(np.all(z['crossing_direction']==1),'Backward crossings')
        w=abs(z['charge_C']);w=w/w.sum();u=z['proper_velocity'].astype(np.longdouble)
        gamma=np.sqrt(1+np.sum(u*u,axis=1))
        if entry is None:entry=gamma.copy()
        mean=lambda v:float(np.sum(w*v))
        particles.append({'plane':name,'count':len(z),'mean_angle_rad':mean(np.arctan2(u[:,0],u[:,2])),
          'mean_x_m':mean(z['position_m'][:,0]),'mean_lab_kinetic_change_eV':mean((gamma-entry)*m_e*c*c/e)})
    from report_discrete import load,calculate
    from relative_ledger import construct,write_tables,plot
    rows,sums,ranks,steps=load(d/'ledger',d/'ledger/energy-ledger.h5');index,q,budget=calculate(rows,sums,steps)
    absolute,pct,norm=construct(rows,sums,index,q);write_tables(out,index,absolute,pct,norm);plot(out,index,pct,norm,'Gaussian dipole: relative energy ledger')
    result={'status':'EXPLORATORY_ABSOLUTE_RADIATION_COMPARISON','physics_accepted':False,
      'raw_lab_total_fields':True,'window':'rectangular, DC retained, no background subtraction, native E/B time layers',
      'frames':n,'time_first_relative_to_magnet_light_s':first,'FFT_period_s':span,'frequency_spacing_eV':float(de),
      'registered_band_eV':band_limits,'included_bin_centers_eV':energy[band].tolist(),'registered_aperture_area_m2':float(roi.sum()*dx*dy),
      'PIC_band_J':Wm,'acceleration_reference_band_J':Wr,'complete_LW_reference_band_J':Wt,
      'velocity_only_band_J':Wv,'velocity_acceleration_interference_band_J':Wc,
      'PIC_relative_difference_from_acceleration':Wm/Wr-1,'PIC_relative_difference_from_complete_LW':Wm/Wt-1,
      'complete_LW_relative_difference_from_acceleration':Wt/Wr-1,
      'reference_GH8_time4N_vs_GH12_time8N_single_pixel_relative_L2':conv_error,
      'PIC_Parseval_relative_error':parseval,'timestamp_jitter_in_samples':jitter,
      'expected_exit_angle_rad':p['theta'],'particles':particles,'energy_ledger':budget,
      'scope':'Finite aperture and finite arc. Continuum coherent radiation only; no physical incoherent shot spectrum. No whole-space power or converged PIC acceptance claimed.'}
    save(out/'RESULT.json',result)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    axes[0].plot(frequencies,sm*1e18,'o-',label='PIC raw total');axes[0].plot(frequencies,st*1e18,label='Complete LW reference')
    axes[0].plot(frequencies,sr*1e18,'--',label='Acceleration reference');axes[0].axvspan(*band_limits,color='grey',alpha=.1)
    axes[0].set(xlabel='Photon energy (eV)',ylabel='Aperture energy (aJ/eV)',title='Registered central 20 x 20 um aperture');axes[0].legend()
    image=np.sum(density[band],axis=0)*de
    im=axes[1].imshow(image*1e12,origin='lower',extent=[(xs[0]-dx/2)*1e6,(xs[-1]+dx/2)*1e6,(ys[0]-dy/2)*1e6,(ys[-1]+dy/2)*1e6],aspect='equal')
    from matplotlib.patches import Rectangle
    axes[1].add_patch(Rectangle((-10,-10),20,20,fill=False,edgecolor='white'))
    axes[1].set(xlabel='x (um)',ylabel='y (um)',title=f'PIC signed {band_limits[0]:g}-{band_limits[1]:g} eV energy density');fig.colorbar(im,ax=axes[1],label='pJ / m2')
    fig.savefig(out/'radiation_comparison.png',dpi=180);fig.savefig(out/'radiation_comparison.pdf');plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)
