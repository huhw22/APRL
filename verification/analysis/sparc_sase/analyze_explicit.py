"""Closed SPARC outputs only; raw laboratory fields are not separated radiation."""
import os,sys,json,csv,gc,math,time
from pathlib import Path
import h5py,numpy as np,yaml
from scipy import fft,ndimage
from scipy.constants import c,h,e,m_e,mu_0,epsilon_0

def need(x,msg):
    if not x:raise RuntimeError(msg)
def save(p,x):
    Path(p).write_text(json.dumps(x,indent=2,allow_nan=False,
      default=lambda q:q.item() if isinstance(q,np.generic) else q.tolist())+'\n')
def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def lagrange(nodes,target):
    x=np.array(nodes,dtype=float)-target
    return np.array([np.prod(-np.delete(x,j)/(x[j]-np.delete(x,j))) for j in range(len(x))])

def particle_read(path):
    with h5py.File(path,'r') as f:
        g=f['particle_plane'];n=int(g['committed_records'][()])
        need(int(g['complete'][()])==1 and n==len(g['records']),'Unclosed particles')
        a=g['records'][:]
    need(np.all(a['crossing_direction']==1),'Upstream crossings require signed treatment')
    for k in a.dtype.names:
        if a[k].dtype.kind=='f':need(np.isfinite(a[k]).all(),'Nonfinite particles')
    a=a[np.argsort(a['source_id'])]
    need(np.all(np.diff(a['source_id'])>0),'Repeated particle IDs')
    return a

def particles(root,out,plan):
    directory=root/'boris/detectors';entry=particle_read(directory/'entry-particles.h5')
    n=plan['input_particles'];need(np.array_equal(entry['source_id'],np.arange(1,n+1)),'Input ID coverage')
    w=abs(entry['charge_C'])/e
    need(np.max(abs(w-entry['weight'])/w)<1e-12,'Weight convention')
    need(np.max(abs(entry['mass_kg']/m_e/w-1))<1e-12,'Macro mass convention')
    planes=yaml.safe_load((root/'boris/run.yaml').read_text())['detectors']['particle_planes']
    common=np.ones(n,bool);coverage=[]
    for p in planes:
        a=particle_read(directory/(p['name']+'.h5'));idx=a['source_id'].astype(int)-1
        seen=np.zeros(n,bool);seen[idx]=True;common&=seen
        missed=~seen
        coverage.append(dict(plane=p['name'],z_m=p['z'],count=len(a),
          charge_fraction=float(w[seen].sum()/w.sum()),
          missing_input_time_ps_min=float(entry['time_s'][missed].min()*1e12) if missed.any() else None,
          missing_input_time_ps_max=float(entry['time_s'][missed].max()*1e12) if missed.any() else None))
    u0=entry['proper_velocity'];g0=np.sqrt(1+(u0*u0).sum(1));omega=2*np.pi*c/plan['physics']['wavelength_m']
    results=[]
    for p,cov in zip(planes,coverage):
        a=particle_read(directory/(p['name']+'.h5'));idx=a['source_id'].astype(int)-1;keep=common[idx]
        a=a[keep];idx=idx[keep];weights=w[idx];weights/=weights.sum()
        u=a['proper_velocity'];gamma=np.sqrt(1+(u*u).sum(1))
        dg=((u-u0[idx])*(u+u0[idx])).sum(1)/(gamma+g0[idx])
        mean=float(weights@dg)
        b=weights@np.exp(1j*omega*(a['time_s']-p['z']/c+plan['geometry']['reference_time_offset_s']))
        x=a['position_m'][:,0];y=a['position_m'][:,1]
        results.append(dict(**cov,common_count=len(a),common_charge_fraction=float(w[common].sum()/w.sum()),
          mean_delta_kinetic_eV=mean*m_e*c*c/e,
          std_delta_kinetic_eV=float(np.sqrt(weights@((dg-mean)**2)))*m_e*c*c/e,
          total_delta_kinetic_J=float(np.sum(w[common]))*mean*m_e*c*c,
          global_b1=float(abs(b)),sigma_x_um=float(np.sqrt(weights@((x-weights@x)**2)))*1e6,
          sigma_y_um=float(np.sqrt(weights@((y-weights@y)**2)))*1e6))
    save(out/'particles.json',results);table(out/'particles.csv',results)
    # Independent transport ledger uses endpoint four-momenta, not gamma of
    # interpolated momenta. Compare the same covered charges before differencing.
    transport=[]
    for name in ['entry-transport','und-entry-transport','exit-transport']:
        with h5py.File(directory/(name+'-energy.h5'),'r') as f:
            g=f['energy_plane'];need(int(g['complete'][()])==1,'Unclosed transport')
            a=g['records'][-1];count=int(a['particle_downstream_count'])
            K=np.longdouble(a['particle_kinetic_lab_J'])+np.longdouble(a['particle_kinetic_lab_J_low'])
            q=float(abs(a['particle_charge_C']))
            other={'entry-transport':'entry-particles','und-entry-transport':'und-entry-particles',
                   'exit-transport':'far-particles'}[name]
            at=particle_read(directory/(other+'.h5'));idx=at['source_id'].astype(int)-1
            need(count==len(idx),'Transport and raw crossing counts differ')
            initial=float(np.sum(w[idx]*(g0[idx]-1))/np.sum(w[idx]))*m_e*c*c/e
            transport.append(dict(plane=name,count=count,charge_C=q,mean_kinetic_eV=float(K/(q/e)/e),
              mean_delta_kinetic_eV=float(K/(q/e)/e)-initial,
              physical_field_energy_J=float(a['physical_field_energy_lab_J']),
              field_sample_count=int(a['field_sample_count']),field_gap_count=int(a['field_gap_count'])))
    save(out/'transport.json',transport)
    return results

def pic_plane(root,out,plan,name,case='boris'):
    started=time.monotonic();path=root/case/'detectors'/(name+'.h5')
    dt=plan['geometry']['dt_s']/plan['physics']['boost_gamma'];dx=plan['geometry']['dx_m'];area=dx*dx
    with h5py.File(path,'r') as f:
        g=f['field_plane'];n=int(g['committed_samples'][()])
        need(int(g['complete'][()])==1 and n==len(g['fields'])==len(g['time_s']),'Incomplete field')
        t=g['time_s'][:];z=float(g.attrs['plane_z_m'])
        ticks=np.rint((t-t[0])/dt).astype(int);need(np.all(np.diff(ticks)>0),'Repeated time')
        need(np.max(abs((t-t[0])/dt-ticks))<1e-4,'Off-grid timestamps')
        N=int(ticks[-1]+1);missing=np.setdiff1d(np.arange(N),ticks)
        need(len(missing)<=1,'More than one missing sample')
        ny,nx=g['fields'].shape[1:];cx=cy=10;px=nx-2*cx;py=ny-2*cy
        cube=np.empty((N,py*px,6))
        for a in range(0,n,32):
            b=min(n,a+32);v=g['fields'][a:b]
            need(np.isfinite(v['electric_V_per_m']).all() and np.isfinite(v['magnetic_T']).all(),'Nonfinite field')
            cube[ticks[a:b],:,:3]=v['electric_V_per_m'][:,cy:ny-cy,cx:nx-cx,:].reshape(b-a,-1,3)
            cube[ticks[a:b],:,3:]=v['magnetic_T'][:,cy:ny-cy,cx:nx-cx,:].reshape(b-a,-1,3)
        for k in missing:
            neighbors=ticks[np.argsort(abs(ticks-k))[:8]]
            cube[k]=np.tensordot(lagrange(neighbors,k),cube[neighbors],axes=(0,0))
        retarded=t[0]+np.arange(N)*dt-z/c
        x=float(g.attrs['x_first_m'])+np.arange(cx,nx-cx)*dx
        y=float(g.attrs['y_first_m'])+np.arange(cy,ny-cy)*dx
    energy=h/e*fft.rfftfreq(N,dt);de=energy[1];one=np.full(len(energy),2.);one[0]=1.
    if N%2==0:one[-1]=1.
    photon=h*c/e/plan['physics']['wavelength_m'];band=abs(energy/photon-1)<=.05
    arr=dict(energy_eV=energy,retarded_time_s=retarded,x_m=x,y_m=y);summ={}
    # Interpolation is diagnosed by held-out observed points around the gap.
    # These values do not overwrite the native files or the observed samples.
    holdout=[]
    for k in missing:
        for j in [int(k-3),int(k-1),int(k+1),int(k+3)]:
            nodes=ticks[np.argsort(abs(ticks-j))[1:9]]
            estimate=np.tensordot(lagrange(nodes,j),cube[nodes],axes=(0,0))
            holdout.append(float(np.linalg.norm(estimate[:,:2]-cube[j,:,:2])))
    for label,window in [('rect',np.ones(N)),('hann',np.hanning(N))]:
        norm=float(np.mean(window*window));window=window/np.sqrt(norm)
        E2=np.zeros(len(energy));flux=np.zeros(len(energy));Ez2=np.zeros(len(energy))
        directE=0.;directS=0.;power=np.zeros(N)
        for a in range(0,px*py,64):
            v=np.array(cube[:,a:a+64,:].transpose(1,2,0),order='C',copy=True)
            v*=window[None,None,:]
            directE+=float(np.sum(v[:,:2,:]**2))*dt*area*epsilon_0*c
            s=(v[:,0,:]*v[:,4,:]-v[:,1,:]*v[:,3,:]).sum(0)*area/mu_0
            directS+=float(s.sum())*dt;power+=s
            F=fft.rfft(v,axis=-1,workers=8,overwrite_x=True)*dt
            E2+=np.sum(abs(F[:,:2,:])**2,axis=(0,1))
            Ez2+=np.sum(abs(F[:,2,:])**2,axis=0)
            flux+=np.sum((F[:,0,:]*F[:,4,:].conj()-F[:,1,:]*F[:,3,:].conj()).real,axis=0)
        E2*=one*area*epsilon_0*c/(N*dt*de);Ez2*=one*area*epsilon_0*c/(N*dt*de)
        flux*=one*area/(mu_0*N*dt*de)
        errE=abs(E2.sum()*de-directE)/max(abs(directE),1e-300)
        errS=abs(flux.sum()*de-directS)/max(np.sum(abs(flux))*de,abs(directS),1e-300)
        need(errE<2e-10 and errS<2e-10,'Parseval failure')
        eb=float(E2[band].sum()*de);sb=float(flux[band].sum()*de)
        peak=int(np.flatnonzero(band)[np.argmax(E2[band])])
        # Cauchy/Parseval gives a conservative perturbation bound from the
        # largest observed leave-one-out error at one reconstructed point.
        dnorm=max(holdout,default=0)*max(window[missing],default=0)
        eps=epsilon_0*c*area*dt*dnorm*dnorm
        bound=2*np.sqrt(max(eps,0)/max(eb,1e-300))+eps/max(eb,1e-300)
        summ[label]=dict(E_transverse_band_proxy_J=eb,native_Poynting_band_J=sb,
          raw_full_frequency_Poynting_J=float(flux.sum()*de),
          peak_energy_eV=float(energy[peak]),
          centroid_eV=float(np.sum(energy[band]*E2[band])/np.sum(E2[band])),
          parseval_E_relative=errE,parseval_flux_relative=errS,
          window_mean_square=norm,empirical_one_missing_sample_energy_bound=float(bound))
        arr[label+'_E_proxy_J_per_eV']=E2;arr[label+'_native_flux_J_per_eV']=flux
        arr[label+'_Ez_proxy_J_per_eV']=Ez2
        if label=='rect':arr['native_total_power_W']=power
    np.savez_compressed(out/(name+'-spectra.npz'),**arr)
    result=dict(plane=name,z_m=z,raw_samples=n,fft_samples=N,dt_s=dt,bin_eV=de,
      missing_ticks=missing.tolist(),missing_retarded_ps=(retarded[missing]*1e12).tolist(),
      missing_sample_method='8-neighbor Lagrange, plus held-out observed-sample diagnostic',
      first_retarded_s=float(retarded[0]),last_retarded_s=float(retarded[-1]),
      physical_aperture_half_width_m=plan['geometry']['physical_half_width_m'],
      results=summ,seconds=time.monotonic()-started)
    save(out/(name+'-summary.json'),result);print('PIC',name,summ['rect'],flush=True)
    del cube;gc.collect();return result

def genesis(root,out,plan,pic,reference_directory=None):
    base=root/'genesis' if reference_directory is None else Path(reference_directory)
    photon=h*c/e/plan['physics']['wavelength_m']
    with h5py.File(base/'reference.out.h5','r') as f:
        z=f['Lattice/z'][:];dz=f['Lattice/dz'][:];aw=f['Lattice/aw'][:]
        gamma0=plan['physics']['gamma'];uz=np.sqrt(gamma0*gamma0-1-aw*aw)
        reference_lag=dz*(1+aw*aw)/(uz*(gamma0+uz))
        dumpsteps=f['Meta/Fielddumps/intstep'][:].astype(int)
        zplot=f['Lattice/zplot'][:];I=f['Beam/current'][:];ga=f['Beam/energy'][:]
        b=f['Beam/bunching'][:]*np.exp(1j*f['Beam/bunchingphase'][:])
        eg=(I*ga).sum(1)/I.sum(1);bg=abs((I*b).sum(1)/I.sum(1))
        np.savez(out/'genesis-beam.npz',z_m=zplot,mean_delta_eV=(eg-eg[0])*m_e*c*c/e,
          global_b1=bg,local_b1_rms=np.sqrt((I*abs(b)**2).sum(1)/I.sum(1)))
    rows=[]
    for p in pic:
        name=p['plane'];position=p['z_m']
        if name=='far':path=base/'far.fld.h5';step=len(z)
        else:
            step=int(dumpsteps[np.argmin(abs(z[dumpsteps]-position))])
            need(abs(z[step]-position)<1e-8,'Wrong dump plane');path=base/f'reference.{step}.fld.h5'
        cache=np.load(out/(name+'-spectra.npz'));x=cache['x_m'];y=cache['y_m'];area=(x[1]-x[0])*(y[1]-y[0])
        with h5py.File(path,'r') as f:
            ns=int(f['slicecount'][0]);ng=int(f['gridpoints'][0]);ds=float(f['slicespacing'][0])
            dx=float(f['gridsize'][0]);s0=float(f['refposition'][0]);dt=ds/c
            axis=(np.arange(ng)-(ng-1)/2)*dx
            xx,yy=np.meshgrid((x-axis[0])/dx,(y-axis[0])/dx);coords=np.array([yy,xx])
            optical=np.empty((ns,len(x)*len(y)),complex);power=np.empty(ns);imgsum=np.zeros((ng,ng))
            for j in range(ns):
                g=f[f'slice{j+1:06d}'];a=g['field-real'][:].reshape(ng,ng)+1j*g['field-imag'][:].reshape(ng,ng)
                need(np.isfinite(a).all(),'Nonfinite reference')
                power[j]=np.sum(abs(a)**2);imgsum+=abs(a)**2
                optical[j]=(ndimage.map_coordinates(a.real,coords,order=3,mode='constant')+
                  1j*ndimage.map_coordinates(a.imag,coords,order=3,mode='constant')).ravel()/dx
            # Independent dump intensity projection is saved before conversion.
            # Its relation to sqrt(W)/cell is explicit in writeFieldHDF5.cpp.
            # Exact constants of installed GENESIS GenMain.cpp, not rounded mc2.
            scale=dx*510999.06/(2*np.pi/plan['physics']['wavelength_m'])/np.sqrt(376.73)
            projection=f['int_xy'][:].reshape(ng,ng)*scale*scale
            normalization_error=float(np.max(abs(imgsum-projection))/max(imgsum.max(),1e-300))
            need(normalization_error<1e-10,'GENESIS field normalization disagreement')
        # G4 slices increase toward the head; positive slippage advances the
        # field to larger index. Convert to increasing physical arrival time.
        tgen=(plan['physics']['bunch_length_m']-s0-(np.arange(ns)+.5)*ds+
              np.sum(reference_lag[:step]))/c-plan['geometry']['reference_time_offset_s']
        tgen=tgen[::-1];optical=optical[::-1];power=power[::-1]
        freqs=fft.fftfreq(ns,dt);energy=photon-h/e*freqs
        order=np.argsort(energy);energy=energy[order];de=h/e/(ns*dt);band=abs(energy/photon-1)<=.05
        full=float(power.sum()*dt);same=float(np.sum(abs(optical)**2)*dt*area)
        arr=dict(energy_eV=energy,retarded_time_s=tgen,power_native_W=power)
        stats={}
        for label in ['rect','hann']:
            if label=='rect':win=np.ones(ns);norm=1.
            else:
                win=.5-.5*np.cos(2*np.pi*(tgen-p['first_retarded_s'])/((p['fft_samples']-1)*p['dt_s']))
                win[(tgen<p['first_retarded_s'])|(tgen>p['last_retarded_s'])]=0.
                norm=p['results']['hann']['window_mean_square']
            # The rectangular primary reference is also clipped to the actual
            # physical observation interval; the full reference stays separate.
            win[(tgen<p['first_retarded_s'])|(tgen>p['last_retarded_s'])]=0.
            F=fft.fft(optical*win[:,None]/np.sqrt(norm),axis=0,workers=8)*dt
            density=np.sum(abs(F)**2,axis=1)[order]*area/(ns*dt*de)
            direct=float(np.sum(abs(optical)**2*win[:,None]**2)*area*dt/norm)
            error=abs(density.sum()*de-direct)/max(direct,1e-300);need(error<2e-10,'Reference Parseval')
            stats[label]=dict(band_J=float(density[band].sum()*de),window_energy_J=direct,
              parseval_relative=error,peak_energy_eV=float(energy[np.flatnonzero(band)[np.argmax(density[band])]]),
              PIC_proxy_over_GENESIS=float(p['results'][label]['E_transverse_band_proxy_J']/(density[band].sum()*de)),
              PIC_native_flux_over_GENESIS=float(p['results'][label]['native_Poynting_band_J']/(density[band].sum()*de)))
            arr[label+'_J_per_eV']=density
        np.savez_compressed(out/(name+'-genesis.npz'),**arr)
        row=dict(plane=name,z_m=position,dump=str(path),native_all_pulse_J=full,
          matched_aperture_all_pulse_J=same,matched_over_native_aperture=same/full,
          captured_pulse_energy_fraction=stats['rect']['window_energy_J']/same,
          normalization_projection_relative_error=normalization_error,results=stats)
        rows.append(row);print('GENESIS',name,row['results']['rect'],flush=True)
    save(out/'genesis.json',rows);return rows

def ledger(root,out,repo):
    helpers=repo/'verification/performance/cases/radiation_method_validation'
    sys.path.insert(0,str(helpers))
    from report_discrete import load,calculate
    from relative_ledger import construct,write_tables,plot
    rows,sums,ranks,steps=load(root/'boris/ledger',root/'boris/ledger/energy-ledger.h5')
    index,q,result=calculate(rows,sums,steps)
    absolute,pct,meta=construct(rows,sums,index,q)
    directory=out/'ledger';write_tables(directory,index,absolute,pct,meta)
    plot(directory,index,pct,meta,'Explicit SPARC short-bunch energy accounting')
    save(out/'ledger.json',dict(summary=result,normalization=meta,mpi_ranks=ranks,steps=steps))
    return result,meta

def main():
    root=Path(sys.argv[1]);out=Path(sys.argv[2]);repo=Path(sys.argv[3]);out.mkdir(exist_ok=True)
    need((root/'boris/PIC_COMPLETE.json').exists(),'Explicit process not closed')
    need(json.loads((root/'genesis/COMPLETE.json').read_text())['run']['exit_code']==0,'Reference not closed')
    plan=json.loads((root/'plan.json').read_text())
    # Optional particle and accounting diagnostics are deliberately opt-in.
    # The current user-facing scope is radiation power.
    extended='--with-particle-ledger' in sys.argv[4:]
    if extended:p=particles(root,out,plan);print('Particles complete',flush=True)
    pic=[pic_plane(root,out,plan,n) for n,z in plan['geometry']['field_planes']]
    save(out/'pic.json',pic);gen=genesis(root,out,plan,pic)
    if extended:
        budget,normalization=ledger(root,out,repo)
        save(out/'ADDITIONAL_DIAGNOSTICS.json',dict(particles=p,ledger_summary=budget,normalization=normalization))
    save(out/'COMPLETE.json',dict(complete=True,physical_acceptance=False,
      pic=pic,genesis=gen,
      scope='Explicit version only, one common noise realization; total lab fields, no near-field subtraction'))
if __name__=='__main__':main()
