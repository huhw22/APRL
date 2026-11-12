"""One closed w64 source: native slicing, common-frequency fields and apertures.

No seed sweep, input changes, physical simulation or amplitude fitting.
"""
import concurrent.futures
import csv
import json
import os
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c,e,h,epsilon_0
from scipy.integrate import simpson
from numpy.polynomial.legendre import leggauss

LD=np.longdouble
TWO=2*np.arccos(LD(-1))
W=Path('/public/home/wfollowm/APRL_Workspace')


def read(path):return json.loads(Path(path).read_text())
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')
def need(value,message):
    if not value:raise RuntimeError(message)
def table(path,rows):
    with Path(path).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def exact_block(omega,timevalue,weights):
    phase=np.asarray(-omega[:,None]*timevalue[None,:],dtype=float)
    values=np.cos(phase)+1j*np.sin(phase)
    return values@weights


def initial_source(root,out,frequency):
    p=read(root/'PLAN.json');d=root/'w64-a';info=read(d/'INPUT_VALIDATED.json');mapping=read(d/'NATIVE_MAPPING.json')
    lam=LD(p['wavelength_m']);omega=TWO*np.asarray(frequency,dtype=LD);omega0=TWO*LD(c)/lam
    anchor=LD(read(d/'GEOMETRY.json')['reference_time_s'])
    with h5py.File(info['input'],'r') as f:
        need(int(f['generation/complete'][()])==1,'Source incomplete');a=f['particles/records'][:]
    need(len(a)==info['macro_records'],'Wrong source population')
    raw=a['arrival_time_offset_s'].astype(LD);timevalue=raw-anchor
    u=a['proper_velocity'];gamma=np.sqrt(1+np.sum(u*u,axis=1))
    weights=np.column_stack([np.ones(len(a)),1/gamma,
        np.exp(-np.sum(a['plane_position_m']**2,axis=1)/(4*.0001**2))/gamma])
    norm=weights.sum(axis=0);exact=np.zeros((len(frequency),3),complex)
    groups=np.array_split(np.arange(len(frequency)),8)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for start in range(0,len(a),32768):
            stop=min(start+32768,len(a))
            futures=[pool.submit(exact_block,omega[index],timevalue[start:stop],weights[start:stop]) for index in groups]
            for index,future in zip(groups,futures):exact[index]+=future.result()
    carrier_phase=np.exp(-1j*np.asarray(omega0*timevalue,dtype=float))
    projected=np.zeros_like(exact);bucket_rows=[];native_ids=[]
    with h5py.File(d/'common.par.h5','r') as f:
        s0=LD(f['refposition'][0]);ds=LD(f['slicespacing'][0])
        for index in mapping['active_groups']:
            g=f[f'slice{index:06d}'];ids=g['source_id'][:];rows=np.searchsorted(a['source_id'],ids)
            need(np.array_equal(a['source_id'][rows],ids),'Native IDs differ')
            z=np.exp(1j*g['theta'][:])*np.exp(1j*float(np.remainder(omega0*anchor,TWO)))
            need(np.max(abs(z-carrier_phase[rows]))<5e-11,'Native carrier phase differs')
            moments=z@weights[rows]
            tb=-(s0+LD(index-1)*ds)/LD(c)
            projected+=np.exp(-1j*np.asarray((omega-omega0)*(tb-anchor),dtype=float))[:,None]*moments
            bucket_rows.append(dict(native_slice=index,records=len(ids),carrier_sum_real=float(moments[0].real),carrier_sum_imag=float(moments[0].imag)))
            native_ids.extend(ids.tolist())
    need(np.array_equal(np.sort(native_ids),a['source_id']),'Native coverage differs')
    carrier=len(frequency)//2
    err=float(np.max(abs(exact[carrier]-projected[carrier])/norm));need(err<5e-11,'Carrier projection mismatch')
    phase0=read(root/'PLAN.json')['source_config']['phase_space'];profile=p['source_config']['profile']
    flight=LD('.01')*LD(phase0['gamma_mean'])/(LD(c)*np.sqrt(LD(phase0['gamma_mean'])**2-1))
    start=LD(profile['window_start_s'])-flight;end=LD(profile['window_stop_s'])-flight
    aligned_origin=float(np.remainder(-LD(c)*start/lam,1))
    alternative=np.zeros_like(exact)
    coordinate=-LD(c)*raw/lam;integers=np.floor(coordinate-LD(aligned_origin)).astype(np.int64)
    for integer in np.unique(integers):
        use=integers==integer;moment=carrier_phase[use]@weights[use]
        tb=-(LD(integer)+LD(aligned_origin))*lam/LD(c)
        alternative+=np.exp(-1j*np.asarray((omega-omega0)*(tb-anchor),dtype=float))[:,None]*moment
    exact/=norm;projected/=norm;alternative/=norm
    # Analytic parabolic mean under a declared reference straight-flight map.
    def smooth_mean(order):
        nodes,w=leggauss(order);pieces=[];moments=[];times=[];total=LD(0)
        bmin=int(np.floor(-LD(c)*end/lam));bmax=int(np.floor(-LD(c)*start/lam))
        for bucket in range(bmin,bmax+1):
            lo=max(start,-LD(bucket+1)*lam/LD(c));hi=min(end,-LD(bucket)*lam/LD(c))
            if hi<=lo:continue
            t=(hi+lo)/2+(hi-lo)/2*nodes.astype(LD)
            tm=t+flight;T=LD(profile['full_duration_s']);density=4*(tm/T)*(1-tm/T)
            mass=w.astype(LD)*(hi-lo)/2*density
            pieces.append(np.column_stack([t,mass]));total+=mass.sum()
            moments.append(np.sum(mass*np.exp(-1j*np.asarray(omega0*(t-anchor),dtype=float))))
            times.append(-LD(bucket)*lam/LD(c))
        points=np.concatenate(pieces)
        exact_mean=np.exp(-1j*np.asarray(omega[:,None]*(points[:,0]-anchor),dtype=float))@np.asarray(points[:,1],float)/float(total)
        projection_mean=np.exp(-1j*np.asarray((omega-omega0)[:,None]*(np.asarray(times)-anchor),dtype=float))@np.asarray(moments,complex)/float(total)
        return exact_mean,projection_mean
    mean,mean_projected=smooth_mean(24);mean2,mean_projected2=smooth_mean(48)
    need(max(np.max(abs(mean-mean2)),np.max(abs(mean_projected-mean_projected2)))<2e-11,'Analytic mean quadrature failed')
    rows=[]
    for j,label in enumerate(['unit','inverse_gamma','Gaussian100um_inverse_gamma']):
        for key,value in [('actual-native',projected),('aligned-buckets-offline',alternative)]:
            q=float(simpson(abs(exact[:,j])**2,x=frequency)/simpson(abs(value[:,j])**2,x=frequency))
            rows.append(dict(weight=label,projection=key,exact_over_projected_band_factor=q,
                scope='Initial source factor only; no trajectory, field Green tensor or output-power correction'))
    table(out/'initial_source_ratios.csv',rows);table(out/'native_buckets.csv',bucket_rows)
    np.savez_compressed(out/'initial_sources.npz',frequency_Hz=frequency,exact=exact,projected=projected,
        aligned_projected=alternative,analytic_mean=mean,analytic_projected_mean=mean_projected)
    save(out/'SOURCE_RESULT.json',dict(complete=True,case='w64-a',records=len(a),carrier_B_difference=err,
        aligned_origin_wavelength=aligned_origin,actual_support_fractional_edges=np.remainder(-LD(c)*np.array([raw.min(),raw.max()])/lam,1).astype(float),
        source_factors=rows,analytic_mean_scope='Parabolic mean with nominal on-axis reference flight; exact random source always uses original actual event times',
        analytic_quadrature_error=float(max(np.max(abs(mean-mean2)),np.max(abs(mean_projected-mean_projected2)))),
        analytic_mean_exact_over_projected_band=float(simpson(abs(mean)**2,x=frequency)/simpson(abs(mean_projected)**2,x=frequency)),
        no_input_changes=True,no_amplitude_fit=True,new_simulations=0))
    print('SOURCE',rows,flush=True)


def measure_fields(root,out,frequency):
    p=read(root/'PLAN.json');d=root/'w64-a';f0=c/p['wavelength_m'];measurements=[];mode_rows=[]
    for plane,dumpname in [('u07','reference.29.fld.h5'),('u77','reference.309.fld.h5')]:
        with h5py.File(d/'detectors'/(plane+'.h5'),'r') as f:
            g=f['field_plane'];need(int(g['complete'][()])==1,'PIC field not complete')
            t=g['time_s'][:];data=g['fields'];dt=float(np.median(np.diff(t)))
            ticks=np.rint((t-t[0])/dt).astype(int);N=int(ticks[-1]+1)
            need(np.max(abs((t-t[0])/dt-ticks))<1e-4,'Off-grid times')
            dx=float(g.attrs['dx_m']);dy=float(g.attrs['dy_m']);nx=data.shape[2];ny=data.shape[1]
            x=float(g.attrs['x_first_m'])+np.arange(nx)*dx;y=float(g.attrs['y_first_m'])+np.arange(ny)*dy
            ix=np.flatnonzero(abs(x)<.00088-dx*.1);iy=np.flatnonzero(abs(y)<.00088-dy*.1)
            cube=np.empty((N,len(iy),len(ix),4))
            for start in range(0,len(t),32):
                block=data[start:start+32];stop=start+len(block)
                cube[ticks[start:stop],:,:,:2]=block['electric_V_per_m'][:,iy][:,:,ix,:2]
                cube[ticks[start:stop],:,:,2:]=c*block['magnetic_T'][:,iy][:,:,ix,:2]
            missing=np.setdiff1d(np.arange(N),ticks);need(len(missing)<=1,'Multiple missing samples')
            for k in missing:
                nodes=ticks[np.argsort(abs(ticks-k))[:8]]
                coefficients=np.array([np.prod([(k-j)/(i-j) for j in nodes if j!=i]) for i in nodes])
                cube[k]=np.tensordot(coefficients,cube[nodes],axes=1)
        kernel=np.exp(-2j*np.pi*frequency[:,None]*(np.arange(N)*dt))*dt
        spectrum=(kernel@cube.reshape(N,-1)).reshape(len(frequency),len(iy),len(ix),4)
        del cube,kernel
        ex,ey,bx,by=[spectrum[:,:,:,j] for j in range(4)]
        plusx=(ex+by)/2;plusy=(ey-bx)/2;minusx=(ex-by)/2;minusy=(ey+bx)/2
        norm=2*epsilon_0*c
        F=norm*(abs(plusx)**2+abs(plusy)**2);B=norm*(abs(minusx)**2+abs(minusy)**2)
        P=norm*(ex*by.conj()-ey*bx.conj()).real
        need(np.max(abs(P-(F-B)))<2e-13*max(F.max(),1e-300),'Forward/backward algebra fails')
        xp=x[ix];yp=y[iy];xx,yy=np.meshgrid(xp,yp)
        for half in [.00016,.00032,.00056,.0008]:
            mask=(abs(xx)<half-dx*.1)&(abs(yy)<half-dy*.1)
            integrand=P[:,mask].sum(axis=1)*dx*dy
            total=float(simpson(integrand,x=frequency));coarse=float(simpson(integrand[::2],x=frequency[::2]))
            measurements.append(dict(plane=plane,program='PIC',half_width_m=half,band_flux_J=total,
                forward_proxy_J=float(simpson(F[:,mask].sum(axis=1)*dx*dy,x=frequency)),
                backward_proxy_J=float(simpson(B[:,mask].sum(axis=1)*dx*dy,x=frequency)),
                common_frequency_quadrature_relative=abs(total-coarse)/max(abs(total),1e-300)))
        for sig in [.0001,.00025,.0004]:
            mode=np.exp(-(xx*xx+yy*yy)/(4*sig*sig))/(np.sqrt(2*np.pi)*sig)
            mask=(abs(xx)<.00056-dx*.1)&(abs(yy)<.00056-dy*.1)
            amplitude=np.sum(np.sqrt(norm)*plusx*mode*mask,axis=(1,2))*dx*dy
            mode_rows.append(dict(plane=plane,program='PIC-forward-x',sigma_m=sig,
                band_mode_energy_J=float(simpson(abs(amplitude)**2,x=frequency))))
        np.savez_compressed(out/(plane+'-PIC.npz'),x_m=xp,y_m=yp,frequency_Hz=frequency,
            forward_x_spectrum=np.sqrt(norm)*plusx,forward_y_spectrum=np.sqrt(norm)*plusy,
            band_forward_x_J_per_m2=simpson(norm*abs(plusx)**2,x=frequency,axis=0))
        del spectrum,ex,ey,bx,by,plusx,plusy,minusx,minusy,F,B,P
        with h5py.File(d/'genesis-off'/dumpname,'r') as f:
            ns=int(f['slicecount'][0]);ng=int(f['gridpoints'][0]);dx=float(f['gridsize'][0]);dt=float(f['slicespacing'][0])/c
            xy=(np.arange(ng)-(ng-1)/2)*dx;index=np.flatnonzero(abs(xy)<.0009)
            values=np.empty((ns,len(index),len(index)),complex)
            for j in range(ns):
                group=f[f'slice{j+1:06d}'];a=group['field-real'][:]+1j*group['field-imag'][:]
                values[j]=a.reshape(ng,ng)[np.ix_(index,index)].conj()/dx
            tb=-np.arange(ns)*dt
        kernel=np.exp(-2j*np.pi*(frequency-f0)[:,None]*tb)*dt
        spec=(kernel@values.reshape(ns,-1)).reshape(len(frequency),len(index),len(index))
        del values,kernel
        xx,yy=np.meshgrid(xy[index],xy[index])
        for half in [.00016,.00032,.00056,.0008]:
            w=np.maximum(0,np.minimum(xy[index]+dx/2,half)-np.maximum(xy[index]-dx/2,-half))/dx
            area=w[:,None]*w[None,:]*dx*dx
            integrand=np.sum(abs(spec)**2*area,axis=(1,2))
            total=float(simpson(integrand,x=frequency));coarse=float(simpson(integrand[::2],x=frequency[::2]))
            measurements.append(dict(plane=plane,program='GENESIS-off',half_width_m=half,band_flux_J=total,
                forward_proxy_J=total,backward_proxy_J=0.,common_frequency_quadrature_relative=abs(total-coarse)/max(abs(total),1e-300)))
        for sig in [.0001,.00025,.0004]:
            mode=np.exp(-(xx*xx+yy*yy)/(4*sig*sig))/(np.sqrt(2*np.pi)*sig)
            half=.00056;w=np.maximum(0,np.minimum(xy[index]+dx/2,half)-np.maximum(xy[index]-dx/2,-half))/dx
            amplitude=np.sum(spec*mode*w[:,None]*w[None,:],axis=(1,2))*dx*dx
            mode_rows.append(dict(plane=plane,program='GENESIS-off',sigma_m=sig,
                band_mode_energy_J=float(simpson(abs(amplitude)**2,x=frequency))))
        np.savez_compressed(out/(plane+'-GENESIS.npz'),x_m=xy[index],y_m=xy[index],frequency_Hz=frequency,
            forward_x_spectrum=spec,band_forward_x_J_per_m2=simpson(abs(spec)**2,x=frequency,axis=0))
        del spec
        print('FIELDS',plane,[x for x in measurements if x['plane']==plane and x['half_width_m']==.00056],flush=True)
    table(out/'field_measurements.csv',measurements);table(out/'gaussian_modes.csv',mode_rows)
    save(out/'FIELD_RESULT.json',dict(complete=True,measurements=measurements,gaussian_modes=mode_rows,
        scope='Same actual continuous frequency points, native spatial quadratures. Forward/backward are algebraic characteristics, not exact near-field separation.',
        no_time_crop=True,no_amplitude_fit=True,zero_padding_or_clock_fit=False,new_simulations=0))


def main(root,out):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use a bounded read-only Slurm allocation')
    root=Path(root);out=Path(out);out.mkdir(exist_ok=False)
    need(read(root/'w64-a/COMPLETE.json')['complete'],'Source campaign incomplete')
    f0=c/read(root/'PLAN.json')['wavelength_m'];frequency=f0*(1+np.linspace(-.05,.05,513))
    initial_source(root,out,frequency);measure_fields(root,out,frequency)
    save(out/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],case='w64-a',new_simulations=0,extra_seeds=0))


if __name__=='__main__':
    import sys
    main(sys.argv[1],sys.argv[2])
