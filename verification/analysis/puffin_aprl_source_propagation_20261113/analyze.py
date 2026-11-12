"""One closed w64 source: current moments, transverse kernels and readout.

Source norms are diagnostics, never substituted for radiation energy.
"""
import concurrent.futures
import csv
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c,epsilon_0
from scipy.integrate import simpson
from scipy.optimize import brentq
from scipy.sparse import coo_matrix

W=Path('/public/home/wfollowm/APRL_Workspace')
PIC=W/'results/sparc_random_startup_20261112_r2/w64-a'
PUF=W/'results/puffin_same_sample_20261113_r3'
PREV=W/'results/sparc_single_seed_gap_20261113_r1'
LD=np.longdouble


def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def need(v,m):
    if not v:raise RuntimeError(m)
def table(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def puffin_dump(kind,step):
    files=[]
    for p in (PUF/'short-u07').glob(f'main_{kind}_*.h5'):
        with h5py.File(p,'r') as f:
            if int(f['runInfo'].attrs['iCsteps'])==step:files.append(p)
    need(bool(files),'Missing physical Puffin step')
    # The first dump at a step precedes matchOut/correctTrans bookkeeping.
    return min(files,key=lambda p:int(p.stem.rsplit('_',1)[1]))


def snapshot(program,plane):
    meta=read(PUF/'INPUT_VALIDATED.json');sc=meta['scale'];ref=read(PIC/'GEOMETRY.json')['reference_time_s']
    z=.01 if plane=='entry' else .206
    if program=='PIC':
        filename='und-entry-particles.h5' if plane=='entry' else 'u07-particles.h5'
        with h5py.File(PIC/'detectors'/filename,'r') as f:
            g=f['particle_plane'];need(int(g['complete'][()])==1,'PIC particles not closed');a=g['records'][:]
        xy=a['position_m'][:,:2];u=a['proper_velocity'];t=a['time_s'];q=abs(a['charge_C'])
    else:
        filename=puffin_dump('electrons',0 if plane=='entry' else 560)
        with h5py.File(filename,'r') as f:
            a=f['electrons'][:];distance=float(f['runInfo'].attrs['zbarTotal'])*float(f['runInfo'].attrs['Lg'])
        need(abs(distance-(z-.01))<1e-10,'Wrong Puffin particle position')
        xy=a[:,:2]*sc['transverse_scale_m'];ux=a[:,3]*sc['aw_peak'];uy=-a[:,4]*sc['aw_peak']
        gamma=a[:,5]*sc['gamma_ref'];uz=np.sqrt(gamma*gamma-1-ux*ux-uy*uy);u=np.column_stack((ux,uy,uz))
        t=meta['arrival_origin_s']-ref+(a[:,2]*sc['cooperation_length_m']+distance)/sc['code_c_m_per_s']
        q=a[:,6]*sc['npkbar']*sc['qe_C']
    need(len(q)==meta['records'] and abs(q.sum()/meta['total_charge_C']-1)<2e-12,'Source population/charge differs')
    return xy,u,np.asarray(t,dtype=LD)-LD(z)/LD(c),q


def source_maps(out,frequency):
    dx=4e-5;n=65;origin=-(n//2)*dx;xygrid=(np.arange(n)-n//2)*dx
    summaries=[];moments=[]
    for plane in ['entry','u07']:
        for program in ['PIC','Puffin']:
            xy,u,time,q=snapshot(program,plane)
            coords=(xy-origin)/dx;ij=np.floor(coords).astype(int);frac=coords-ij
            need(np.all(ij>=0) and np.all(ij<n-1),'Source outside audit grid')
            columns=np.arange(len(q));row=[];col=[];val=[]
            for a in (0,1):
                for b in (0,1):
                    row.append((ij[:,1]+b)*n+ij[:,0]+a);col.append(columns)
                    val.append((frac[:,0] if a else 1-frac[:,0])*(frac[:,1] if b else 1-frac[:,1]))
            deposition=coo_matrix((np.concatenate(val),(np.concatenate(row),np.concatenate(col))),shape=(n*n,len(q))).tocsr()
            # Density proxy, quiver-density proxy and actual transverse J/vz.
            weights=np.column_stack((q,q/u[:,2],q*u[:,0]/u[:,2],q*u[:,1]/u[:,2]))
            gaussian=np.column_stack((np.ones(len(q)),np.exp(-np.sum(xy**2,axis=1)/(4*.0001**2)),
                                      np.exp(-np.sum(xy**2,axis=1)/(4*.00025**2)),
                                      np.exp(-np.sum(xy**2,axis=1)/(4*.0004**2))))
            direct_weights=np.stack([weights[:,j,None]*gaussian for j in range(4)],axis=-1).reshape(len(q),-1)
            maps=np.zeros((len(frequency),n*n,4),complex);direct=np.zeros((len(frequency),4,4),complex)
            def block(indices):
                phase=np.exp(-2j*np.pi*np.asarray(frequency[indices,None].astype(LD)*time[None,:],float)).T
                result=np.stack([deposition@(phase*weights[:,j,None]) for j in range(4)],axis=-1)
                exact=phase.T@direct_weights
                return indices,result.transpose(1,0,2),exact.reshape(len(indices),4,4)
            indices=[np.arange(i,min(i+8,len(frequency))) for i in range(0,len(frequency),8)]
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                for index,value,exact in pool.map(block,indices):maps[index]=value;direct[index]=exact
            need(np.max(abs(maps.sum(axis=1)-direct[:,0,:]))<1e-12*q.sum(),'Spatial source sum not conserved')
            np.savez_compressed(out/f'source-{program}-{plane}.npz',frequency_Hz=frequency,xy_m=xygrid,
                                CIC_maps_C=maps.reshape(len(frequency),n,n,4),direct_modes_C=direct,
                                components=np.array(['charge','inverse_uz','Jx_over_vz','Jy_over_vz']))
            mid=len(frequency)//2
            summaries.append(dict(program=program,plane=plane,charge_C=float(q.sum()),
                mean_ux=float(np.average(u[:,0],weights=q)),mean_uy=float(np.average(u[:,1],weights=q)),
                mean_uz=float(np.average(u[:,2],weights=q)),sigma_x_m=float(np.std(xy[:,0])),sigma_y_m=float(np.std(xy[:,1])),
                global_bf=float(abs(direct[mid,0,0])/q.sum()),
                CIC_Jx_band_norm=float(simpson((abs(maps[:,:,2])**2).sum(axis=1)/dx**2,x=frequency)),
                CIC_quiver_density_band_norm=float(simpson((abs(maps[:,:,1])**2).sum(axis=1)/dx**2,x=frequency))))
            for mode,label in enumerate(['global','Gaussian100um','Gaussian250um','Gaussian400um']):
                for component,name in enumerate(['charge','inverse_uz','Jx_over_vz','Jy_over_vz']):
                    moments.append(dict(program=program,plane=plane,mode=label,source_component=name,
                        band_source_norm=float(simpson(abs(direct[:,mode,component])**2,x=frequency)),
                        carrier_real=float(direct[mid,mode,component].real),carrier_imag=float(direct[mid,mode,component].imag)))
            print('SOURCE_DONE',summaries[-1],flush=True)
    table(out/'source_summary.csv',summaries);table(out/'direct_source_modes.csv',moments)


def integrate(spec,x,y,frequency,half=.00056):
    dx=x[1]-x[0];dy=y[1]-y[0]
    wx=np.maximum(0,np.minimum(x+dx/2,half)-np.maximum(x-dx/2,-half))/dx
    wy=np.maximum(0,np.minimum(y+dy/2,half)-np.maximum(y-dy/2,-half))/dy
    density=(abs(spec)**2*wy[None,:,None,None]*wx[None,None,:,None]).sum(axis=(1,2,3))*dx*dy
    return float(simpson(density,x=frequency)),density


def measure(out,frequency):
    sc=read(PUF/'INPUT_VALIDATED.json')['scale'];d=PUF/'short-u07'
    with h5py.File(puffin_dump('aperp',560),'r') as f:
        r=dict(f['runInfo'].attrs);a=f['aperp'][:];nz=a.shape[1];ng=a.shape[2]
        need(abs(float(r['zbarTotal'])*float(r['Lg'])-.196)<1e-11,'Wrong field plane')
        dx=float(r['sLengthOfElmX'])*sc['transverse_scale_m'];dy=float(r['sLengthOfElmY'])*sc['transverse_scale_m']
        dt=float(r['sLengthOfElmZ2'])*float(r['Lc'])/sc['code_c_m_per_s']
    need(np.isfinite(a).all(),'Puffin field nonfinite')
    kernel=np.exp(-2j*np.pi*frequency[:,None]*np.arange(nz)*dt)*dt
    # Imaginary A is -Ey, consistent with the code's px-i py convention.
    p=np.empty((len(frequency),ng,ng,2),complex)
    for j,sign in [(0,1),(1,-1)]:p[:,:,:,j]=sign*(kernel@a[j].reshape(nz,-1)).reshape(len(frequency),ng,ng)
    p*=np.sqrt(2*sc['eps0']*sc['code_c_m_per_s'])*sc['field_scale_V_per_m']
    del a,kernel
    x=(np.arange(ng)-(ng-1)/2)*dx;y=(np.arange(ng)-(ng-1)/2)*dy
    with np.load(PREV/'analysis/u07-PIC.npz') as f:
        xp=f['x_m'];yp=f['y_m'];fp=f['frequency_Hz'];select=np.searchsorted(fp,frequency)
        need(np.array_equal(fp[select],frequency),'Frequency grid changed')
        pic=np.stack((f['forward_x_spectrum'][select],f['forward_y_spectrum'][select]),axis=-1)
    rows=[];density_rows=[];mode_rows=[]
    ep,sp=integrate(p,x,y,frequency);rows.append(dict(stage='Puffin native',band_energy_J=ep))
    ei,si=integrate(pic,xp,yp,frequency);rows.append(dict(stage='PIC stored detector',band_energy_J=ei))
    kx=2*np.pi*np.fft.fftfreq(ng,dx);ky=2*np.pi*np.fft.fftfreq(ng,dy)
    mx=np.exp(1j*(xp-x[0])[:,None]*kx[None,:])/ng
    my=np.exp(1j*(yp-y[0])[:,None]*ky[None,:])/ng
    spectrum=np.fft.fft2(p,axes=(1,2));variants={}
    for name,kind in [('Fourier positions only','none'),('PIC transverse detector averaging','detector'),
                      ('Both-axis linear sampling','bilinear')]:
        q=spectrum.copy()
        if kind=='detector':
            q[:,:,:,0]*=np.cos(ky*dy/2)[None,:,None]
            q[:,:,:,1]*=np.cos(kx*dx/2)[None,None,:]
        if kind=='bilinear':q*=np.cos(ky*dy/2)[None,:,None,None]*np.cos(kx*dx/2)[None,None,:,None]
        v=np.empty_like(pic)
        for j in range(2):v[:,:,:,j]=(my@q[:,:,:,j])@mx.T
        energy,density=integrate(v,xp,yp,frequency)
        rows.append(dict(stage=name,band_energy_J=energy));variants[name]=v
        density_rows.append((name,density))
    for name,v,xx,yy in [('Puffin native',p,x,y),('PIC detector',pic,xp,yp),
                        ('Puffin with PIC transverse readout',variants['PIC transverse detector averaging'],xp,yp)]:
        X,Y=np.meshgrid(xx,yy);ddx=xx[1]-xx[0];ddy=yy[1]-yy[0]
        mask=(abs(X)<.00056-1e-12)&(abs(Y)<.00056-1e-12)
        for sig in [.0001,.00025,.0004,1000.]:
            mode=np.ones_like(X) if sig>1 else np.exp(-(X*X+Y*Y)/(4*sig*sig))/(np.sqrt(2*np.pi)*sig)
            amp=np.sum(v[:,:,:,0]*mode*mask,axis=(1,2))*ddx*ddy
            mode_rows.append(dict(program=name,sigma_m=sig,band_projected_amplitude_squared=float(simpson(abs(amp)**2,x=frequency))))
    np.savez_compressed(out/'fields.npz',frequency_Hz=frequency,x_m=x,y_m=y,Puffin=p,PIC_x_m=xp,PIC_y_m=yp,
                        PIC=pic,Puffin_same_pixels=variants['Fourier positions only'],
                        Puffin_transverse_detector=variants['PIC transverse detector averaging'])
    table(out/'field_readout_stages.csv',rows);table(out/'field_modes.csv',mode_rows)
    save(out/'MEASUREMENT_RESULT.json',dict(complete=True,stages=rows,modes=mode_rows,
        scope='Known transverse readout applied to the same stored Puffin field, not a new simulation or output correction.',
        omitted=['Longitudinal interpolation and native B half-time','Source deposition','Particle feedback','Different propagators']))
    print('MEASUREMENT_DONE',rows,flush=True)


def propagation(out):
    plan=read(PIC.parent/'PLAN.json');geo=read(PIC/'GEOMETRY.json');G=plan['boost_gamma'];beta=np.sqrt(1-1/G**2)
    dz=geo['dz'];dx=dy=4e-5;rx=(dz/dx)**2;ry=rx;den=rx*ry+rx+ry
    b=np.array([rx,ry,1])*(1-rx*ry/den)/8
    dxy=rx*ry*(1/16-rx*ry/(8*den));dyz=ry*(1/16-ry/(8*den));dzx=rx*(1/16-rx/(8*den))
    alpha=np.array([1-2*b[1]-2*b[2]-4*dyz,1-2*b[2]-2*b[0]-4*dzx,1-2*b[0]-2*b[1]-4*dxy])
    klab=2*np.pi/plan['wavelength_m'];a=klab*dz/G;rows=[]
    for axis in ['x','y','diagonal']:
        for kt in np.linspace(0,60000,121):
            kx=kt if axis=='x' else (0 if axis=='y' else kt/np.sqrt(2))
            ky=kt if axis=='y' else (0 if axis=='x' else kt/np.sqrt(2))
            Kx=2*np.sin(kx*dx/2)/dx;Ky=2*np.sin(ky*dy/2)/dy;cx=np.cos(kx*dx);cy=np.cos(ky*dy)
            def dispersion(q):
                cz=np.cos(q)
                symbols=[alpha[0]+2*b[1]*cy+2*b[2]*cz+4*dyz*cy*cz,
                         alpha[1]+2*b[2]*cz+2*b[0]*cx+4*dzx*cz*cx,
                         alpha[2]+2*b[0]*cx+2*b[1]*cy+4*dxy*cx*cy]
                return np.sin((a-beta*q)/2)**2-dz*dz/4*(Kx*Kx*symbols[0]+Ky*Ky*symbols[1]+(2*np.sin(q/2)/dz)**2*symbols[2])
            q_exact=G*(np.sqrt(klab*klab-kt*kt)-beta*klab)*dz
            q=brentq(dispersion,q_exact-.4,q_exact+.4,xtol=1e-14)
            kz_num=beta*klab+q/(G*dz);kz_exact=np.sqrt(klab*klab-kt*kt);kz_par=klab-kt*kt/(2*klab)
            rows.append(dict(axis=axis,k_transverse_per_m=float(kt),boosted_phase_per_cell=float(q),
                PIC_minus_exact_phase_7period_rad=float((kz_num-kz_exact)*.196),
                Puffin_paraxial_minus_exact_phase_7period_rad=float((kz_par-kz_exact)*.196),
                F3_radiation_source_power=float(np.cos(q/2)**4),
                longitudinal_linear_source_power=float(np.sinc(q/(2*np.pi))**4)))
    table(out/'vacuum_dispersion.csv',rows)
    save(out/'PROPAGATION_SCOPE.json',dict(complete=True,
        scope='Exact homogeneous infinite-grid Cowan/Yee Fourier dispersion at fixed physical lab frequency, no boundary/source residue included.',
        source_files=['src/eb_field.cpp:forZDispersion','src/eb_field.cpp:smoothedElectric','src/eb_field.cpp:advanceMagneticCowan','src/eb_field.cpp:advanceElectric'],
        amplitude_attribution=False))


def main(root):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use bounded read-only allocation')
    need(read(PUF/'COMPLETE.json')['complete'] and read(PIC/'COMPLETE.json')['complete'],'Parent not complete')
    out=Path(root)/'analysis';out.mkdir()
    with np.load(PREV/'analysis/initial_sources.npz') as f:freq=f['frequency_Hz'][::2]
    measure(out,freq);source_maps(out,freq[::2]);propagation(out)
    save(out/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],new_seeds=0,new_simulations=0))


if __name__=='__main__':main(sys.argv[1])
