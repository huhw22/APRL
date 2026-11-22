"""Generate native four-side OSIRIS startup for the actual THz beam and mesh."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from scipy import fft
import yaml

from startup_operator import CASE,CFG,NX,NY,NZ,L,H,DT,C,BETA,ALPHA,B,D,values,Grid,E_OFF,B_OFF
from startup_collar import masks,solve,NAMES

HERE=CASE/'startup'
SCRATCH=HERE/'scratch'
REF=CASE/'preflight/reference'
WINDOW=CFG['mesh']['moving_window']['cells']
K0=int(re.search(r'initial_offset=(\d+)',(CASE/'preflight/aprl.log').read_text())[1])
PAD=256
FIRST=K0-PAD
N=fft.next_fast_len(WINDOW+2*PAD+2)
assert FIRST>=1 and FIRST+N<NZ


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def taper(dist,width):
    return .5-.5*np.cos(np.pi*np.clip(dist/width,0,1))


def layouts():
    g=Grid(2*np.pi/N);me,mb=masks(g)
    mask={k:v.reshape(sh).T for k,v,sh in zip(NAMES,me+mb,g.eshapes+g.bshapes)}
    dims={k:int(v.sum()) for k,v in mask.items()}
    for idx,u in enumerate(g.F+g.A):
        if u.axis<2:
            sh=(g.bshapes if idx<6 else g.eshapes)[u.target]
            dims['psi'+str(idx)]=2*L*sh[1-u.axis]
    return mask,dims


def donor(phi,theta):
    def select(a,b):return phi[1+a:1+a+NX+1,1+b:1+b+NY+1]
    cen=select(0,0);sx=select(1,0)+select(-1,0);sy=select(0,1)+select(0,-1)
    sxy=select(1,1)+select(1,-1)+select(-1,1)+select(-1,-1)
    cz=2*np.cos(theta);shift,=values('translation')
    s0=(ALPHA[0]+B[2]*cz)*cen+(B[1]+D[0]*cz)*sy
    s1=(ALPHA[1]+B[2]*cz)*cen+(B[0]+D[1]*cz)*sx
    s2=(ALPHA[2]-shift)*cen+B[0]*sx+B[1]*sy+D[2]*sxy
    ex=-np.diff(s0,axis=0)/H[0];ey=-np.diff(s1,axis=1)/H[1]
    ez=-2j*np.sin(theta/2)/H[2]*s2
    w0,w1,scale=values('magnetic_interpolation')
    fac=C*scale*(w0+w1*np.exp(1j*theta))*np.exp(-.5j*(theta+BETA*theta))
    return [ex,ey,ez,-fac*ey,fac*ex,np.zeros((NX,NY),complex)]


def prepare():
    ready=json.loads((CASE/'preflight/REFERENCE_READY.json').read_text())
    assert ready['reference_complete'] and ready['expected_duration_guard'] and not ready['physical_run_complete']
    HERE.mkdir();SCRATCH.mkdir()
    p=np.memmap(REF/'phi.bin',dtype='<f8',mode='r',shape=(NZ-1,NY-1,NX-1))
    side=np.memmap(REF/'side-potential.bin',dtype='<f8',mode='r',shape=(NZ+1,4*(NX+NY+6)))
    out=np.lib.format.open_memmap(SCRATCH/'potential.npy',mode='w+',dtype='complex128',shape=(N//2+1,NX+3,NY+3))
    k=slice(FIRST,FIRST+N);m=np.arange(N);edge=taper(np.minimum(m,N-1-m),PAD//2)
    for j0 in range(-1,NY+2,8):
        y=np.arange(j0,min(j0+8,NY+2));block=np.zeros((N,len(y),NX+3))
        for iy,j in enumerate(y):
            if 0<j<NY:block[:,iy,2:-2]=p[FIRST-1:FIRST+N-1,j-1,:]
            else:
                slot=j+1 if j<=0 else j-NY+2
                off=4*(NY+3)+slot*(NX+3)
                block[:,iy,:]=side[k,off:off+NX+3]
            for slot,i in enumerate((-1,0,NX,NX+1)):
                block[:,iy,i+1]=side[k,slot*(NY+3)+j+1]
        out[:,:,y+1]=fft.rfft(block*edge[:,None,None],axis=0).transpose(0,2,1)
    out.flush()
    _,dims=layouts()
    for key,width in {**dims,**{k+'-donor':dims[k] for k in NAMES}}.items():
        x=np.lib.format.open_memmap(SCRATCH/(key+'.npy'),mode='w+',dtype='complex128',shape=(N//2+1,width))
        del x
    print('STARTUP FFT READY',N,'modes',N//2+1,'native offset',K0,flush=True)


def worker():
    rank=int(os.environ['SLURM_PROCID']);count=int(os.environ['SLURM_NTASKS'])
    ph=np.load(SCRATCH/'potential.npy',mmap_mode='r');mask,dims=layouts()
    specs={k:np.load(SCRATCH/(k+'.npy'),mmap_mode='r+') for k in list(dims)+[k+'-donor' for k in NAMES]}
    rows=[]
    for m in range(rank,N//2+1,count):
        tic=time.monotonic();theta=2*np.pi*m/N
        # Same finite conjugate-pair DC continuation as the 2.1 reference.
        solve_theta=theta if m else 2*np.pi/N/64
        g=Grid(solve_theta);vac=donor(ph[m],solve_theta)
        fld,mem,stats=solve(g,vac)
        normalizer=max(max(float(abs(v).max()) for v in fld),1e-99)
        assert stats['ampere_residual']/normalizer<1e-8,(m,stats)
        if m==0:
            neg,nmem,_=solve(Grid(-solve_theta),donor(ph[m],-solve_theta))
            assert max(float(abs(a-b.conj()).max()) for a,b in zip(fld,neg))/normalizer<1e-7
        for idx,name in enumerate(NAMES):
            z=(E_OFF+B_OFF)[idx][2]
            phase=np.exp(1j*solve_theta*z+(.5j*g.omega if idx>=3 else 0))
            unit=1 if idx<3 else 1/C
            out=phase*unit*fld[idx].T[mask[name]]
            raw=phase*unit*donor(ph[m],theta)[idx].T[mask[name]]
            specs[name][m]=out if m else out.real
            specs[name+'-donor'][m]=raw
        for idx,u in enumerate(g.F+g.A):
            if u.axis==2:continue
            ar=mem['psi'+str(idx)];sides=np.r_[0:L,ar.shape[u.axis]-L:ar.shape[u.axis]]
            slots=ar[sides,:] if u.axis==0 else ar[:,sides].T
            target=u.target+(3 if idx<6 else 0)
            phase=np.exp(1j*solve_theta*(E_OFF+B_OFF)[target][2])
            out=(slots*phase*(1 if idx<6 else 1/C)).ravel()
            specs['psi'+str(idx)][m]=out if m else out.real
        rows.append(dict(mode=m,seconds=time.monotonic()-tic,relative_residual=stats['ampere_residual']/normalizer))
        if len(rows)==1:print('FIRST MODE',rank,rows[-1],flush=True)
    (SCRATCH/f'worker-{rank}.json').write_text(json.dumps(rows)+'\n')
    print('MODES COMPLETE',rank,len(rows),flush=True)


def finish(count):
    rows=sum([json.loads((SCRATCH/f'worker-{r}.json').read_text()) for r in range(count)],[])
    assert sorted(r['mode'] for r in rows)==list(range(N//2+1))
    output=HERE/'arrays';output.mkdir();g=Grid(2*np.pi/N);_,dims=layouts();files={}
    for key,width in dims.items():
        hist=key.startswith('psi')
        idx=int(key[3:]) if hist else NAMES.index(key)
        target=(g.F+g.A)[idx].target+(3 if idx<6 else 0) if hist else idx
        z=(E_OFF+B_OFF)[target][2];nk=WINDOW+int(z==0)
        label=key if hist else key[0].upper()+key[1]
        spec=np.load(SCRATCH/(key+'.npy'),mmap_mode='r')
        base=None if hist else np.load(SCRATCH/(key+'-donor.npy'),mmap_mode='r')
        out=np.memmap(output/(label+'.bin'),mode='w+',dtype='<f8',shape=(nk,width))
        taper_z=taper(np.minimum(np.arange(nk)+z,WINDOW-np.arange(nk)-z),48)[:,None]
        keep=slice(K0-FIRST,K0-FIRST+nk)
        for low in range(0,width,128):
            sl=slice(low,min(low+128,width));x=fft.irfft(spec[:,sl],n=N,axis=0)[keep]
            if base is None:out[:,sl]=x*taper_z
            else:
                control=fft.irfft(base[:,sl],n=N,axis=0)[keep]
                out[:,sl]=control+taper_z*(x-control)
        out.flush();assert np.isfinite(out).all()
        files[label]=dict(shape=[nk,width],sha256=sha(output/(label+'.bin')),peak=float(abs(out).max()))
        print('EXPORTED',label,files[label],flush=True)
    boundary=CFG['boundary'];beam=CFG['beam'];mesh=CFG['mesh']
    names=list(files)
    fnv=CASE.parents[0]/'fnv1a64'
    hashes=subprocess.check_output([str(fnv),beam['input']['file']]+[str(output/(k+'.bin')) for k in names],text=True).splitlines()
    assert len(hashes)==len(names)+1
    manifest=dict(version=1,method='osiris',reference_cells=[NX,NY,NZ],reference_origin_m=values('origin').tolist(),
       spacing_m=H.tolist(),dt_s=float(DT),reference_beta=float(BETA),frame_gamma=mesh['boost_gamma'],
       window_cells=[NX,NY,WINDOW],reference_offset=K0,cpml_cells=boundary['cells'],
       target_reflection=boundary['target_reflection'],polynomial_order=boundary['polynomial_order'],
       kappa_max=boundary['kappa_max'],alpha_fraction=boundary['alpha_fraction'],
       beam_position_offset_m=beam['input']['position_offset'],beam_electrons=beam['input']['electrons'],
       beam_reference=beam['reference'],transverse_x_binomial3=True,longitudinal_shape=mesh['longitudinal_shape'],
       beam_fnv1a64=hashes[0],arrays_fnv1a64=dict(zip(names,hashes[1:])))
    (output/'startup.yaml').write_text(yaml.safe_dump(manifest,sort_keys=False))
    (output/'READY.json').write_text(json.dumps(dict(complete=True,modes=N//2+1,samples=N,native_offset=K0,
       maximum_relative_modal_residual=max(r['relative_residual'] for r in rows),files=files),indent=2)+'\n')
    print('FOUR SIDE STARTUP COMPLETE',flush=True)


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='worker':worker()
    elif sys.argv[1]=='finish':finish(int(sys.argv[2]))
