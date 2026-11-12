"""Transverse-basis probe on the same closed PIC crossing events."""
import os
import sys
from pathlib import Path
import numpy as np
import analyze
from analyze import snapshot,save,table


def main(root):
    out=Path(root)/'kernel-analysis';out.mkdir()
    xy,u,t,q=snapshot('PIC','u07');dx=4e-5;n=65;origin=-(n//2)*dx
    co=(xy-origin)/dx;index=np.floor(co).astype(int);fraction=co-index
    lam=analyze.read(analyze.PIC.parent/'PLAN.json')['wavelength_m'];f0=analyze.c/lam
    weight=q*u[:,0]/u[:,2]*np.exp(-2j*np.pi*np.asarray(t*f0,float))
    maps={}
    for kind in ['CIC_xy','Whitney_Jx_NGP_x_CIC_y']:
        value=np.zeros(n*n,complex)
        for a in ([0,1] if kind=='CIC_xy' else [0]):
            for b in [0,1]:
                w=(fraction[:,1] if b else 1-fraction[:,1])
                if kind=='CIC_xy':w=w*(fraction[:,0] if a else 1-fraction[:,0])
                cell=(index[:,1]+b)*n+index[:,0]+a
                value+=np.bincount(cell,weights=w*weight.real,minlength=n*n)+1j*np.bincount(cell,weights=w*weight.imag,minlength=n*n)
        maps[kind]=value.reshape(n,n)
    k=2*np.pi*np.fft.fftfreq(n,dx);kx,ky=np.meshgrid(k,k);kt=np.hypot(kx,ky)
    a=np.fft.fft2(maps['CIC_xy']);b=np.fft.fft2(maps['Whitney_Jx_NGP_x_CIC_y'])*np.exp(-1j*kx*dx/2)
    rows=[]
    for upper in [10000,20000,30000,40000,60000,1e9]:
        use=kt<=upper;pa=float(np.sum(abs(a[use])**2));pb=float(np.sum(abs(b[use])**2));cross=np.vdot(a[use],b[use])/np.sqrt(pa*pb)
        rows.append(dict(k_perp_max_m_inverse=upper,NGP_CIC_over_CIC_CIC_source_power=pb/pa,complex_coherence=float(abs(cross))))
    table(out/'transverse_source_basis.csv',rows)
    save(out/'RESULT.json',dict(complete=True,rows=rows,F3_on_axis_source_power=float(np.cos(np.pi/20)**4),
        scope='Exact transverse basis evaluation in the infinitesimal-chord limit on the same stored PIC events, carrier only.',
        limitations=['Not actual per-Maxwell-step 4D current; finite particle chords, longitudinal/time deposition and feedback are not reconstructed.',
                     'Source power ratios are not radiation energy correction factors.']))


if __name__=='__main__':
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Use bounded read-only allocation')
    main(sys.argv[1])
