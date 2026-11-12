"""Closed-source timing and native four-node field interpolation diagnostics."""
import os
import sys
from pathlib import Path
import numpy as np
from scipy.constants import c
from scipy.integrate import quad,simpson
from analyze import read,save,table,snapshot,PIC,PUF


def main(root):
    out=Path(root)/'phase-analysis';out.mkdir();rows=[]
    for program in ['PIC','Puffin']:
        for plane in ['entry','u07']:
            xy,u,t,q=snapshot(program,plane)
            rows.append(dict(program=program,plane=plane,
                mean_retarded_time_s=float(np.sum(t*q.astype(np.longdouble))/q.sum()),
                mean_x_m=float(np.average(xy[:,0],weights=q)),mean_y_m=float(np.average(xy[:,1],weights=q)),
                mean_gamma=float(np.average(np.sqrt(1+np.sum(u*u,axis=1)),weights=q))))
    lam=read(PIC.parent/'PLAN.json')['wavelength_m'];gamma=298.2608843802298;K=2.060074987322757;ku=2*np.pi/.028
    def ux(z):
        h=.5*(1-np.cos(np.pi*z/.028)) if z<.028 else 1.
        return -K*h*np.cos(ku*z)
    def retard(z):
        transverse=1+ux(z)**2;uz=np.sqrt(gamma*gamma-transverse)
        return transverse/(uz*(gamma+uz))/c
    theoretical=quad(retard,0,.028,epsabs=1e-28)[0]+quad(retard,.028,.196,epsabs=1e-28,limit=300)[0]
    delays={}
    for program in ['PIC','Puffin']:
        a=next(x for x in rows if x['program']==program and x['plane']=='entry')
        b=next(x for x in rows if x['program']==program and x['plane']=='u07')
        delays[program]=b['mean_retarded_time_s']-a['mean_retarded_time_s']
    measured_phase=-2*np.pi*c/lam*(delays['Puffin']-delays['PIC'])
    phases=[]
    for plane in ['entry','u07']:
        with np.load(Path(root)/'analysis'/f'source-PIC-{plane}.npz') as a, \
             np.load(Path(root)/'analysis'/f'source-Puffin-{plane}.npz') as b:
            for j,label in enumerate(['charge','inverse_uz','Jx','Jy']):
                x=a['CIC_maps_C'][64,:,:,j];y=b['CIC_maps_C'][64,:,:,j]
                cross=np.vdot(x,y)/np.sqrt(np.vdot(x,x).real*np.vdot(y,y).real)
                phases.append(dict(plane=plane,component=label,coherence=float(abs(cross)),phase_rad=float(np.angle(cross))))
    table(out/'arrival_moments.csv',rows);table(out/'source_phase.csv',phases)
    save(out/'TIMING.json',dict(complete=True,delays_s=delays,
        on_axis_external_field_exact_reference_delay_s=theoretical,
        delta_source_phase_from_mean_arrival_rad=float(measured_phase),
        scope='Weighted real arrival moments; on-axis zero-emittance external-only timing is a reference, not an ensemble PIC trajectory prediction.'))
    with np.load(Path(root)/'analysis/fields.npz') as f:
        freq=f['frequency_Hz'];p=f['Puffin'];x=f['x_m'];y=f['y_m'];xp=f['PIC_x_m'];yp=f['PIC_y_m'];pic=f['PIC'];d=f['Puffin_transverse_detector']
    dx=x[1]-x[0];i=np.searchsorted(x,xp)-1;j=np.searchsorted(y,yp)-1
    fx=(xp-x[i])/dx;fy=(yp-y[j])/dx
    # Directly use the four native node values and actual interpolation weights.
    v=np.zeros_like(pic)
    for a in (0,1):
        for b in (0,1):
            v+=p[:,j[:,None]+b,i[None,:]+a,:]*(fy if b else 1-fy)[None,:,None,None]*(fx if a else 1-fx)[None,None,:,None]
    k=np.flatnonzero(abs(xp)<.00056-1e-12);l=np.flatnonzero(abs(yp)<.00056-1e-12)
    def energy(a):return float(simpson((abs(a[:,l[:,None],k[None,:],:])**2).sum(axis=(1,2,3))*dx**2,x=freq))
    ec=energy(v);ep=energy(pic);ef=energy(p[:,j[:,None],i[None,:]])
    coherence=[]
    for name,a in [('Puffin bilinear native evaluation',v),('Puffin PIC transverse readout',d)]:
        aa=a[:,l[:,None],k[None,:],0].reshape(len(freq),-1);bb=pic[:,l[:,None],k[None,:],0].reshape(len(freq),-1)
        z=np.sum(aa.conj()*bb,axis=1)/np.sqrt(np.sum(abs(aa)**2,axis=1)*np.sum(abs(bb)**2,axis=1))
        coherence.append(dict(stage=name,mean_spatial_coherence=float(np.mean(abs(z))),
                              carrier_spatial_coherence=float(abs(z[len(freq)//2]))))
    save(out/'INTERPOLATION.json',dict(complete=True,bilinear_native_Puffin_energy_J=ec,
        PIC_energy_J=ep,Puffin_over_PIC_after_same_point_native_evaluation=ec/ep,
        max_target_fraction_distance_from_half=float(max(np.max(abs(fx-.5)),np.max(abs(fy-.5)))),
        coherence=coherence,
        scope='Literal Puffin four-node transverse field evaluation at PIC detector coordinates; not an arbitrary smoothing or amplitude fit.',
        limitations=['Two solvers store different staggered grids; agreement of this scalar energy is not agreement of the full field.',
                     'Longitudinal/time sampling and propagation remain distinct.']))
    np.savez_compressed(out/'bilinear_field.npz',frequency_Hz=freq,x_m=xp,y_m=yp,Puffin_bilinear=v)
    print((out/'TIMING.json').read_text(),flush=True);print((out/'INTERPOLATION.json').read_text(),flush=True)


if __name__=='__main__':
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Use a bounded read-only allocation')
    main(sys.argv[1])
