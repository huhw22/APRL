"""Compare a closed Puffin short run with the existing same-sample fields."""
import csv
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson

from prepare import W,C,read,save
from validate import require,h5attr


def main(root):
    require(read(root/'TRACKING_CLOSED.json')['complete'],'Tracking is not closed')
    info=read(root/'INPUT_VALIDATED.json');scale=info['scale'];out=root/'analysis';out.mkdir()
    d=root/'short-u07';files=list(d.glob('main_aperp_*.h5'))
    field=max(files,key=lambda p:float(h5attr(p,'runInfo','zbarTotal')))
    with h5py.File(field,'r') as f:
        run=dict(f['runInfo'].attrs);g=f['aperp'];nx=int(run['nX']);ny=int(run['nY']);nz=int(run['nZ2'])
        lg=float(run['Lg']);lc=float(run['Lc']);trans=np.sqrt(lg*lc)
        dx=float(run['sLengthOfElmX'])*trans;dy=float(run['sLengthOfElmY'])*trans
        dt=float(run['sLengthOfElmZ2'])*lc/scale['code_c_m_per_s'];z=float(run['zbarTotal'])*lg+.01
        require(abs(z-.206)<1e-10,'Wrong physical observation plane')
        require(g.shape==(2,nz,ny,nx),'Unexpected Puffin field storage order')
        # Actual interpolation positions use (N-1)*dx/2; the old vizSchema
        # globalLimits can extend to N*dx/2 and are not used to shift data.
        x=(np.arange(nx)-(nx-1)/2)*dx;y=(np.arange(ny)-(ny-1)/2)*dy
        ix=np.flatnonzero(abs(x)<=.00056+1e-12);iy=np.flatnonzero(abs(y)<=.00056+1e-12)
        values=np.stack([g[j,:,:,ix[0]:ix[-1]+1][:,iy,:] for j in range(2)],axis=-1)
    require(np.isfinite(values).all(),'Nonfinite raw Puffin field')
    wx=np.maximum(0,np.minimum(x[ix]+dx/2,.00056)-np.maximum(x[ix]-dx/2,-.00056))/dx
    wy=np.maximum(0,np.minimum(y[iy]+dy/2,.00056)-np.maximum(y[iy]-dy/2,-.00056))/dy
    area=wy[:,None]*wx[None,:]*dx*dy
    norm=scale['eps0']*scale['code_c_m_per_s']*scale['field_scale_V_per_m']**2
    cached=W/'results/sparc_single_seed_gap_20261113_r1'
    with np.load(cached/'analysis/initial_sources.npz') as f:frequency=f['frequency_Hz']
    f0=frequency[len(frequency)//2]
    kernel=np.exp(-2j*np.pi*frequency[:,None]*np.arange(nz)*dt)*dt
    spec=(kernel@values.reshape(nz,-1)).reshape(len(frequency),len(iy),len(ix),2)
    density=2*norm*np.sum(abs(spec)**2*area[None,:,:,None],axis=(1,2,3))
    energy=float(simpson(density,x=frequency));coarse=float(simpson(density[::2],x=frequency[::2]))
    fft=np.fft.rfft(values,axis=0)*dt;factors=np.full(len(fft),2.)
    factors[0]=1
    if nz%2==0:factors[-1]=1
    transformed=float(norm*np.sum(abs(fft)**2*area[None,:,:,None]*factors[:,None,None,None])/(nz*dt))
    direct=float(norm*np.sum(values**2*area[None,:,:,None])*dt)
    require(abs(transformed/direct-1)<2e-12,'One-sided real-field Parseval failed')
    time_energy=np.sum(values**2*area[None,:,:,None],axis=(1,2,3))
    edge=float((time_energy[:40].sum()+time_energy[-40:].sum())/time_energy.sum())
    summary=[dict(program='Puffin matched external magnet',energy_J=energy)]
    spectra=[('Puffin',density)]
    for program in ['PIC','GENESIS']:
        with np.load(cached/'analysis'/('u07-'+program+'.npz')) as f:
            xp=f['x_m'];yp=f['y_m'];ddx=xp[1]-xp[0];ddy=yp[1]-yp[0]
            qx=np.maximum(0,np.minimum(xp+ddx/2,.00056)-np.maximum(xp-ddx/2,-.00056))/ddx
            qy=np.maximum(0,np.minimum(yp+ddy/2,.00056)-np.maximum(yp-ddy/2,-.00056))/ddy
            sd=abs(f['forward_x_spectrum'])**2
            if 'forward_y_spectrum' in f:sd+=abs(f['forward_y_spectrum'])**2
            sd=(sd*qy[None,:,None]*qx[None,None,:]).sum(axis=(1,2))*ddx*ddy
        spectra.append((program,sd));summary.append(dict(program=program,energy_J=float(simpson(sd,x=frequency))))
    half=read(cached/'half-bucket/RESULT.json')['comparisons'][0]
    summary.append(dict(program='GENESIS half-bucket',energy_J=half['half_bucket_band_J']))
    particle=max(d.glob('main_electrons_*.h5'),key=lambda p:float(h5attr(p,'runInfo','zbarTotal')))
    with h5py.File(particle,'r') as f:a=f['electrons'][:]
    require(len(a)==info['records'],'Particles removed or duplicated in tracking')
    charge=float(a[:,6].sum()*scale['npkbar']*scale['qe_C'])
    require(abs(charge/info['total_charge_C']-1)<2e-12,'Tracking changed total charge')
    bf=float(abs(np.sum(a[:,6]*np.exp(2j*np.pi*f0*a[:,2]*lc/scale['code_c_m_per_s']))/a[:,6].sum()))
    result=dict(complete=True,physical_acceptance=False,plane_m=z,band_half_fraction=.05,aperture_half_m=.00056,
        comparisons=summary,common_frequency_quadrature_relative=abs(energy-coarse)/energy,
        real_field_Parseval_relative=abs(transformed/direct-1),first_last_40_nodes_energy_fraction=edge,
        particle_count=len(a),charge_C=charge,Puffin_full_pulse_global_b1=bf,
        scope='Same sampled events, explicit matched external magnet; raw forward transverse Puffin fields vs PIC algebraic forward fields. No amplitude fitting.',
        limitations=['Single short run; no mesh or step convergence acceptance.',
                    'No longitudinal self-field, backward wave, or exact velocity field in Puffin.',
                    'Puffin global BF is not the earlier fixed-central-ID local RMS BF.'])
    save(out/'RESULT.json',result)
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['program','energy_J']);writer.writeheader();writer.writerows(summary)
    np.savez_compressed(out/'spectra.npz',frequency_Hz=frequency,**{name:sd for name,sd in spectra})
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    colors={'Puffin':'#32936f','PIC':'#1769a5','GENESIS':'#c46b10'}
    for name,sd in spectra:ax[0].plot(100*(frequency/f0-1),sd*f0*1e10,label=name,color=colors[name])
    ax[0].set(xlabel='Offset from comparison carrier (%)',ylabel='Band energy density (pJ / %)',title='Same w64-a sample, 7 periods')
    ax[0].legend();ax[0].grid(alpha=.2)
    names=['Puffin','PIC','GENESIS\noriginal','GENESIS\nhalf-bucket']
    bar=ax[1].bar(names,[x['energy_J']*1e12 for x in summary],color=['#32936f','#1769a5','#c46b10','#8c63a9'])
    ax[1].bar_label(bar,fmt='%.2f',padding=3);ax[1].set(ylabel='Rect ±5% band energy (pJ)',title='Same ±560 µm aperture')
    ax[1].margins(y=.16)
    fig.supxlabel('Puffin: matched external magnet, unaveraged forward paraxial field, no space charge. No amplitude matching.',fontsize=9)
    fig.savefig(out/'comparison.png',dpi=180);plt.close(fig)


if __name__=='__main__':main(Path(sys.argv[1]))
