"""Fixed, absolutely normalized LW acceleration mode; not a fitted correction."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import c,e,hbar
from scipy.special import roots_hermitenorm
from theory import arc_spectrum
from prepare import save

def main():
    ws=Path('/public/home/wfollowm/APRL_Workspace')
    root=ws/'results/radiation_method_validation_20261103'
    base=root/'synchrotron-r1'
    plan=json.loads((base/'energy/PLAN.json').read_text())
    data=np.load(base/'energy/analysis/field_spectra.npz')
    en=data['energy_eV'];band=(en>=8)&(en<=12);freq=en[band]
    xx,yy=np.meshgrid(data['x_m'],data['y_m'])
    nodes,w=roots_hermitenorm(8);w=w/np.sqrt(2*np.pi)
    sx,sy=np.meshgrid(nodes*plan['sigma_x_m'],nodes*plan['sigma_y_m'])
    weights=np.outer(w,w).ravel()
    model=np.empty((len(freq),xx.size,3),complex)
    for start in range(0,xx.size,16):
        stop=min(xx.size,start+16)
        points=np.stack(np.broadcast_arrays(
          xx.ravel()[start:stop,None]-sx.ravel(),
          yy.ravel()[start:stop,None]-sy.ravel(),
          np.full((stop-start,64),plan['parameters']['screen'])),axis=-1).reshape(-1,3)
        E,_=arc_spectrum(points,freq,plan['parameters'],512)
        model[:,start:stop]=np.einsum('fpqd,q->fpd',E.reshape(len(freq),stop-start,64,3),weights)
    form=np.exp(-.5*(freq*e/hbar*plan['sigma_z_m']/(plan['parameters']['beta']*c))**2)
    model=model.reshape(len(freq),len(data['y_m']),len(data['x_m']),3)*form[:,None,None,None]*plan['represented_electrons']
    window=np.hanning(len(data['y_m']))[:,None]*np.hanning(len(data['x_m']))[None,:]
    model=model[...,:2]
    denominator=np.sum(abs(model*window[None,:,:,None])**2,axis=(1,2,3))
    result={}
    for name in ['boris','energy','momentum']:
        a=np.load(base/name/'analysis/field_spectra.npz')
        F=a['measured_Fourier_fields'][band,:,:,:2]
        alpha=np.sum(np.conj(model)*F*window[None,:,:,None]**2,axis=(1,2,3))/denominator
        result[name]=dict(fixed_reference_mode_band_power_ratio=float(np.sum(abs(alpha)**2*denominator)/sum(denominator)),
          complex_amplitude_RMS_error=float(np.sqrt(np.sum(abs(alpha-1)**2*denominator)/sum(denominator))),
          amplitude_real=alpha.real.tolist(),amplitude_imag=alpha.imag.tolist())
    save(root/'GAUSSIAN_MODE_AUDIT.json',dict(energy_eV=freq.tolist(),methods=result,
      scope='Projection onto a fixed absolute theoretical mode; no theory rescaling and no claim of general near-field separation',
      reference='Acceleration only; full LW velocity/interference correction separately quantified in primary aperture results'))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
