"""Compare measured transverse grid-scale content with independent LW fields."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import c,e,h,epsilon_0
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from theory import arc_spectrum
from prepare import save

def main():
    ws=Path('/public/home/wfollowm/APRL_Workspace')
    root=ws/'results/radiation_method_validation_20261103'
    base=root/'single-electron-r1'
    cfg=json.loads((base/'energy/PLAN.json').read_text())
    cache=np.load(base/'energy/analysis/field_spectra.npz')
    x,y=cache['x_m'],cache['y_m'];en=cache['energy_eV']
    band=(en>=50)&(en<=100);freq=en[band];de=en[1]-en[0]
    xx,yy=np.meshgrid(x,y);p=cfg['parameters']
    pts=np.c_[xx.ravel(),yy.ravel(),np.full(xx.size,p['screen'])]
    E,B=arc_spectrum(pts,freq,p,512)
    theory=E.reshape(len(freq),len(y),len(x),3)[...,:2]
    dx=x[1]-x[0];dy=y[1]-y[0]
    kx=np.fft.fftfreq(len(x),dx);ky=np.fft.fftfreq(len(y),dy)
    win=np.hanning(len(y))[:,None]*np.hanning(len(x))[None,:]
    high=(abs(kx)[None,:]>.7/(2*dx))|(abs(ky)[:,None]>.7/(2*dy))
    arrays={'LW acceleration':theory}
    for label in ['boris','energy','momentum']:
        a=np.load(base/label/'analysis/field_spectra.npz')
        arrays[label]=a['measured_Fourier_fields'][band,:,:,:2]
    metrics={};powers={}
    for label,a in arrays.items():
        F=np.fft.fft2(a*win[None,:,:,None],axes=(1,2))
        power=np.sum(abs(F)**2,axis=(0,3));powers[label]=power
        metrics[label]=dict(high_transverse_wavenumber_fraction=float(power[high].sum()/power.sum()))
        denominator=np.sum(abs(theory*win[None,:,:,None])**2,axis=(1,2,3))
        overlap=np.sum(np.conj(theory)*a*win[None,:,:,None]**2,axis=(1,2,3))/denominator
        metrics[label].update(
          fixed_LW_mode_amplitude_real=overlap.real.tolist(),
          fixed_LW_mode_amplitude_imag=overlap.imag.tolist(),
          fixed_LW_mode_band_power_ratio=float(np.sum(abs(overlap)**2*denominator)/np.sum(denominator)),
          mode_amplitude_complex_relative_RMS_error=float(np.sqrt(np.sum(abs(overlap-1)**2*denominator)/np.sum(denominator))))
    factor=2*epsilon_0*c*e/h*de
    measured=factor*np.sum(abs(arrays['energy'])**2,axis=(0,3))
    expected=factor*np.sum(abs(theory)**2,axis=(0,3))
    out=root/'figures';out.mkdir(exist_ok=True)
    with plt.rc_context({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42}):
        fig,axes=plt.subplots(2,2,figsize=(7.5,6.3),layout='constrained')
        extent=[(x[0]-dx/2)*1e6,(x[-1]+dx/2)*1e6,(y[0]-dy/2)*1e6,(y[-1]+dy/2)*1e6]
        norm=LogNorm(vmin=expected.max()*.02,vmax=max(measured.max(),expected.max()))
        for ax,data,title in [(axes[0,0],expected,'LW acceleration reference'),
                              (axes[0,1],measured,'Energy iteration: recorded field')]:
            im=ax.imshow(data,origin='lower',extent=extent,norm=norm,cmap='viridis',aspect='equal')
            ax.set(xlabel='x (um)',ylabel='y (um)',title=title)
        fig.colorbar(im,ax=list(axes[0]),label='50-100 eV transverse E energy proxy (J/m2)',shrink=.88)
        colors={'LW acceleration':'black','boris':'.5','energy':'#1764a5','momentum':'#c76b1d'}
        for label,power in powers.items():
            sx=np.fft.fftshift(power.sum(axis=0));k=np.fft.fftshift(kx)*2*dx
            axes[1,0].semilogy(k,np.maximum(sx/sx.sum(),1e-16),label=label,color=colors[label])
        axes[1,0].set(xlabel='kx / kx,Nyquist',ylabel='Fraction per transverse k bin',
                      ylim=(1e-9,1),xlim=(-1,1),title='Identical spatial Hann window')
        axes[1,0].legend(fontsize=7,frameon=False)
        names=list(arrays);values=[100*metrics[q]['high_transverse_wavenumber_fraction'] for q in names]
        axes[1,1].bar(range(4),values,color=[colors[q] for q in names])
        axes[1,1].set_xticks(range(4),['LW','A','B','C'])
        axes[1,1].set(ylabel='Energy fraction above 0.7 Nyquist (%)',ylim=(0,60),
                      title='Grid-scale content')
        for i,v in enumerate(values):axes[1,1].text(i,v+1,f'{v:.3g}%',ha='center',fontsize=8)
        fig.savefig(out/'single_electron_spatial_audit.png',dpi=450)
        fig.savefig(out/'single_electron_spatial_audit.pdf')
        fig.savefig(out/'single_electron_spatial_audit.svg')
        plt.close(fig)
    np.savez(root/'SPATIAL_REFERENCE.npz',energy_eV=freq,x_m=x,y_m=y,acceleration_E=theory)
    save(root/'SPATIAL_AUDIT.json',dict(metrics=metrics,energy_eV=freq.tolist(),
       definition='Same finite aperture, spatial Hann, high if |kx| or |ky| exceeds 0.7 corresponding Nyquist',
       reference='LW acceleration term; complete LW band velocity/interference correction is separately reported',
       scope='Fixed reference mode projection is diagnostic only: the theory amplitude is never rescaled to match PIC; origin of grid contamination is not uniquely identified'))
    print(json.dumps({label:{k:v for k,v in row.items() if 'amplitude_real' not in k and 'amplitude_imag' not in k} for label,row in metrics.items()},indent=2))

if __name__=='__main__':main()
