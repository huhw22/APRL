"""Radiation-only handoff: same band, aperture, and observation interval."""
import sys,json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import h,c,e

def main(out,root):
    out=Path(out);root=Path(root)
    pic=json.loads((out/'pic.json').read_text());gen=json.loads((out/'genesis.json').read_text())
    plan=json.loads((root/'plan.json').read_text());photon=h*c/e/plan['physics']['wavelength_m']
    rows=[]
    for p,g in zip(pic,gen):
        assert p['plane']==g['plane']
        T=p['fft_samples']*p['dt_s']
        a=p['results']['rect'];b=g['results']['rect']
        rows.append(dict(plane=p['plane'],z_lab_m=p['z_m'],record_ps=T*1e12,
          PIC_native_band_nJ=a['native_Poynting_band_J']*1e9,
          PIC_E_proxy_band_nJ=a['E_transverse_band_proxy_J']*1e9,
          GENESIS_band_nJ=b['band_J']*1e9,
          PIC_native_record_mean_W=a['native_Poynting_band_J']/T,
          PIC_E_proxy_record_mean_W=a['E_transverse_band_proxy_J']/T,
          GENESIS_record_mean_W=b['band_J']/T,
          PIC_native_over_GENESIS=a['native_Poynting_band_J']/b['band_J'],
          PIC_E_proxy_over_GENESIS=a['E_transverse_band_proxy_J']/b['band_J'],
          PIC_peak_eV=a['peak_energy_eV'],GENESIS_peak_eV=b['peak_energy_eV'],
          PIC_Hann_over_GENESIS=g['results']['hann']['PIC_native_flux_over_GENESIS'],
          GENESIS_time_window_capture=g['captured_pulse_energy_fraction']))
    with (out/'radiation_power.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,
      'legend.fontsize':10,'xtick.labelsize':11,'ytick.labelsize':11,
      'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white',
      'savefig.facecolor':'white','lines.linewidth':1.8})
    colors=['#24648c','#d47928','#668d6b']
    fig,ax=plt.subplots(1,2,figsize=(12.4,4.8),layout='constrained')
    z=np.array([p['z_lab_m'] for p in rows])
    ax[0].semilogy(z,[p['PIC_native_record_mean_W'] for p in rows],'o-',color=colors[0],label='PIC explicit: native E x B')
    ax[0].semilogy(z,[p['PIC_E_proxy_record_mean_W'] for p in rows],':',color=colors[2],label='PIC transverse E-only proxy')
    ax[0].semilogy(z,[p['GENESIS_record_mean_W'] for p in rows],'s-',color=colors[1],label='GENESIS radiation')
    ax[0].set_ylabel('Record-mean band power (W)');ax[0].set_title('(a) Near-resonant power along the lattice')
    ax[0].legend(loc='best')
    ratios=np.array([p['PIC_native_over_GENESIS'] for p in rows])
    ax[1].plot(z,ratios,'o-',color=colors[0],label='PIC / GENESIS')
    ax[1].plot(z,[p['PIC_Hann_over_GENESIS'] for p in rows],'s--',color=colors[2],label='Matched Hann window')
    ax[1].axhspan(.95,1.05,color='#cccccc',alpha=.45,label='5% screening band')
    ax[1].axhline(1,color='#777777',lw=.9);ax[1].set_ylabel('Power ratio');ax[1].set_title('(b) Difference from the reference')
    if ratios.max()/max(ratios.min(),1e-100)>100:ax[1].set_yscale('log')
    ax[1].legend(loc='best')
    for a in ax:
        a.set_xlabel('Fixed laboratory plane z (m)');a.grid(alpha=.20)
    fig.suptitle('SPARC short-bunch SASE: explicit version',fontsize=15)
    fig.savefig(out/'radiation_power.png',dpi=210);plt.close(fig)
    a=np.load(out/'far-spectra.npz');b=np.load(out/'far-genesis.npz')
    T=pic[-1]['fft_samples']*pic[-1]['dt_s']
    fig,axes=plt.subplots(2,2,figsize=(12.4,8.2),layout='constrained')
    for i,(label,title) in enumerate([('rect','Rectangular window'),('hann','Energy-normalized Hann')]):
        xp=a['energy_eV'];xg=b['energy_eV'];ep=a[label+'_E_proxy_J_per_eV']/T
        sp=a[label+'_native_flux_J_per_eV']/T;gp=b[label+'_J_per_eV']/T
        for j in range(2):
            ax=axes[i,j]
            ax.semilogy(xp,np.where(ep>0,ep,np.nan),color=colors[0],label='PIC transverse E-only proxy')
            ax.semilogy(xp,np.where(sp>0,sp,np.nan),'--',color=colors[2],label='PIC native E x B')
            ax.semilogy(xg,np.where(gp>0,gp,np.nan),color=colors[1],label='GENESIS radiation')
            ax.axvline(photon,color='#777777',ls=':',lw=1)
            ax.set_xlim((0,5) if j==0 else (photon*.95,photon*1.05))
            sel=(xp>=ax.get_xlim()[0])&(xp<=ax.get_xlim()[1]);selg=(xg>=ax.get_xlim()[0])&(xg<=ax.get_xlim()[1])
            hi=max(np.nanmax(ep[sel]),np.nanmax(gp[selg]))
            ax.set_ylim(max(hi*(1e-10 if j==0 else 1e-5),1e-20),hi*2)
            ax.set_xlabel('Photon energy (eV)');ax.set_ylabel('Record-mean spectrum (W/eV)')
            ax.set_title(f'({chr(97+2*i+j)}) {title}: '+('wide view' if j==0 else 'resonant band'))
            ax.grid(alpha=.18)
            if i==0 and j==0:ax.legend(loc='upper right')
    fig.suptitle('Outlet at 15.246 m: full physical aperture',fontsize=15)
    fig.savefig(out/'explicit_spectra.png',dpi=210);plt.close(fig)
    result=dict(plane_rows=rows,photon_energy_eV=photon,band_eV=[photon*.95,photon*1.05],
      physical_aperture_mm=[2.4,2.4],power_definition='Band pulse energy divided by the same PIC record duration, not peak power',
      raw_fields=True,near_field_subtraction=False,physical_acceptance=False,
      same_noise_seed=True,genesis_extra_particle_comparison_cancelled_by_user=True)
    (out/'RADIATION_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(rows[-1],indent=2))
if __name__=='__main__':main(sys.argv[1],sys.argv[2])
