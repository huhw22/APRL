"""Render the completed raw-field power review without changing its data."""
import os
from pathlib import Path
import json
import sys

import numpy as np

os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-critical-power-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(root):
    root = Path(root).resolve()
    review = json.loads((root/'RESULT.json').read_text())
    if review.get('complete') is not True:
        raise RuntimeError('Require completed spectral review')
    photon = review['target_photon_eV']
    colors = {'off':'#24648c','on':'#b18428'}
    labels = {'off':'Unfiltered','on':'Seven-point kernel'}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.2,8),layout='constrained')
    for mode in ['off','on']:
        data=np.load(root/mode/'p1-spectra.npz')
        duration=review['cases'][mode]['fft_samples']*review['cases'][mode]['dt_s']
        time=data['retarded_time_s'];time=(time-time[0])*1e15
        axes[0,0].plot(time,data['native_total_power_W']/1e6,color=colors[mode],label=labels[mode])
        energy=data['energy_eV']
        for window,style in [('rect','-'),('hann','--')]:
            density=data[window+'_native_flux_J_per_eV']/duration
            axes[0,1].plot(energy,density,style,color=colors[mode],lw=1.3,label=labels[mode]+' / '+window)
            zoom=(energy>=.9*photon)&(energy<=1.1*photon)
            axes[1,0].plot(energy[zoom],density[zoom],style,color=colors[mode],lw=1.3,label=labels[mode]+' / '+window)
        values=[next(q['native_resonant_band_record_mean_power_W'] for q in review['rows'] if q['mode']==mode and q['window']==window) for window in ['rect','hann']]
        x=np.arange(2)+(-.18 if mode=='off' else .18)
        bars=axes[1,1].bar(x,values,width=.34,color=colors[mode],label=labels[mode])
        axes[1,1].bar_label(bars,fmt='%.1f',padding=3,fontsize=9)
    axes[0,0].set(title='Raw total-field power through the exit plane',xlabel='Time from first sample (fs)',ylabel='Signed power (MW)')
    axes[0,1].set(title='Full spectrum: original stored E and B',xlabel='Photon energy (eV)',ylabel='Signed mean-power density (W/eV)')
    axes[0,1].set_yscale('symlog',linthresh=.01)
    axes[0,1].set_xlim(0,6*photon)
    axes[0,1].axvspan(.95*photon,1.05*photon,color='#777777',alpha=.13)
    axes[1,0].set(title='Resonance near 491.3 nm',xlabel='Photon energy (eV)',ylabel='Signed mean-power density (W/eV)')
    axes[1,0].set_xlim(.9*photon,1.1*photon)
    axes[1,0].axvspan(.95*photon,1.05*photon,color='#777777',alpha=.13)
    axes[1,1].set(title='Integrated resonance band (±5%)',ylabel='Finite-record mean power (W)',xticks=[0,1],xticklabels=['Rectangular','Normalized Hann'])
    axes[1,1].set_ylim(bottom=0)
    for ax in axes.flat:
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle('First-segment field-power control — same 2.4 × 2.4 mm aperture\nRaw total lab fields; stored Yee half-time B; no near-field subtraction',fontsize=12)
    fig.savefig(root/'power_comparison.png',dpi=200)
    fig.savefig(root/'power_comparison.pdf')
    plt.close(fig)


if __name__=='__main__':
    plot(sys.argv[1])
