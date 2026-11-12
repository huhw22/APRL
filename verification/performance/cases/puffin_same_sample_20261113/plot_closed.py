"""Redraw the closed small caches using stable program colors in both panels."""
import sys
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import read,save,W

root=Path(sys.argv[1]);result=read(root/'analysis/RESULT.json')
if not result['complete']:raise RuntimeError('Analysis is incomplete')
with np.load(root/'analysis/spectra.npz') as f:
    freq=f['frequency_Hz'];f0=freq[len(freq)//2]
    spectra={name:f[name] for name in ['Puffin','PIC','GENESIS']}
with (W/'results/sparc_random_startup_20261112_r2/reviewed-r1/pic_fixed_cohorts.csv').open() as stream:
    row=next(r for r in csv.DictReader(stream) if r['case']=='w64-a' and r['plane']=='u07'
             and r['scope']=='all' and r['group_cells']=='1')
save(root/'analysis/ADDITIONAL_CLOSED_DIAGNOSTICS.json',dict(
    PIC_full_pulse_global_b1=float(row['global_b1']),
    Puffin_full_pulse_global_b1=result['Puffin_full_pulse_global_b1'],
    BF_scope='All sampled particles, full-pulse GLOBAL BF; not central fixed-ID local RMS.',
    Puffin_over_PIC_band_energy=result['comparisons'][0]['energy_J']/result['comparisons'][1]['energy_J'],
    Puffin_over_PIC_carrier_spectral_density=float(spectra['Puffin'][len(freq)//2]/spectra['PIC'][len(freq)//2]),
    spectral_Pearson_correlation=float(np.corrcoef(spectra['Puffin'],spectra['PIC'])[0,1]),
    no_amplitude_fit=True))
colors=['#32936f','#1769a5','#c46b10','#8c63a9']
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for (name,sd),color in zip(spectra.items(),colors):
    axes[0].plot(100*(freq/f0-1),sd*f0*1e10,label=name,color=color)
axes[0].set(xlabel='Offset from comparison carrier (%)',ylabel='Band energy density (pJ / %)',title='Same w64-a sample, 7 periods')
axes[0].legend();axes[0].grid(alpha=.2)
bar=axes[1].bar(['Puffin','PIC','GENESIS\noriginal','GENESIS\nhalf-bucket'],
                [r['energy_J']*1e12 for r in result['comparisons']],color=colors)
axes[1].bar_label(bar,fmt='%.2f',padding=3)
axes[1].set(ylabel='Rect ±5% band energy (pJ)',title='Same ±560 µm aperture');axes[1].margins(y=.16)
fig.supxlabel('Puffin: matched external magnet, unaveraged forward paraxial field, no space charge. No amplitude matching.',fontsize=9)
fig.savefig(root/'analysis/comparison.png',dpi=180)
