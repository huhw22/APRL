#!/usr/bin/env python3
"""Standalone scientific figures from the reviewed, inspectable CSV only."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-space-charge-review-mpl')
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(out):
    with (out/'comparison_all_planes.csv').open() as f:
        rows=list(csv.DictReader(f))
    names={'uniform':'Cold uniform sphere','converging':'Converging uniform sphere','gaussian':'Cold spherical Gaussian'}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'axes.titleweight':'semibold','axes.labelsize':10,'legend.fontsize':9})
    fig,axes=plt.subplots(4,3,figsize=(13.8,12.2),sharex=True,sharey='row',layout='constrained')
    for col,(kind,title) in enumerate(names.items()):
        data=[r for r in rows if r['case']==kind]
        x=np.array([float(r['mean_rest_tau']) for r in data])
        for row,(pic,ref,scale,ylabel) in enumerate([
            ('PIC_sigma_x_m','reference_sigma_x_m',1e6,'Transverse rms size (um)'),
            ('PIC_mean_rest_gain_eV','reference_mean_rest_gain_eV',1,'Mean rest-frame gain (eV/electron)'),
            ('PIC_std_rest_gain_eV','reference_std_rest_gain_eV',1,'Std of rest-frame gains (eV/electron)'),
            ('PIC_mean_lab_gain_eV','reference_mean_lab_gain_eV',1e-3,'Mean lab-frame gain (keV/electron)')]):
            ax=axes[row,col]
            ax.plot(x,[float(r[ref])*scale for r in data],color='#bc7c29',ls='--',marker='x',lw=1.6,ms=6,label='Spherical theory')
            ax.plot(x,[float(r[pic])*scale for r in data],color='#356f9b',marker='o',lw=1.6,ms=4.5,label='APRL')
            ax.grid(alpha=.18)
            if row in (1,3):ax.axhline(0,color='#555555',lw=.6)
            if col==0:ax.set_ylabel(ylabel)
            if row==3:ax.set_xlabel('Mean rest-frame crossing time / t0')
            if row==0:ax.set_title(title,pad=10)
    axes[0,0].legend(loc='upper left')
    fig.suptitle('Space-charge comparison at fixed laboratory planes',fontsize=15,fontweight='semibold')
    fig.supxlabel('Same particle IDs and actual event times; gain relative to each entry event. Lines connect sampled planes.',fontsize=10)
    for ext in ('png','pdf'):fig.savefig(out/f'space_charge_comparison.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13.8,4.6),sharex=True,sharey=True,layout='constrained')
    for ax,(kind,title) in zip(axes,names.items()):
        data=[r for r in rows if r['case']==kind]
        x=[float(r['mean_rest_tau']) for r in data]
        for field,color,style,mark,label in [
            ('mean_lab_gain_error_eV','#333333','-','o','Measured lab gain error'),
            ('rest_longitudinal_momentum_error_contribution_lab_eV','#356f9b','--','x','From longitudinal momentum'),
            ('rest_energy_error_contribution_lab_eV','#bc7c29',':','s','From rest-frame energy')]:
            ax.plot(x,[float(r[field])*1e-3 for r in data],color=color,ls=style,marker=mark,lw=1.6,ms=4,label=label)
        ax.axhline(0,color='#555555',lw=.6);ax.grid(alpha=.18);ax.set_title(title)
        ax.set_xlabel('Mean rest-frame crossing time / t0')
    axes[0].set_ylabel('APRL minus theory (keV/electron)')
    axes[0].legend(loc='lower left')
    fig.suptitle('Lorentz decomposition of the laboratory mean-energy discrepancy',fontsize=14,fontweight='semibold')
    fig.supxlabel('Exact particle energy transformation: error = Gamma * rest-energy error + Gamma * beta * mc^2 * longitudinal-momentum error',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(out/f'lab_energy_error_decomposition.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    render(p.parse_args().output)
