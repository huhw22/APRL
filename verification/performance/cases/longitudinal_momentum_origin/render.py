#!/usr/bin/env python3
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-longitudinal-origin-mpl')
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(path):
    with path.open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]

def render(root,existing):
    def analysis(name):
        reviewed=root/'review'/name
        return reviewed if reviewed.exists() else root/name/'analysis'
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.7),layout='constrained')
    for kind,color,label in [('uniform','#356f9b','Cold uniform'),('converging','#ba7d31','Converging uniform'),('gaussian','#a55c81','Cold Gaussian')]:
        data=load(existing/kind/'common_boost_time.csv')
        x=[r['step'] for r in data]
        axes[0].plot(x,[r['mean_rest_uz_error']*1e4 for r in data],color=color,label=label)
        axes[1].plot(x,[r['lab_energy_error_from_same_box_slice_eV']*1e-3 for r in data],color=color,label=label)
    for ax in axes:ax.axhline(0,color='#777777',lw=.7);ax.grid(alpha=.2);ax.set_xlabel('Maxwell step');ax.legend(fontsize=9)
    axes[0].set_ylabel('Mean rest-frame longitudinal momentum error / (1e-4 mc)')
    axes[1].set_ylabel('Lab energy error of same boost-time cohort (keV/electron)')
    fig.suptitle('Existing runs: all particles released, before the first particle exit')
    for ext in ('png','pdf'):fig.savefig(root/f'common_time_existing.{ext}',dpi=180)
    plt.close(fig)
    if not all((analysis(n)/'ANALYSIS.json').exists() for n in ('base','inlet_guard_256','transverse_240')):return
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True,layout='constrained')
    summary={}
    cases=[('base','#356f9b','-','Original 160 x 160, inlet guard 64'),('inlet_guard_256','#ba7d31','--','Same transverse grid, inlet guard 256'),('transverse_240','#a55c81','-.','240 x 240, inlet guard 64')]
    if (analysis('matched_domain_thick_cpml')/'ANALYSIS.json').exists():
        cases.append(('matched_domain_thick_cpml','#78843f',':','Same reference aperture, thicker CPML'))
    for name,color,style,label in cases:
        data=load(analysis(name)/'common_boost_time.csv');force=load(analysis(name)/'measured_force.csv')
        axes[0].plot([r['step'] for r in data],[r['mean_rest_uz_error']*1e4 for r in data],color=color,ls=style,label=label)
        axes[1].plot([r['step'] for r in force],[r['mean_Ez_error_V_per_m']*1e-3 for r in force],color=color,ls=style,label=label)
        summary[name]=json.loads((analysis(name)/'ANALYSIS.json').read_text())
    for ax in axes:ax.axhline(0,color='#777777',lw=.7);ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[0].set_ylabel('Mean rest-frame momentum error / (1e-4 mc)')
    axes[1].set_ylabel('Mean sampled Ez minus reference (kV/m)')
    axes[1].set_xlabel('Maxwell step')
    fig.suptitle('Boundary-distance contrasts; identical input, pusher and field update')
    for ext in ('png','pdf'):fig.savefig(root/f'boundary_contrasts.{ext}',dpi=180)
    plt.close(fig)
    first=root.parent/'longitudinal_momentum_origin_20261030_r4/first_step'
    if (first/'COMPLETE').exists():
        fig,axes=plt.subplots(1,2,figsize=(11.8,4.7),sharey=True,layout='constrained')
        for ax,name,wall in zip(axes,('base','transverse_240'),(350.,550.)):
            with (first/(name+'.csv')).open() as f:rows=list(csv.DictReader(f))
            for stage,color,style,marker,label in [('electric_curl','#356f9b','-',None,'After original electric curl'),('electric_cpml','#ba7d31','--','o','After CPML correction')]:
                r=[x for x in rows if x['stage']==stage and abs(float(x['x_m'])*1e6-wall)<=35.01]
                ax.plot([float(x['x_m'])*1e6-wall for x in r],[(float(x['Ez_V_per_m'])-float(x['translated_reference_Ez_V_per_m']))*1e-3 for x in r],color=color,ls=style,marker=marker,ms=4,mfc='none',label=label)
            ax.axvline(0,color='#777777',lw=.8);ax.axhline(0,color='#777777',lw=.7);ax.grid(alpha=.2)
            ax.set_xlabel('Distance from reference / CPML interface (um)')
            ax.set_title(f'Interface at x = {wall:g} um');ax.legend(fontsize=8)
        axes[0].set_ylabel('One-step Ez minus translated reference (kV/m)')
        fig.suptitle('First field step, no current at this interface: the electric curl creates the pulse')
        for ext in ('png','pdf'):fig.savefig(root/f'first_step_boundary.{ext}',dpi=180)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.6),layout='constrained')
    force=summary['base']['force'];x=[r['step'] for r in force]
    for name,color,label in [('electric_impulse_error','#356f9b','Electric Ez'),('magnetic_impulse_error','#ba7d31','Magnetic rotation'),('net_field_impulse_error','#444444','Net impulse')]:
        axes[0].plot(x,[r[name]*1e7 for r in force],color=color,label=label)
    axes[1].plot(x,[r['mean_Ez_error_V_per_m']*1e-3 for r in force],color='#356f9b')
    for ax in axes:ax.axhline(0,color='#777777',lw=.7);ax.grid(alpha=.2);ax.set_xlabel('Maxwell step')
    axes[0].legend(fontsize=9);axes[0].set_ylabel('Mean longitudinal impulse error / (1e-7 mc)')
    axes[1].set_ylabel('Mean sampled Ez minus ideal reference (kV/m)')
    fig.suptitle('Actual particle-gathered fields and independent Boris replay')
    for ext in ('png','pdf'):fig.savefig(root/f'force_decomposition.{ext}',dpi=180)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('existing',type=Path);a=p.parse_args();render(a.root,a.existing)
