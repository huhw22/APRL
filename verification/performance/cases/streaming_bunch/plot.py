#!/usr/bin/env python3
"""Export scientific comparison figures from review.py's checked report."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=json.loads(args.review.read_text()); cases=report['cases']; out=args.output
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25,
                         'figure.dpi':120,'savefig.dpi':180})
    colors={'clamp':'#1565c0','mur':'#d84315','mur-gauss':'#9c27b0','zero':'#777777','full':'#222222'}
    fig,axes=plt.subplots(2,3,figsize=(13,7.4),layout='constrained')
    for row,beta in enumerate((0.,.23)):
        for mode,label in [('clamp','Reference clamp'),('mur','Mur residual (unstable at beta=0.23)')]:
            selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and
                             c['metrics']['mode']==mode],key=lambda c:c['metrics']['guard_cells'])
            x=[c['metrics']['guard_cells'] for c in selected]
            ys=[[max(abs(c['expansion_increment_relative_difference'])*100,1e-12) for c in selected],
                [c['particle_difference']['delta_energy_lab_eV']['rms_abs'] for c in selected],
                [c['field_difference']['E']['relative_l2']*100 for c in selected]]
            for col,y in enumerate(ys):
                axes[row,col].plot(x,y,'o-',color=colors[mode],label=label)
        for col in range(3):
            axes[row,col].set_yscale('log');axes[row,col].set_xticks([16,24,40])
            axes[row,col].set_xlabel('Guard per side (cells; dz = 1 um)')
            axes[row,col].set_title(f'Box beta = {beta:g}')
        axes[row,0].set_ylabel('Error in radial expansion increment (%)')
        axes[row,1].set_ylabel('Per-particle lab energy difference, RMS (eV)')
        axes[row,2].set_ylabel('Whole-plane native E difference, L2 (%)')
    axes[0,0].legend(fontsize=8,loc='best')
    axes[0,0].text(.03,.03,'Display floor: 1e-12%',transform=axes[0,0].transAxes,fontsize=8)
    fig.suptitle('Charged-bunch moving-window prototype vs complete domain\n'
                 'Same particles and prescribed release; field comparison includes all recorded times',fontsize=13)
    fig.savefig(out/'guard_comparison.png');fig.savefig(out/'guard_comparison.pdf');plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(11,7.2),layout='constrained')
    for row,beta in enumerate((0.,.23)):
        base=next(c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']=='full')
        for mode in ('full','clamp','mur','mur-gauss'):
            case=next(c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']==mode and c['metrics']['guard_cells']==16)
            hist=np.genfromtxt(Path(case['path'])/'history.csv',names=True,delimiter=',')
            axes[row,0].semilogy(hist['step'],np.maximum(hist['total_gauss_after_shift'],1e-16),
                                 color=colors[mode],label=mode)
            field=np.genfromtxt(Path(case['path'])/'field-plane.csv',names=True,delimiter=',')
            steps=np.unique(field['step']); magnitude=[]
            for step in steps:
                at=field[field['step']==step]
                magnitude.append(np.sqrt(np.mean(sum(at[k]**2 for k in ('Ex_V_m','Ey_V_m','Ez_V_m')))))
            axes[row,1].semilogy(steps,magnitude,color=colors[mode],label=mode)
        for col in range(2):
            axes[row,col].set_xlabel('PIC step'); axes[row,col].set_title(f'Box beta = {beta:g}; guard = 16 cells')
        axes[row,0].set_ylabel('Gauss residual / fixed initial source L2')
        axes[row,1].set_ylabel('Native E magnitude, plane RMS (V/m)')
    axes[0,0].legend(fontsize=9,ncol=2)
    fig.suptitle('A small Gauss residual does not establish stable physical fields',fontsize=13)
    fig.savefig(out/'boundary_stability.png');fig.savefig(out/'boundary_stability.pdf');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for i,beta in enumerate((0.,.23)):
        selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']=='clamp'],
                        key=lambda c:c['metrics']['guard_cells'])
        labels=[f"W={c['metrics']['window_nz']}" for c in selected]
        bottom=np.zeros(len(selected))
        for key,label in [('initialization_wall_s','Initial full Cowan solve'),('cache_preparation_wall_s','Reference cache'),
                          ('evolution_wall_s','Window evolution + all diagnostics')]:
            values=np.array([c['metrics'][key] for c in selected])
            axes[i].bar(labels,values,bottom=bottom,label=label);bottom+=values
        axes[i].set_title(f'Box beta = {beta:g}');axes[i].set_ylabel('One probe run, wall time (s)')
    axes[0].legend(fontsize=8)
    fig.suptitle('Measured prototype cost: one MPI rank / one physical core\n'
                 'Full reference cache remains; no production or asymptotic speedup claim',fontsize=12)
    fig.savefig(out/'prototype_cost.png');fig.savefig(out/'prototype_cost.pdf');plt.close(fig)


if __name__=='__main__':
    main()
