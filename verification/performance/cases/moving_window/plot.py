#!/usr/bin/env python3
"""Scientific figures for the unchanged-kernel moving-window checks."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('review',type=Path)
    p.add_argument('--controls',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=json.loads(args.review.read_text());args.output.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.grid':True,'axes.axisbelow':True,'grid.alpha':.25,'savefig.dpi':180})
    widths=[c['width'] for c in r['costs']];x=np.arange(len(widths))
    fig,ax=plt.subplots(1,2,figsize=(10.8,4.2),layout='constrained')
    for dx,mode,label,color in [(-.18,'fixed','Original update, fixed window','#555555'),(.18,'compact','Same update + integer movement','#1565c0')]:
        med=np.array([c['modes'][mode]['total_s']/64*1000 for c in r['costs']])
        lo=np.array([c['modes'][mode]['total_s_min']/64*1000 for c in r['costs']])
        hi=np.array([c['modes'][mode]['total_s_max']/64*1000 for c in r['costs']])
        ax[0].bar(x+dx,med,width=.34,label=label,color=color,yerr=np.array([med-lo,hi-med]),capsize=3)
    for i,c in enumerate(r['costs']):
        ax[0].annotate(f"+{100*c['same_width_total_overhead_fraction']:.2f}%",(i+.18,c['modes']['compact']['total_s']/64*1000),
                       xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
    bottom=np.zeros(3)
    for key,label,color in [('field_move_s','Native E/B/J movement','#1565c0'),('history_move_s','CPML history movement / zeroing','#e69f00'),
                          ('incoming_s','New-layer provider + fill','#009e73')]:
        y=np.array([c['modes']['compact'][key]/64*1000 for c in r['costs']])
        ax[1].bar(x,y,bottom=bottom,label=label,color=color,width=.55);bottom+=y
    for a in ax:a.set_xticks(x,widths);a.set_xlabel('Window z cells');a.legend(fontsize=8)
    ax[0].set_ylabel('Time per step (ms)');ax[1].set_ylabel('Movement cost per step (ms)')
    fig.suptitle('Unchanged Cowan + CPML kernels; 40 x 40 transverse cells\n'
                 'One MPI rank / one physical core; medians of 3 runs, error bars show range',fontsize=12)
    fig.savefig(args.output/'movement_cost.png');fig.savefig(args.output/'movement_cost.pdf');plt.close(fig)
    if args.controls:
        fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
        for row,(moving,fixed,depth) in enumerate([('wave','wave-static',8),('wave20','wave-static20',20)]):
            a=np.genfromtxt(args.controls/moving/'wave-plane.csv',delimiter=',',names=True)
            b=np.genfromtxt(args.controls/fixed/'wave-plane.csv',delimiter=',',names=True)
            n=a['step'];norm=np.max(a['combined_reference_rms_V_m'])
            ax[row,0].plot(n,a['Ey_reference_center_V_m'],color='#333333',label='Large original domain')
            ax[row,0].plot(n,a['Ey_candidate_center_V_m'],'--',color='#1565c0',label='Moving CPML window')
            ax[row,0].set_ylabel('Center native Ey (V/m)');ax[row,0].set_title(f'Pulse trace; z CPML = {depth} cells')
            ax[row,1].semilogy(n,np.maximum(a['combined_difference_rms_V_m']/norm,1e-16),color='#1565c0',label='Moving minus large')
            ax[row,1].semilogy(n,np.maximum(b['combined_difference_rms_V_m']/norm,1e-16),color='#e69f00',label='Fixed nearby CPML minus large')
            ax[row,1].axvspan(100,176,color='#999999',alpha=.1,label='Predefined late interval')
            ax[row,1].set_ylabel('Plane difference RMS / main pulse peak')
            ax[row,1].text(.02,.03,'Zero values shown at 1e-16',transform=ax[row,1].transAxes,fontsize=8)
            ax[row,1].set_title(f'Common large-domain reference; z CPML = {depth}')
            for col in range(2):ax[row,col].set_xlabel('Step');ax[row,col].legend(fontsize=8)
        fig.suptitle('Independent vacuum-wave comparison with original CPML formulas\n'
                     'Fixed and moving end locations differ; error norms are not isolated reflection coefficients',fontsize=12)
        fig.savefig(args.output/'cpml_wave_comparison.png');fig.savefig(args.output/'cpml_wave_comparison.pdf');plt.close(fig)


if __name__=='__main__':
    main()
