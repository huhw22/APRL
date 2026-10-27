#!/usr/bin/env python3
"""Scientific figures from the checked residual-CPML comparison."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def save(fig, output, name):
    fig.savefig(output/(name+'.png'))
    fig.savefig(output/(name+'.pdf'))
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=json.loads(args.review.read_text());cases=report['cases']
    require=report['cohort_charge_and_input_consistency_status']=='PASS'
    if not require:raise RuntimeError('Use a reviewed cohort/charge report')
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    betas=sorted({c['metrics']['beta_box'] for c in cases})
    plt.rcParams.update({'font.size':10,'axes.grid':True,'axes.axisbelow':True,
                         'grid.alpha':.23,'savefig.dpi':180})
    guard_colors={16:'#56b4e9',24:'#1565c0',40:'#009e73'}
    fig,axes=plt.subplots(len(betas),4,figsize=(15,3.8*len(betas)),squeeze=False,layout='constrained')
    for row,beta in enumerate(betas):
        selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']=='residual'],
                        key=lambda c:c['metrics']['guard_cells'])
        x=[c['metrics']['guard_cells'] for c in selected]
        values=[[100*c['expansion_increment_relative_difference'] for c in selected],
                [c['particle_difference']['delta_energy_lab_eV']['rms_abs'] for c in selected],
                [c['moment_difference']['mean_delta_energy_lab_eV'] for c in selected],
                [100*c['sigma_energy_relative_difference'] for c in selected]]
        labels=['Radial expansion increment difference (%)','Per-particle lab energy difference RMS (eV)',
                'Mean particle delta-E difference (eV)','Lab energy spread difference (%)']
        for col,(y,label) in enumerate(zip(values,labels)):
            axes[row,col].plot(x,y,'o-',color='#1565c0')
            axes[row,col].axhline(0,color='#444444',linewidth=.8)
            axes[row,col].set_xticks([16,24,40]);axes[row,col].set_xlabel('Nonabsorbing guard per side (cells)')
            axes[row,col].set_ylabel(label);axes[row,col].set_title(f'Particle box beta = {beta:g}')
    fig.suptitle('Residual-input CPML window versus full-domain control: fixed IDs and prescribed release\n'
                 'Longitudinal CPML adds 20 cells outside each guard; original transverse PEC and update kernels',fontsize=12)
    save(fig,out,'particle_guard_comparison')

    fig,axes=plt.subplots(len(betas),3,figsize=(12,3.7*len(betas)),squeeze=False,layout='constrained')
    for row,beta in enumerate(betas):
        selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']=='residual'],
                        key=lambda c:c['metrics']['guard_cells'])
        x=[c['metrics']['guard_cells'] for c in selected]
        for col,vector in enumerate(('E','B')):
            axes[row,col].plot(x,[100*c['field_difference'][vector]['relative_l2'] for c in selected],
                               'o-',color='#1565c0')
            axes[row,col].set_ylabel(f'Whole-plane native {vector} relative L2 difference (%)')
        for wavelength,color in [('1e-06','#e69f00'),('5e-06','#009e73'),('1e-05','#cc79a7')]:
            axes[row,2].plot(x,[c['arrival_current_bunching_amplitude_difference'][wavelength] for c in selected],
                               'o-',color=color,label=f'lambda = {float(wavelength)*1e6:g} um')
        axes[row,2].set_ylabel('Arrival-current bunching amplitude difference')
        axes[row,2].legend(fontsize=8)
        for ax in axes[row]:
            ax.set_xticks([16,24,40]);ax.set_xlabel('Nonabsorbing guard per side (cells)')
            ax.set_title(f'Particle box beta = {beta:g}')
    fig.suptitle('Native plane-field differences and fixed-plane arrival-current bunching\n'
                 'No simultaneous lab E/B or pure-radiation interpretation',fontsize=12)
    save(fig,out,'field_and_bunching_guard_comparison')

    fig,axes=plt.subplots(len(betas),2,figsize=(11,3.7*len(betas)),squeeze=False,layout='constrained')
    for row,beta in enumerate(betas):
        selected=[c for c in cases if c['metrics']['beta_box']==beta]
        for case in selected:
            mode=case['metrics']['mode'];guard=case['metrics']['guard_cells']
            label='Full domain' if mode=='full' else ('Total CPML, g24' if mode=='cpml' else f'Residual CPML, g{guard}')
            color='#222222' if mode=='full' else ('#d55e00' if mode=='cpml' else guard_colors[guard])
            for col,vector in enumerate(('E','cB')):
                frames=case['field_by_time']
                values=[max(f[vector]['candidate_rms_V_m'],1e-16) for f in frames]
                axes[row,col].semilogy([f['step'] for f in frames],values,color=color,
                                       linestyle='--' if mode=='cpml' else '-',label=label)
        for col,vector in enumerate(('E','cB')):
            axes[row,col].set_ylabel(f'Native {vector}, plane RMS (V/m)')
            axes[row,col].set_xlabel('PIC step');axes[row,col].set_title(f'Particle box beta = {beta:g}')
            axes[row,col].legend(fontsize=8)
    fig.suptitle('Plane fields during independent full-domain and moving-window evolution\n'
                 'Native staggered field times; zero values displayed at 1e-16 V/m',fontsize=12)
    save(fig,out,'native_field_history')

    fig,axes=plt.subplots(len(betas),3,figsize=(12,3.6*len(betas)),squeeze=False,layout='constrained')
    for row,beta in enumerate(betas):
        selected=[c for mode in ('full','cpml','residual') for c in cases if
                  c['metrics']['beta_box']==beta and c['metrics']['mode']==mode and c['metrics']['guard_cells']==24]
        labels=['Full' if c['metrics']['mode']=='full' else ('Residual g24' if c['metrics']['mode']=='residual' else 'Total g24') for c in selected]
        for col,(key,label) in enumerate([('radius_growth_percent','Radial expansion (%)'),
                                          ('mean_delta_energy_lab_eV','Mean lab delta E (eV; positive = gain)'),
                                          ('sigma_energy_lab_eV','Lab energy spread (eV)')]):
            axes[row,col].bar(labels,[c['moments'][key] for c in selected],
                               color=['#333333' if c['metrics']['mode']=='full' else '#1565c0' if c['metrics']['mode']=='residual' else '#d55e00' for c in selected])
            axes[row,col].set_ylabel(label);axes[row,col].set_title(f'Particle box beta = {beta:g}')
    fig.suptitle('Same guard24 geometry: total versus residual CPML derivative input\n'
                 'Original total-field particle advance and supplied initial reference remain unchanged',fontsize=12)
    save(fig,out,'cpml_input_particle_comparison')

    fig,axes=plt.subplots(1,len(betas),figsize=(5.7*len(betas),4.5),squeeze=False,layout='constrained')
    for col,beta in enumerate(betas):
        selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode'] in ('full','cpml','residual')],
                        key=lambda c:({'full':0,'cpml':1,'residual':2}[c['metrics']['mode']],c['metrics']['guard_cells']))
        labels=['Full\nW=424' if c['metrics']['mode']=='full' else f"{'Total' if c['metrics']['mode']=='cpml' else 'Residual'} g{c['metrics']['guard_cells']}\nW={c['metrics']['window_nz']}" for c in selected]
        bottom=np.zeros(len(selected))
        for key,label in [('initialization_wall_s','Initial Cowan solve'),('cache_preparation_wall_s','Frozen reference cache'),
                          ('evolution_wall_s','Evolution including diagnostics')]:
            values=np.array([c['metrics'][key] for c in selected])
            axes[0,col].bar(labels,values,bottom=bottom,label=label);bottom+=values
        axes[0,col].set_ylabel('Wall time for one probe run (s)');axes[0,col].set_title(f'Particle box beta = {beta:g}')
        axes[0,col].legend(fontsize=8)
    fig.suptitle('Different-domain workload observations; 1 MPI rank / 1 physical core\n'
                 'One run per case, full initial cache retained, research diagnostics included; not production speedup',fontsize=12)
    save(fig,out,'prototype_workload_cost')

    fig,axes=plt.subplots(len(betas),3,figsize=(12.5,3.7*len(betas)),squeeze=False,layout='constrained')
    for row,beta in enumerate(betas):
        comparisons=sorted([c for c in report['residual_comparisons'] if c['beta_box']==beta],
                           key=lambda c:c['guard_cells'])
        guards=[c['guard_cells'] for c in comparisons]
        selections=[next(c for c in cases if c['path']==comparison['residual_path']) for comparison in comparisons]
        series=[('New residual CPML','#1565c0','-',selections)]
        for label,color,style in [('previous total-field CPML','#d55e00','--'),
                                  ('historical reference clamp','#8e44ad',':')]:
            series.append((label,color,style,[next(c for c in comparison['secondary_historical_controls']
                                                    if c['label']==label) for comparison in comparisons]))
        for label,color,style,entries in series:
            current=label=='New residual CPML'
            e=[100*c['field_difference' if current else 'field_difference_vs_own_full']['E']['relative_l2'] for c in entries]
            mean=[c['moment_difference' if current else 'moment_difference_vs_own_full']['mean_delta_energy_lab_eV'] for c in entries]
            spread=[100*c['sigma_energy_relative_difference'] for c in entries]
            for col,y in enumerate((e,mean,spread)):
                axes[row,col].plot(guards,y,marker='o',linestyle=style,color=color,label=label)
        for col,label in enumerate(('Native E relative L2 difference (%)','Mean lab delta-E difference (eV)',
                                     'Energy spread relative difference (%)')):
            axes[row,col].set_ylabel(label);axes[row,col].set_xlabel('Nonabsorbing guard per side (cells)')
            axes[row,col].set_xticks([16,24,40]);axes[row,col].set_title(f'Particle box beta = {beta:g}')
            axes[row,col].legend(fontsize=7)
    fig.suptitle('Primary residual comparison plus explicitly historical boundary-method references\n'
                 'Each uses its own full control; historical clamp lacks the additional CPML layers',fontsize=12)
    save(fig,out,'historical_boundary_comparison')

    fig,axes=plt.subplots(1,2,figsize=(10.8,4.2),layout='constrained')
    for beta,color in zip(betas,('#1565c0','#d55e00')):
        selected=sorted([c for c in cases if c['metrics']['beta_box']==beta and c['metrics']['mode']=='residual'],
                        key=lambda c:c['metrics']['guard_cells'])
        guards=[c['metrics']['guard_cells'] for c in selected]
        axes[0].plot(guards,[c['metrics']['residual_reference_work_bytes']/1024 for c in selected],
                     'o-',color=color,label=f'Particle beta = {beta:g}')
        axes[1].plot(guards,[100*c['residual_work_accounting']['input_wall_fraction_of_field_advance'] for c in selected],
                     'o-',color=color,label=f'Particle beta = {beta:g}')
    axes[0].set_ylabel('Reference-input buffers (KiB; within grid total)')
    axes[1].set_ylabel('Reference-input time / field-advance time (%)')
    for ax in axes:
        ax.set_xlabel('Nonabsorbing guard per side (cells)');ax.set_xticks([16,24,40]);ax.legend(fontsize=8)
    fig.suptitle('Residual-reference work is a component of existing memory and timing totals\n'
                 'Do not add these quantities again; one run per case, no repeated-throughput claim',fontsize=12)
    save(fig,out,'residual_reference_work')


if __name__=='__main__':
    main()
