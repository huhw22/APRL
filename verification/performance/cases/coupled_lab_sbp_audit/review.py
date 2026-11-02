"""Plot and summarize completed scalar budgets; never read running HDF5."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    p=argparse.ArgumentParser();p.add_argument('result',type=Path);a=p.parse_args();r=a.result.resolve();result=json.loads((r/'RESULT.json').read_text());summaries={}
    fig,axes=plt.subplots(2,2,figsize=(12,7.5),layout='constrained')
    for col,name in enumerate(['undulator-on','drift-on']):
        ax=axes[0,col];detail=axes[1,col]
        d=r/name;x=np.genfromtxt(d/'whole_steps_J.csv',delimiter=',',names=True,dtype=np.longdouble);old=np.genfromtxt(d/'lab_balance_paired.csv',delimiter=',',names=True,dtype=np.longdouble);Q=np.longdouble(json.loads((d/'RESULT.json').read_text())['charge_abs_C']);D=np.max(abs(old['particle_delta_lab_eV']))*Q;assert D>0;mask=x['clean']==1
        values={'Particle kinetic change':x['particle_change'],'Cowan operator + Gauss':x['Cowan_energy_operator']+x['Cowan_momentum_operator']+x['Gauss_B_momentum'],'Force/current coupling difference':x['energy_coupling_gap']+x['electric_impulse_gap']+x['magnetic_impulse_gap'],'Total numerical discrepancy':x['total_numerical_discrepancy']};columns=[x['step'][mask]]
        for color,(label,v) in enumerate(values.items()):
            y=100*np.cumsum(v[mask])/D;steps=np.r_[x['step'][mask][0]-1,x['step'][mask]];curve=np.r_[0,y];ax.plot(steps,curve,label=label,color='C'+str(color));columns.append(y)
            if color:detail.plot(steps,curve,label=label,color='C'+str(color))
        ax.set(title=name,xlabel='Field step (all particles active)',ylabel='% of fixed full-record particle-change scale');ax.grid(alpha=.2);ax.axhline(0,color='black',lw=.5)
        detail.set(title='Numerical terms (same denominator)',xlabel='Field step (all particles active)',ylabel='% of same fixed particle-change scale');detail.grid(alpha=.2);detail.axhline(0,color='black',lw=.5)
        np.savetxt(d/'clean_relative.csv',np.column_stack(columns),delimiter=',',header='step,particle_percent,operator_percent,coupling_percent,total_percent',comments='')
        v=result['cases'][name]['clean_interval']['eV_per_electron'];src=sum(v[k+'_source'] for k in ['magnetic_supply','magnetic_cpml','magnetic_sync','electric_hold','electric_cpml','electric_inlet','electric_sync','reference_current']);peak=float(np.max(abs(np.cumsum(x['total_numerical_discrepancy'][mask])))/Q)
        summaries[name]=dict(normalization=dict(denominator_J=float(D),denominator_eV_per_electron=float(D/Q),definition='max absolute full-record stable-cohort lab particle change including actual physical removals; shared by all curves; clean interval cumulative curves reset at first clean step'),clean_numerical_max_abs_eV=peak,clean_final_numerical_eV=v['total_numerical_discrepancy'],ratio_to_abs_net_particle_change=abs(v['total_numerical_discrepancy']/v['particle_change']),relative_to_full_record_particle_scale_percent=float(100*abs(v['total_numerical_discrepancy'])*Q/D),Cowan_operator_plus_Gauss_eV=v['Cowan_energy_operator']+v['Cowan_momentum_operator']+v['Gauss_B_momentum'],force_current_coupling_eV=v['energy_coupling_gap']+v['electric_impulse_gap']+v['magnetic_impulse_gap'],budget_eV=dict(particle_change=v['particle_change'],field_reference_change=v['field_reference_change'],net_outward_including_remap=v['native_propagation_out']-v['remap_exchange'],measured_field_sources=src,numerical_discrepancy=v['total_numerical_discrepancy']))
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside upper center',ncol=2,fontsize=9)
    fig.savefig(r/'clean_relative.png',dpi=180);fig.savefig(r/'clean_relative.pdf');plt.close(fig);(r/'REVIEW.json').write_text(json.dumps(summaries,indent=2)+'\n');print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
