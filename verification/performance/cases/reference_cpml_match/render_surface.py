#!/usr/bin/env python3
"""Compare fixed-plane shell theory and common boost-time momentum separately."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(p.read_text())
def save(fig,p):
    fig.savefig(p.with_suffix('.png'),dpi=180,bbox_inches='tight')
    fig.savefig(p.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig)
def main(root):
    parent=root.parent;out=root/'figures';out.mkdir(exist_ok=True)
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    cases=[('Original truncated reference',parent/'longitudinal_momentum_origin_20261030_r3/review/base','#67778d'),
           ('Extend reference into CPML (failed)',parent/'reference_cpml_match_20261030/outer_500/momentum-analysis','#bf6962'),
           ('Prescribed surface current',root/'surface_500/momentum-analysis','#16886d')]
    for label,d,c in cases:
        x=np.genfromtxt(d/'common_boost_time.csv',names=True,delimiter=',')
        ax[0].plot(x['step'],1e3*x['mean_rest_uz_error'],label=label,color=c)
        ax[1].plot(x['step'],x['lab_energy_error_from_same_box_slice_eV']/1e3,label=label,color=c)
    for a in ax:a.axhline(0,color='k',lw=.6);a.grid(alpha=.2);a.set_xlabel('Field step')
    ax[0].set_ylabel('Mean rest-frame $u_z$ error ($10^{-3}$)')
    ax[1].set_ylabel('Lab mean-energy error on same boost slice (keV)')
    ax[0].legend(fontsize=8);fig.suptitle('Complete released cohort, before any particle exits: uniform sphere')
    save(fig,out/'pilot_boundary_comparison')
    if not all((root/(k+'_full')/'COMPLETION.json').exists() for k in ('uniform','converging','gaussian')):return
    fig,ax=plt.subplots(3,3,figsize=(13,10));summary={}
    metrics=[('mean_lab_gain_eV',1e-3,'Mean lab gain (keV)'),('mean_rest_gain_eV',1.,'Mean rest-frame gain (eV)'),('sigma_x_m',1e6,'Transverse RMS size ($\mu$m)')]
    for row,k in enumerate(('uniform','converging','gaussian')):
        new=read(root/(k+'_full')/'analysis/RESULT.json');oldroot=Path(read(root/'plan.json')['cases'][k+'_full']['original'])
        old=read(oldroot/'analysis/RESULT.json');summary[k]={'original':old['planes'][-1],'surface':new['planes'][-1],'ledger':read(root/(k+'_full')/'LEDGER_CHECK.json')}
        for col,(key,scale,ylabel) in enumerate(metrics):
            a=ax[row,col];x=[v['lab_z_m'] for v in new['planes']]
            a.plot(x,[v['reference_'+key]*scale for v in new['planes']],color='k',label='Independent shell theory')
            a.plot(x,[v['PIC_'+key]*scale for v in old['planes']],'s--',color='#67778d',label='Original reference boundary',ms=4)
            a.plot(x,[v['PIC_'+key]*scale for v in new['planes']],'o-',color='#16886d',label='Prescribed surface current',ms=4)
            a.set_ylabel(ylabel);a.set_title(k);a.grid(alpha=.2)
            if row==2:a.set_xlabel('Fixed lab detector z (m)')
    ax[0,0].legend(fontsize=8);fig.suptitle('Same input events, mesh, Cowan/Boris/deposition, and original CPML')
    fig.tight_layout();save(fig,out/'full_space_charge_comparison')
    (root/'PHYSICS_COMPARISON.json').write_text(json.dumps(summary,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root)
