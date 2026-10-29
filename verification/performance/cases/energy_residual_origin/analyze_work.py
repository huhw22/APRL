#!/usr/bin/env python3
"""Decompose measured free-drift work without fitting a correction to a residual."""
import argparse
import csv
import json
import os
from pathlib import Path
import h5py
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-matplotlib-free-work')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD=np.longdouble


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    a=p.parse_args(); root=a.root.resolve()
    plan=json.loads((root/'plan.json').read_text()); case=Path(plan['case'])
    assert json.loads((root/'COMPLETION.json').read_text())['checks_pass']
    with h5py.File(case/'ledger/energy-ledger.h5') as f: r=f['energy_ledger/records'][:]
    assert len(r)==258
    files=sorted((case/'mesh-work').glob('mesh-work-rank-*.csv'))
    assert [x.name for x in files]==[f'mesh-work-rank-{k:05d}.csv' for k in range(16)]
    values={}
    for file in files:
        with file.open() as f:
            reader=csv.DictReader(f); names=reader.fieldnames; rows=list(reader)
        assert len(rows)==257 and [int(x['step']) for x in rows]==list(range(1,258))
        for n in names[1:]:
            data=np.array([LD(x[n]) for x in rows],dtype=LD)
            assert np.isfinite(data).all()
            if n not in values: values[n]=np.zeros(257,dtype=LD)
            values[n]+=data
    for prefix in ('native','center'):
        assert np.allclose(values[prefix+'_mid_J'],.5*(values[prefix+'_old_J']+values[prefix+'_new_J']),rtol=1e-15,atol=1e-35)
    cumulative={n:np.r_[LD(0),np.cumsum(v,dtype=LD)] for n,v in values.items()}
    def ledger(name): return np.array(r[name],dtype=LD)
    grid=ledger('particle_grid_work_J')
    constraints=ledger('frozen_constraint_work_J')+ledger('carrier_constraint_work_J')
    source=grid-constraints
    field=ledger('field_exchange_residual_J')
    deltaU=ledger('delta_field_energy_J')
    flux=ledger('outward_field_energy_J')
    remap=ledger('moving_remap_field_change_J')
    inlet=ledger('inlet_overwrite_work_J')
    boundary=ledger('boundary_source_work_J')
    release=int(np.flatnonzero(r['released_charge_C']==r['released_charge_C'][-1])[0])
    first_exit=int(np.flatnonzero(r['removed_macroparticles']>0)[0])
    all_exit=int(np.flatnonzero(r['active_macroparticles']==0)[0])
    stages=[('injection',0,release),('released_no_exit',release,first_exit-1),
            ('particle_exit',first_exit-1,all_exit),('after_all_particle_exit',all_exit,257),('whole',0,257)]
    gap={}; propagation={}; summaries={}
    for name,v in cumulative.items():
        if 'pairing' in name: continue
        gap[name]=source-v
        propagation[name]=deltaU+flux-remap-inlet-boundary+v
        diff=np.max(np.abs(field-gap[name]-propagation[name]))
        assert diff<1e-30, (name,float(diff))
        summaries[name]={}
        for label,i,j in stages:
            sums={'particle_work':grid,'constraints':constraints,'particle_minus_constraints':source,
                  'measured_mesh_work':v,'coupling_gap':gap[name],
                  'propagation_and_boundary_gap':propagation[name],'field_exchange_residual':field}
            item={n+'_J':float(x[j]-x[i]) for n,x in sums.items()}
            per_step_gap=np.diff(gap[name][i:j+1])
            per_step_particle=np.diff(grid[i:j+1])
            per_step_source=np.diff(source[i:j+1])
            item['sum_abs_step_coupling_gap_J']=float(np.sum(abs(per_step_gap)))
            item['sum_abs_step_source_work_J']=float(np.sum(abs(per_step_source)))
            item['sum_abs_step_coupling_gap_over_abs_step_source_work']=float(np.sum(abs(per_step_gap))/np.sum(abs(per_step_source))) if np.any(per_step_source) else None
            item['maximum_abs_step_coupling_gap_J']=float(np.max(abs(per_step_gap)))
            item['maximum_step_gap_over_particle_work']=float(np.max(abs(per_step_gap/per_step_particle))) if np.all(per_step_particle!=0) else None
            item.update(start_step=i,end_step=j)
            denom=abs(grid[j]-grid[i])
            item['coupling_gap_over_abs_particle_work']=float((gap[name][j]-gap[name][i])/denom) if denom else None
            summaries[name][label]=item
    output=root/'analysis';output.mkdir(exist_ok=True)
    result={'checks_pass':True,'frame':'boost','unit':'J',
            'all_released_step':release,'first_exit_step':first_exit,'all_exit_step':all_exit,
            'summaries':summaries,'measured_inlet_endpoint_pairing_correction_J':
              {n:float(v[-1]) for n,v in cumulative.items() if 'pairing' in n},
            'definitions':{
              'coupling_gap':'W_Boris_grid - W_frozen - W_carrier - independently measured mesh work',
              'propagation_and_boundary_gap':'delta_U + Q_six - W_remap - W_inlet - W_boundary + measured mesh work',
              'native':'Volume quadrature of native edge J*E; shared edges have adjacent-cell weights',
              'center':'Volume quadrature of separately collocated J and E',
              'mid':'E endpoint average after MPI electric interface completion, BEFORE inlet overwrite',
              'clean_interval':'All input particles released; no particle has yet exited; constrained-source work is constant'},
            'limits':['Mesh work is measured independently, never inferred from a residual.',
                      'Neither norm is claimed to be the exact Cowan discrete invariant.',
                      'Whole-run frozen/carrier work uses existing midpoint path estimates; use the released/no-exit interval to isolate feedback.',
                      'One small case; no production error extrapolation.']}
    (output/'WORK_DECOMPOSITION.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (output/'work_values.csv').open('w') as f:
        writer=csv.writer(f); names=list(cumulative)
        writer.writerow(['step','particle_grid_work_J','constraints_J','field_exchange_residual_J']+names)
        for i in range(258): writer.writerow([i,str(grid[i]),str(constraints[i]),str(field[i])]+[str(cumulative[n][i]) for n in names])
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True,layout='constrained')
    for n,label,color in [('native_old_J','Native J · old E','#377eb8'),('native_mid_J','Native J · mean E','#e68613'),('center_mid_J','Centred J · mean E','#984ea3')]:
        axes[0].plot(r['step'],cumulative[n]*1e18,label=label,color=color)
        axes[1].plot(r['step'],gap[n]*1e18,label=label,color=color)
    axes[0].plot(r['step'],source*1e18,label='Boris work − constrained-source work',color='black',linestyle='--')
    axes[1].plot(r['step'],field*1e18,label='Full field-accounting residual',color='#555555',linestyle='--')
    for ax in axes:
        ax.axvspan(release,first_exit-1,color='#dddddd',alpha=.3)
        for x in (release,first_exit,all_exit): ax.axvline(x,color='#aaaaaa',linewidth=.7)
        ax.axhline(0,color='#cccccc',linewidth=.7);ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[0].set_ylabel('Cumulative work (aJ)');axes[0].set_title('Free drift: independently measured particle and mesh work')
    axes[0].text(.99,.97,'Grey interval: all particles released, no exits (steps 61–180)',transform=axes[0].transAxes,ha='right',va='top',fontsize=9)
    axes[1].set_ylabel('Particle minus mesh work (aJ)');axes[1].set_xlabel('Maxwell step')
    fig.savefig(output/'free_drift_work.png',dpi=180);fig.savefig(output/'free_drift_work.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True,layout='constrained')
    for n,label,color in [('native_mid_J','Native-edge work','#e68613'),('center_mid_J','Cell-centred work','#984ea3')]:
        x=np.arange(release+1,first_exit)
        dg=np.diff(gap[n][release:first_exit])
        dw=np.diff(grid[release:first_exit])
        axes[0].plot(x,dg*1e21,label=label,color=color)
        axes[1].plot(x,100*dg/abs(dw),label=label,color=color)
    for ax in axes: ax.axhline(0,color='#999999',linewidth=.7);ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[0].set_title('Free drift: every step in the released/no-exit interval')
    axes[0].set_ylabel('Particle − mesh work per step (zJ)')
    axes[1].set_ylabel('Gap / |particle work per step| (%)');axes[1].set_xlabel('Maxwell step')
    fig.savefig(output/'clean_interval_steps.png',dpi=180);fig.savefig(output/'clean_interval_steps.pdf');plt.close(fig)
    print(json.dumps({'released_no_exit':summaries['native_mid_J']['released_no_exit'],
                      'whole':summaries['native_mid_J']['whole']},indent=2))


if __name__=='__main__':main()
