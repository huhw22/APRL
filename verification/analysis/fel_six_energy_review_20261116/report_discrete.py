#!/usr/bin/env python3
"""Combine a completed production/small-run discrete ledger; no control run needed."""
import argparse,csv,json,os
from pathlib import Path
import h5py
import numpy as np
from review_discrete import validate_columns,require,save,LD
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-discrete-ledger-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(directory,physical):
    with h5py.File(physical) as f:
        group=f['energy_ledger'];rows=group['records'][:]
        require(int(group.attrs['format_version'])==5,'Need unchanged v5 physical ledger')
        require(int(group['complete'][()])==1 and int(group['committed_records'][()])==len(rows),'Physical ledger incomplete')
        ranks=int(group.attrs['mpi_size']);frame=group.attrs['frame']
    require(b'boosted computational frame' in frame if isinstance(frame,bytes) else 'boosted computational frame' in frame,'Wrong frame')
    files=sorted(directory.glob('discrete-field-rank-*.csv'))
    require([p.name for p in files]==[f'discrete-field-rank-{k:05d}.csv' for k in range(ranks)],'Missing/extra rank outputs')
    steps=int(rows['step'][-1]);sums={}
    for path in files:
        require(path.with_suffix('.complete').read_text().split()==['1',str(steps)],'Rank output incomplete')
        with path.open() as f:reader=csv.DictReader(f);records=list(reader);names=reader.fieldnames
        require([int(x['step']) for x in records]==list(range(1,steps+1)),'Missing/repeated step')
        v={n:np.array([LD(x[n]) for x in records]) for n in names if n!='step'}
        require(all(np.isfinite(x).all() for x in v.values()),'Nonfinite field ledger')
        validate_columns(v)
        for n,x in v.items():
            if n not in sums:sums[n]=np.zeros_like(x)
            sums[n]=np.maximum(sums[n],x) if n.startswith('max_') or '_prediction_scale_' in n else sums[n]+x
    return rows,sums,ranks,steps

def calculate(rows,sums,steps):
    index=np.asarray(rows['step'],dtype=np.int64)
    require(index[0]==0 and index[-1]==steps and np.all(np.diff(index)>0),'Invalid sample cadence')
    a=lambda n:np.asarray(rows[n],dtype=LD)
    full=lambda n:np.r_[LD(0),np.cumsum(sums[n],dtype=LD)]
    c=lambda n:full(n)[index]
    native=np.r_[LD(0),np.cumsum(sums['native_end_J']-sums['native_begin_J'],dtype=LD)][index]
    field=c('field_operator_residual_J');mid=c('applied_current_mid_J')
    time=c('time_electric_J')+c('time_magnetic_J')+c('current_time_pairing_J')
    cowan=c('cowan_operator_J');cpml=c('magnetic_cpml_J')+c('electric_cpml_J')
    total=a('delta_particle_kinetic_box_J')+native+c('discrete_curl_exchange_J')-time-cowan-cpml-c('inlet_overwrite_J')-c('remap_J')-a('prescribed_boris_work_J')-a('boundary_source_work_J')-a('frozen_constraint_work_J')-a('carrier_constraint_work_J')
    coupling=a('particle_grid_work_J')-a('boundary_source_work_J')-a('frozen_constraint_work_J')-a('carrier_constraint_work_J')-mid
    spatial=(a('delta_field_energy_J')-a('inlet_overwrite_work_J')-a('moving_remap_field_change_J'))-(native-c('inlet_overwrite_J')-c('remap_J'))
    flux=a('outward_field_energy_J')-c('discrete_curl_exchange_J')
    old=a('balance_residual_boris_J')
    scale=abs(native)+abs(c('discrete_curl_exchange_J'))+abs(time)+abs(cowan)+abs(cpml)+abs(c('inlet_overwrite_J'))+abs(c('remap_J'))+abs(a('delta_particle_kinetic_box_J'))+abs(a('prescribed_boris_work_J'))+abs(a('frozen_constraint_work_J'))+abs(a('carrier_constraint_work_J'))+abs(a('boundary_source_work_J'))+abs(old)+abs(total)
    bound=512*LD(np.finfo(float).eps)*np.maximum(scale,LD('1e-100'))
    split_error=abs(total-field-coupling-a('particle_work_balance_residual_J')-a('magnetic_rotation_delta_K_J'))
    bridge_error=abs(old-total-spatial-flux-time-cowan-cpml)
    require(np.all(split_error<=bound),'Particle/field split differs')
    require(np.all(bridge_error<=bound),'Physical/discrete bridge differs')
    quantities={'old_monitor_residual_J':old,'discrete_total_residual_J':total,'field_operator_residual_J':field,
      'particle_grid_coupling_gap_J':coupling,'spatial_definition_J':spatial,'flux_definition_J':flux,
      'paired_time_term_J':time,'Cowan_operator_term_J':cowan,'physical_cpml_term_J':cpml,
      'unapplied_current_mid_work_J':c('unapplied_current_mid_J')}
    denominator=abs(a('delta_particle_kinetic_box_J')[-1])
    summary={'final':{n:float(x[-1]) for n,x in quantities.items()},
      'max_split_identity_error_J':float(np.max(split_error)),'max_bridge_identity_error_J':float(np.max(bridge_error)),
      'max_serialization_roundoff_bound_J':float(np.max(bound)),
      'max_abs_cumulative_field_operator_residual_J':float(np.max(abs(full('field_operator_residual_J')))),
      'final_discrete_total_over_abs_net_particle_change':float(total[-1]/denominator) if denominator else None,
      'old_final_over_abs_net_particle_change':float(old[-1]/denominator) if denominator else None}
    return index,quantities,summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ledger-directory',type=Path,required=True)
    p.add_argument('--physical-ledger',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();directory=a.ledger_directory.resolve();out=a.output.resolve()
    require(out!=directory,'Use a separate analysis directory');out.mkdir(parents=True,exist_ok=True)
    require(not list(out.iterdir()),'Output directory must be empty')
    physical=(a.physical_ledger or directory/'energy-ledger.h5').resolve()
    rows,sums,ranks,steps=load(directory,physical);index,q,summary=calculate(rows,sums,steps)
    save(out/'summary.json',{'checks_pass':True,'frame':'boost','unit':'J','mpi_ranks':ranks,'field_steps':steps,
      'physical_ledger':str(physical),'physical_sample_count':len(rows),'scope':'Discrete numerical accounting; time/Cowan/curl-exchange terms are not radiated energy or physical external work.',
      'physical_conservation_accepted':False,**summary})
    with (out/'discrete_energy.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['step']+list(q));w.writerows([[int(step)]+[str(v[k]) for v in q.values()] for k,step in enumerate(index)])
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True,layout='constrained')
    for n,label,color in [('old_monitor_residual_J','Original collocated monitor','#555555'),('discrete_total_residual_J','Discrete accounting + particle gap','#e68613'),('field_operator_residual_J','Field operator residual','#984ea3')]:
        axes[0].plot(index,q[n]*1e18,label=label,color=color)
    for n,label,color in [('spatial_definition_J','Spatial definition','#377eb8'),('flux_definition_J','Flux definition','#4daf4a'),('paired_time_term_J','Time/current pairing','#ff7f00'),('Cowan_operator_term_J','Cowan operator','#984ea3')]:
        axes[1].plot(index,q[n]*1e18,label=label,color=color)
    axes[0].set_title('Boost-frame energy ledger: physical monitor and discrete accounting')
    axes[1].set_title('Explicit numerical / definition contributions')
    for ax in axes:ax.axhline(0,color='#aaaaaa',linewidth=.7);ax.grid(alpha=.2);ax.set_ylabel('Cumulative energy (aJ)');ax.legend(fontsize=9)
    axes[1].set_xlabel('Maxwell step');fig.savefig(out/'discrete_energy.png',dpi=180);fig.savefig(out/'discrete_energy.pdf');plt.close(fig)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
