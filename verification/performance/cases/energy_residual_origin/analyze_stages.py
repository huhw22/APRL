#!/usr/bin/env python3
"""Locate field-accounting changes by observing actual, unchanged operations."""
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
STAGES=['begin','magnetic_curl','magnetic_cpml','electric_curl_local',
        'electric_cpml','inlet_overwrite','electric_mpi','remap']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve()
    plan=json.loads((root/'plan.json').read_text());case=Path(plan['case'])
    assert plan['field_stages'] and json.loads((root/'COMPLETION.json').read_text())['checks_pass']
    with h5py.File(case/'ledger/energy-ledger.h5') as f:r=f['energy_ledger/records'][:]
    files=sorted((case/'field-stages').glob('field-stages-rank-*.csv'))
    assert [p.name for p in files]==[f'field-stages-rank-{i:05d}.csv' for i in range(16)]
    summed={};maxima={'upper_tangential_max_delta_V_per_m','max_abs_J_any_A_per_m2'}
    for path in files:
        with path.open() as f: reader=csv.DictReader(f); rows=list(reader)
        assert len(rows)==257*8
        assert [(int(x['step']),x['stage']) for x in rows]==[(i,s) for i in range(1,258) for s in STAGES]
        for name in reader.fieldnames:
            if name in ('step','stage'):continue
            values=np.array([LD(x[name]) for x in rows],dtype=LD).reshape(257,8)
            assert np.isfinite(values).all()
            if name not in summed:summed[name]=np.zeros((257,8),dtype=LD)
            if name in maxima:summed[name]=np.maximum(summed[name],values)
            elif name=='current_sampled':summed[name]=np.maximum(summed[name],values)
            else:summed[name]+=values
    def col(name,stage):return summed[name][:,STAGES.index(stage)]
    def ledger(name):return np.array(r[name],dtype=LD)
    u=summed['center_energy_J'];du=np.diff(u,axis=1)
    old=ledger('field_energy_interior_J')
    scale=np.maximum(abs(old[:-1])+abs(old[1:]),LD('1e-40'))
    bound=256*np.finfo(float).eps*scale
    checks={}
    def check(name,error,tol):
        checks[name]={'pass':bool(np.all(abs(error)<=tol)),'maximum_absolute_error':float(np.max(abs(error))),
                      'maximum_roundoff_bound':float(np.max(tol))}
        assert checks[name]['pass'],(name,checks[name])
    check('initial_stage_matches_previous_ledger_inventory',u[:,0]-old[:-1],bound)
    check('remap_stage_matches_next_ledger_inventory',u[:,-1]-old[1:],bound)
    dinlet=du[:,4];dremap=du[:,6]
    check('inlet_actual_energy_change_matches_ledger',dinlet-np.diff(ledger('inlet_overwrite_work_J')),bound)
    check('remap_actual_energy_change_matches_ledger',dremap-np.diff(ledger('moving_remap_field_change_J')),bound)
    power_names=[f'power_{axis}{side}_W' for axis in 'xyz' for side in ('minus','plus')]
    power=np.stack([summed[n] for n in power_names],axis=2)
    dt=LD(str(r['time_box_s'][1]-r['time_box_s'][0]))
    q=.5*dt*(power[:,0,:]+power[:,6,:]-(power[:,5,:]-power[:,4,:]))
    check('stage_flux_matches_six_face_ledger',q-np.diff(ledger('outward_field_energy_face_J'),axis=0),bound[:,None])
    mesh=np.zeros(257,dtype=LD)
    for path in sorted((case/'mesh-work').glob('mesh-work-rank-*.csv')):
        with path.open() as f:rows=list(csv.DictReader(f))
        mesh+=np.array([LD(x['center_mid_J']) for x in rows])
    source=ledger('particle_grid_work_J')-ledger('frozen_constraint_work_J')-ledger('carrier_constraint_work_J')
    coupling=np.diff(source)-mesh
    physical_cpml=du[:,1]+du[:,3]
    curl_and_mpi=du[:,0]+du[:,2]+du[:,5]
    propagation=curl_and_mpi+physical_cpml+np.sum(q,axis=1)+mesh
    reconstructed=propagation+coupling
    check('all_observed_operations_reproduce_residual_step_changes',
          reconstructed-np.diff(ledger('field_exchange_residual_J')),bound)
    current=col('max_abs_J_any_A_per_m2','electric_mpi')
    begin_particles=col('particles','begin');begin_carriers=col('carriers','begin')
    zero=(current==0)&(begin_particles==0)&(begin_carriers==0)
    zero_start=next((i for i in range(257) if np.all(zero[i:])),None)
    phases=[('injection',0,61),('released_no_exit',61,180),('particle_exit',180,239),('after_all_physical_particles_exit',239,257)]
    if zero_start is not None:phases.append(('no_particles_no_carriers_exact_zero_current',zero_start,257))
    summaries={}
    for label,i,j in phases:
        summaries[label]={'start_ledger_step':i,'end_ledger_step':j,
          'magnetic_curl_change_J':float(np.sum(du[i:j,0])),
          'magnetic_cpml_direct_change_J':float(np.sum(du[i:j,1])),
          'electric_curl_local_change_J':float(np.sum(du[i:j,2])),
          'electric_cpml_direct_change_J':float(np.sum(du[i:j,3])),
          'inlet_overwrite_direct_change_J':float(np.sum(dinlet[i:j])),
          'electric_mpi_completion_change_J':float(np.sum(du[i:j,5])),
          'moving_remap_direct_change_J':float(np.sum(dremap[i:j])),
          'six_face_outward_J':float(np.sum(q[i:j])),
          'center_mesh_work_J':float(np.sum(mesh[i:j])),
          'particle_mesh_coupling_gap_J':float(np.sum(coupling[i:j])),
          'propagation_and_flux_gap_J':float(np.sum(propagation[i:j])),
          'ledger_residual_change_J':float(r['field_exchange_residual_J'][j]-r['field_exchange_residual_J'][i]),
          'begin_carriers_first':int(begin_carriers[i]) if j>i else None,
          'begin_carriers_last':int(begin_carriers[j-1]) if j>i else None,
          'maximum_current_A_per_m2':float(np.max(current[i:j])) if j>i else None}
    result={'checks_pass':all(x['pass'] for x in checks.values()),'checks':checks,
      'stages':STAGES,'source':str(case),'frame':'boost','energy_unit':'J',
      'first_complete_step_with_no_particles_no_carriers_exact_zero_current':zero_start+1 if zero_start is not None else None,
      'no_current_suffix_starting_ledger_step':zero_start,
      'maximum_physical_cpml_direct_energy_change_per_step_J':float(np.max(abs(physical_cpml))),
      'sum_abs_physical_cpml_direct_energy_change_J':float(np.sum(abs(du[:,1]))+np.sum(abs(du[:,3]))),
      'maximum_upper_tangential_field_change_before_remap_V_per_m':float(np.max(summed['upper_tangential_max_delta_V_per_m'][:,:-1])),
      'phases':summaries,
      'interpretation_limits':['Energy changes of operations are not individually errors; curl/current/MPI terms must be summed.',
        'CPML result covers the physical ledger aperture; CPML can still affect incoming waves indirectly.',
        'A residual in propagation+flux does not distinguish numerical field evolution from an approximate energy/flux diagnostic.',
        'The last tangential E plane being fixed during a field step is an observed boundary rule, not by itself proof of contamination inside the trusted collar.']}
    output=root/'stage-analysis';output.mkdir(exist_ok=True)
    (output/'FIELD_STAGE_REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (output/'operation_values.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['step']+STAGES[1:]+['outward_J','center_mesh_work_J','coupling_gap_J','propagation_gap_J','max_J_A_per_m2','begin_particles','begin_carriers'])
        for i in range(257):w.writerow([i+1]+[str(x) for x in du[i]]+[str(np.sum(q[i])),str(mesh[i]),str(coupling[i]),str(propagation[i]),str(current[i]),int(begin_particles[i]),int(begin_carriers[i])])
    fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained')
    steps=np.arange(1,258)
    ax[0].plot(steps,np.cumsum(propagation)*1e18,label='Curl/MPI + physical CPML + flux + mesh work',color='#555555')
    ax[0].plot(steps,np.cumsum(coupling)*1e18,label='Particle–mesh work gap',color='#e68613')
    ax[0].set_ylabel('Cumulative gap (aJ)');ax[0].legend(fontsize=8)
    ax[0].set_title('Free drift: observed field operations and remaining energy gap')
    ax[1].plot(steps,begin_particles,label='Physical particles at start of step',color='#377eb8')
    ax[1].plot(steps,begin_carriers,label='Current carriers at start of step',color='#984ea3')
    ax[1].set_ylabel('Count');ax[1].set_xlabel('Maxwell step');ax[1].legend(fontsize=8)
    for x in ax:
        x.axvline(239,color='#999999',linestyle=':');x.grid(alpha=.2)
        if zero_start is not None:x.axvspan(zero_start+1,257,color='#cde8d0',alpha=.45)
    if zero_start is not None:ax[0].set_title('Free drift: observed field operations and remaining energy gap\n'+f'Green: exactly zero sources, steps {zero_start+1}–257',fontsize=11)
    fig.savefig(output/'field_stage_audit.png',dpi=180);fig.savefig(output/'field_stage_audit.pdf');plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
