"""Derive full-source moving-window coverage using the accepted r2 mapping."""
import copy
import math
from pathlib import Path
import h5py
import numpy as np
import yaml
from common import C,read,save,need,sha


def geometry(a,plan):
    G=plan['boost_gamma'];bb=math.sqrt(1-1/G**2)
    dz=G*(1+bb)*plan['wavelength_m']/20;dt=dz/C
    u=a['proper_velocity'];gamma=np.sqrt(1+(u*u).sum(axis=1));beta=u[:,2]/gamma
    w=a['macro_weight'].astype(np.longdouble);t=a['arrival_time_offset_s'].astype(np.longdouble)
    mb=(w*beta).sum()/w.sum();mbt=(w*beta*t).sum()/w.sum()
    light=(np.hypot(1.,u[:,2])+u[:,2])/(G+math.sqrt(G-1)*math.sqrt(G+1))
    hu=.5*(light-1/light);held=hu/np.hypot(1.,hu)
    center=float(C*(mb*t.max()-mbt)+2*dz/G)
    ref=float((np.longdouble(center)/C+mbt)/mb)
    snap=beta*C*(ref-a['arrival_time_offset_s']);inp=-snap/(C*beta)
    origin=float(np.min(inp-2*dz/(G*C*(bb+held))))
    z=-G*C*(bb+held)*(inp-origin);translation=-.5*float(z.min()+z.max())
    last=plan['last_plane_m']
    # 128 extra optical cycles are observation padding, not extra bunch charge.
    # They also cover the old six-segment late particle/field-cap cutoff.
    span=float(t.max()-t.min())+(6*77+128)*plan['wavelength_m']/C+(1-bb)*last/C+4e-15
    steps=math.ceil(span*G/dt);active=math.ceil(last/(G*dz))+64+128+2
    extent=max((1+held.max())*(steps+2)+66-translation/dz,translation/dz-64+active+2)
    full=math.ceil(2*extent/256)*256
    stop=G*G*(center+bb*C*origin)+math.sqrt(G*G-1)*dz*steps
    need(held.min()>0,'Held motion must be forward')
    return dict(center=center,reference_time_s=ref,origin_time_s=origin,translation_m=translation,
        dz=dz,dt=dt,steps=steps,active_nz=active,reference_nz=full,stop=stop,
        input_span_s=float(t.max()-t.min()),extra_observation_padding_cycles=128,
        grid_entry_fraction=float((translation/dz+full/2)%1))


def prepare(root):
    plan=read(root/'PLAN.json');d=root/'six-aprl';parent=Path(plan['aprl_parent'])
    old=read(parent/'PLAN.json');puf=read(root/'six/PLAN.json');info=read(root/'six/INPUT_VALIDATED.json')
    with h5py.File(info['input']) as f:a=f['particles/records'][:]
    geo=geometry(a,dict(boost_gamma=old['physics']['boost_gamma'],wavelength_m=puf['wavelength_m'],last_plane_m=puf['outputs'][-1]['actual_lab_z_m']))
    base=yaml.safe_load((parent/'energy-binomial3/run.yaml').read_text())
    need(base['sources']['magnetic_elements']==puf['magnetic_lattice'],'APRL/Puffin six-segment magnet lists differ')
    nx=56;dx=base['mesh']['cell_size'][0]
    need(float(abs(a['plane_position_m']).max())+2*dx<(nx/2-6)*dx,'Gaussian tails need more transverse space; do not clip')
    planes=[['p'+str(x['segment']),x['actual_lab_z_m']] for x in puf['outputs']]
    entries=[]
    for name,shape in [('energy-binomial3','linear-binomial3'),('energy-momentum','quintic-compensated')]:
        case=d/name;case.mkdir()
        for sub in ['detectors','ledger','lab-audit','worldtube']:(case/sub).mkdir()
        card=copy.deepcopy(base)
        card['output']['manifest']=str(case/'run-manifest.yaml')
        card['mesh'].update(cells=[nx,nx,geo['reference_nz']],duration=(geo['steps']+.5)*geo['dt'],longitudinal_shape=shape)
        card['mesh']['moving_window'].update(cells=geo['active_nz'],reference_directory=str(case/'reference'),particle_collar_cells=8)
        card['beam']['reference']['initial_center_z']=geo['center']
        card['beam']['input'].update(file=info['input'],electrons=info['physical_electrons'])
        card['stop']=dict(mode='reference-center-z',z=geo['stop'])
        all_planes=[['entry',0.],['und-entry',.01],['audit-entry',.03],*planes]
        card['detectors'].update(directory=str(case/'detectors'),substep_particle_planes=True,
            field_planes=[dict(name=n,z=z,rhythm=geo['dt']/old['physics']['boost_gamma'],buffer_samples=8,compression=1,particle_background={'enabled':False}) for n,z in planes],
            particle_planes=[dict(name=n+'-particles',z=z,buffer_records=32768,compression=0) for n,z in all_planes],
            energy_planes=[dict(name=n+'-transport',z=z,rhythm=geo['dt']/old['physics']['boost_gamma'],buffer_samples=64,compression=0) for n,z in all_planes])
        card['energy_ledger'].update(directory=str(case/'ledger'),sample_interval_steps=1,discrete_operator_audit=True)
        card['runtime']['resource_monitor']['progress_interval_steps']=100
        (case/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
        entries.append(dict(name=name,shape=shape,directory=str(case),expected_steps=geo['steps']))
    preflight=[]
    for name,method,observer in [('binomial3-off','energy-binomial3',False),('binomial3-on','energy-binomial3',True),('momentum-on','energy-momentum',True)]:
        case=d/'preflight'/name;case.mkdir(parents=True)
        for sub in ['detectors','ledger','lab-audit','worldtube']:(case/sub).mkdir()
        card=yaml.safe_load((d/method/'run.yaml').read_text())
        card['output']['manifest']=str(case/'run-manifest.yaml')
        card['mesh']['moving_window']['reference_directory']=str(case/'reference')
        card['detectors']['directory']=str(case/'detectors');card['energy_ledger']['directory']=str(case/'ledger')
        card['mesh']['duration']=(40-.25)*geo['dt'];card['energy_ledger']['discrete_operator_audit']=observer
        (case/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
        preflight.append(dict(name=name,directory=str(case),observer=observer,expected_steps=40))
    out=copy.deepcopy(old)
    out.update(root=str(d),previous_campaign=str(parent),binary=plan['aprl_binary'],
        input_file=info['input'],input_sha256=info['input_sha256'],input_particles=info['records'],
        input_represented_electrons=info['physical_electrons'],cases=entries,preflight=preflight,
        steps=geo['steps'],active_nz=geo['active_nz'],reference_nz=geo['reference_nz'],
        stop_reference_z_m=geo['stop'],physics_timeout_s=70*3600,
        worldtube_exit_planes=planes,worldtube_half_width_m=(nx/2-6-1)*dx,
        source_exclusion='Prescribed transverse current sheets at +/-880um; worldtube faces +/-840um, one cell inward.',
        comparison_aperture_half_width_m=.00052,original_input_unchanged=False,
        independent_random_source=True,target_weight=355,actual_weight=info['actual_weight'],transverse_cells=nx,
        input_shared_with_Puffin=str(root/'six'),magnetic_lattice_identical=True,
        field_and_particle_planes_match_Puffin_actual_steps=True,
        lab_field_raw_upper_bytes_per_case=geo['steps']*geo['active_nz']*6*8,
        original_binary_reused=True,observer_source_changes=[])
    out['physics']['charge_C']=info['total_charge_C']
    save(d/'PLAN.json',out);save(d/'GEOMETRY.json',geo)
    save(d/'INPUT_VALIDATED.json',info)
    return out
