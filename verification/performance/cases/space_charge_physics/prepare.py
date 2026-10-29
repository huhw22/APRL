#!/usr/bin/env python3
"""Prepare bounded physics comparisons; never build or submit jobs."""
import argparse,copy,hashlib,json,math,shutil
from pathlib import Path
import h5py,numpy as np,yaml
from model import Sphere,G_BOOST,DZ,DT,c

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def require(ok,message):
    if not ok:raise RuntimeError(message)

def geometry(data,last):
    u=np.asarray(data['proper_velocity'],dtype=np.longdouble);arrival=np.asarray(data['arrival_time_offset_s'],dtype=np.longdouble)
    w=np.asarray(data['macro_weight'],dtype=np.longdouble);g=np.sqrt(1+np.sum(u*u,axis=1));bz=u[:,2]/g
    boost=np.longdouble(G_BOOST);gb=np.sqrt((boost-1)*(boost+1));beta=gb/boost
    mean_beta=np.sum(w*bz)/np.sum(w)
    mean_beta_time=np.sum(w*bz*arrival)/np.sum(w)
    required_center=np.longdouble(c)*(mean_beta*np.max(arrival)-mean_beta_time)
    # The common LAB snapshot must not back-propagate any plane event. This
    # selects a coordinate/time origin, not an extra physical propagation age.
    roundoff_margin=256*np.finfo(float).eps*max(1.,abs(float(required_center)))
    initial_center=max(1e-7,float(required_center)+roundoff_margin)
    ref=(np.longdouble(initial_center)/c+mean_beta_time)/mean_beta
    assert ref>=np.max(arrival)
    # Held transverse trajectories omit transverse momentum in their velocity.
    # heldLongitudinalBetaZ first removes LAB transverse momentum, then boosts.
    plus=(np.sqrt(1+u[:,2]**2)+u[:,2])/(boost+gb)
    uz=(plus-1/plus)/2;held=uz/np.sqrt(1+uz*uz)
    assert np.min(held)>=0 and np.max(held)<1
    input_time=arrival-ref
    frame_t0=np.min(input_time-2*DZ/(boost*c*(beta+held)))
    source=-boost*c*(beta+held)*(input_time-frame_t0)
    translation=-.5*(source.min()+source.max())
    initial_ref=boost*boost*(np.longdouble(initial_center)+beta*c*frame_t0)
    active=math.ceil(last/(G_BOOST*DZ))+64+128+2
    # Full reference is streamed, and sized only from its actual query bounds.
    full=256
    for _ in range(100):
        origin=-.5*full*DZ;offset=math.floor((float(translation)-origin)/DZ)-64
        outer=G_BOOST*(origin+(offset+active)*DZ-float(translation))
        stop=outer+6*G_BOOST*DZ
        steps=math.ceil((stop-float(initial_ref))/(float(gb)*DZ));duration=(steps+.5)*DT
        bound=math.ceil(duration/DT)+1
        minimum=offset-(1+float(np.mean(held)))*bound-1
        maximum=offset+active
        need=max((1+float(np.mean(held)))*bound+66-float(translation)/DZ,
                 float(translation)/DZ-64+active+2)
        required=math.ceil(2*need/256)*256
        if full>=required and minimum>=0 and maximum<=full-1:break
        full=max(full+256,required)
    else:raise RuntimeError('Reference sizing did not converge')
    return {'reference_cells':full,'active_cells':active,'initial_offset':offset,
      'duration_s':duration,'expected_stop_step':steps,'stop_reference_lab_z_m':stop,
      'reference_initial_lab_z_m':float(initial_ref),'input_reference_time_offset_s':float(ref),
      'initial_center_lab_z_m':initial_center,'minimum_forward_snapshot_center_z_m':float(required_center),
      'initial_frame_lab_time_s':float(frame_t0),'box_translation_m':float(translation),
      'minimum_reference_query_cells':minimum,'maximum_reference_query_cells':maximum,
      'held_beta_range':[float(held.min()),float(held.max())],
      'origin_reference_m':origin,'outer_lab_z_m':outer}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--cases',nargs='+',choices=('uniform','converging','gaussian'),default=['uniform','converging','gaussian'])
    a=p.parse_args();ws=a.workspace.resolve();root=a.root.resolve();root.mkdir(parents=True,exist_ok=False)
    here=Path(__file__).resolve().parent;binary=ws/'build/discrete-energy-ledger-20261029/aprl'
    assert sha(binary)=='c48c376929fbdabde4f9cbf08145aa11ea83f86e11839aec24e5f26f3601977b'
    base=yaml.safe_load((ws/'results/energy_ledger_test_record_20261030/tests/01_single_node/on/run.yaml').read_text())
    scripts=root/'scripts';scripts.mkdir()
    for name in ('model.py','prepare.py','run.py','analyze.py','launcher.py','run.sbatch'):shutil.copy2(here/name,scripts/name)
    shutil.copy2(here.parent/'discrete_energy_ledger/report.py',scripts/'report_discrete.py')
    shutil.copy2(here.parent/'discrete_energy_ledger/review_discrete.py',scripts/'review_discrete.py')
    shutil.copy2(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py',scripts/'relative_ledger.py')
    cases=[]
    for kind in a.cases:
        directory=root/kind;directory.mkdir();(directory/'reference_input').mkdir()
        model=Sphere(kind);entry=model.event(0.);count=len(model.weights)
        dtype=[('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')]
        data=np.zeros(count,dtype=dtype);data['plane_position_m']=entry['xy_m'];data['arrival_time_offset_s']=entry['lab_time_s'];data['proper_velocity']=entry['lab_u'];data['source_id']=np.arange(1,count+1);data['macro_weight']=model.weights*count
        inputfile=directory/'reference_input/particles.h5'
        with h5py.File(inputfile,'x') as h:
            g=h.create_group('particles');g.attrs['format_version']=4
            g.attrs['scope']='Theory worldline crossings of fixed lab plane; weighted quiet-start continuum quadrature'
            g.create_dataset('records',data=data)
        # Covers the coherent inward-chirp turning point (~1.71 t0), with
        # a post-turn observation. The first pilot measured the audit cost;
        # 2.2 t0 keeps each physical comparison within its bounded budget.
        target_taus=[.6,1.,1.4,1.75,2.2]
        planes=[{'name':'entry','z':0.}]+[{'name':f'probe-{i+1}','z':float(model.manifest()['beam_gamma_lab']*np.sqrt(1-model.manifest()['beam_gamma_lab']**-2)*c*(t*model.t0-model.entry_center_time))} for i,t in enumerate(target_taus)]
        assert all(p['z']>0 for p in planes[1:])
        expected={}
        for plane in planes:
            event=model.event(plane['z'])
            for key,v in event.items():expected[plane['name']+'__'+key]=v
        expected['weights']=model.weights;expected['shell']=model.shell;expected['directions']=model.directions
        np.savez(directory/'reference_input/expected_events.npz',**expected)
        geo=geometry(data,planes[-1]['z']);cfg=copy.deepcopy(base)
        cfg['output']['manifest']=str(directory/'run-manifest.yaml')
        cfg['mesh']['cells']=[160,160,geo['reference_cells']];cfg['mesh']['duration']=geo['duration_s']
        cfg['mesh']['moving_window']['reference_directory']=str(directory/'reference')
        cfg['mesh']['moving_window']['cells']=geo['active_cells']
        cfg['beam']['reference'].update(input_plane_z=0.,initial_center_z=geo['initial_center_lab_z_m'],activation_plane_z=0.)
        cfg['beam']['input'].update(file=str(inputfile),electrons=model.electrons)
        cfg['sources']={'incident_waves':[],'magnetic_elements':[]}
        cfg['detectors']['directory']=str(directory/'detectors');cfg['detectors']['field_planes']=[]
        cfg['detectors']['particle_planes']=[dict(p,buffer_records=32768,compression=0) for p in planes]
        cfg['trajectory'].update(directory=str(directory/'trajectory'),rhythm=20*DT,buffer_records=8192)
        cfg['energy_ledger'].update(directory=str(directory/'ledger'),enabled=True,discrete_operator_audit=True)
        cfg['stop']={'mode':'reference-center-z','z':geo['stop_reference_lab_z_m']}
        (directory/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        manifest={'model':model.manifest(),'geometry':geo,'planes':planes,'macroparticles':count,
           'weights_sum':float(np.sum(model.weights)),'input_sha256':sha(inputfile),'timeout_s':1200,
           'initial_transverse_centroid_m':np.sum(data['plane_position_m']*model.weights[:,None],axis=0).tolist(),
           'reference_input_note':'Each particle has its own analytic crossing time, position, momentum. Common cold/chirped origin is not reset at injection.',
           'scientific_acceptance':'Exploratory comparison; lifecycle/input checks do not imply agreement with the physical reference.'}
        assert abs(manifest['weights_sum']-1)<1e-12
        save(directory/'case.json',manifest);cases.append({'name':kind,'directory':str(directory),'timeout_s':1200})
    plan={'root':str(root),'workspace':str(ws),'binary':str(binary),'binary_sha256':sha(binary),
       'build':str(binary.parent),'binding_probe':str(ws/'build/cluster-tools/binding_probe'),
       'cases':cases,'ranks':16,'threads':8,'nodes':1,'active_budget_s':3700,
       'scope':'First exploratory spherical space-charge physical comparisons; no core changes, no acceptance inferred from execution completion'}
    save(root/'plan.json',plan)
    paths=[root/'plan.json',binary,Path(plan['binding_probe'])]
    paths += [q for q in scripts.rglob('*') if q.is_file()]
    for item in cases:
        d=Path(item['directory']);paths += [d/'run.yaml',d/'case.json',*list((d/'reference_input').iterdir())]
    save(root/'frozen.json',{str(q):sha(q) for q in paths})
    print(json.dumps({item['name']:json.loads((Path(item['directory'])/'case.json').read_text())['geometry'] for item in cases},indent=2))
if __name__=='__main__':main()
