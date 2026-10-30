#!/usr/bin/env python3
"""Freeze charge-scaling cases at fixed physical planes, mesh and cold radius."""
import argparse,copy,hashlib,json,math,shutil,sys
from pathlib import Path
import h5py,numpy as np,yaml

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
def save(p,value):Path(p).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def main(workspace,root):
    here=Path(__file__).resolve().parent;common=here.parent/'space_charge_physics'
    sys.path.insert(0,str(common))
    from model import Sphere,G_BOOST,DZ,DT,c,G_BEAM,B_BEAM,m_e,e
    from prepare import geometry
    workspace=workspace.resolve();root=root.resolve()
    bench=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    baseline=bench/'charge_matched_inlet_20261030/uniform_full'
    original=bench/'space_charge_pic_20261030_r2/uniform'
    base_case=json.loads((original/'case.json').read_text())
    base=yaml.safe_load((baseline/'run.yaml').read_text())
    binary=workspace/'build/charge-matched-inlet-20261030/aprl'
    assert sha(binary)=='79d643bc8ade4db4d9e4488a367cb1bf5ef74e7c602dc3a2094462709aed9435'
    root.mkdir(parents=True,exist_ok=False);scripts=root/'scripts';scripts.mkdir()
    for p in here.iterdir():
        if p.is_file():shutil.copy2(p,scripts/p.name)
    for name in ['model.py','prepare.py','analyze.py','launcher.py']:shutil.copy2(common/name,scripts/name)
    for source,target in [('report.py','report_discrete.py'),('review_discrete.py','review_discrete.py')]:
        shutil.copy2(here.parent/'discrete_energy_ledger'/source,scripts/target)
    shutil.copy2(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py',scripts/'relative_ledger.py')
    plan=dict(root=str(root),workspace=str(workspace),build=str(binary.parent),binary=str(binary),binary_sha256=sha(binary),
              binding_probe=str(workspace/'build/cluster-tools/binding_probe'),baseline=str(baseline),
              original_baseline=str(original),cases={},nominal_dt_s=DT,
              scope='Fixed initial rest radius and cold velocities, lab gamma, physical planes, grid and CPML; vary represented charge only. Plane events follow each charge-specific sphere worldline, so input plane momenta are not forced identical.',
              hypothesis='Measure absolute error and error/theory separately. Do not assume constant percentage, constant eV, or a power law.',
              caveat='Sphere-family scaling is not an error bound for a long nonspherical production bunch.')
    dtype=[('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')]
    theory=[];mc2=np.longdouble(m_e)*np.longdouble(c)**2/np.longdouble(e)
    for factor in [1.,.1,.01,.001]:
        model=Sphere('uniform',electrons=1e8*factor);entry=model.event(0.);count=len(model.weights)
        data=np.zeros(count,dtype=dtype);data['plane_position_m']=entry['xy_m'];data['arrival_time_offset_s']=entry['lab_time_s'];data['proper_velocity']=entry['lab_u'];data['source_id']=np.arange(1,count+1);data['macro_weight']=model.weights*count
        end=model.event(base_case['planes'][-1]['z'])
        u=entry['rest_u'].astype(np.longdouble);v=end['rest_u'].astype(np.longdouble)
        aa=np.sum(u*u,axis=1);bb=np.sum(v*v,axis=1)
        dk=mc2*(bb/(np.sqrt(1+bb)+1)-aa/(np.sqrt(1+aa)+1))
        gain=np.longdouble(G_BEAM)*(dk+np.longdouble(B_BEAM)*mc2*(v[:,2]-u[:,2]))
        theory.append(dict(charge_factor=factor,electrons=model.electrons,charge_pC=model.electrons*e*1e12,
                           fixed_plane_lab_mean_gain_eV=float(np.sum(model.weights*gain)),
                           fixed_plane_rest_mean_gain_eV=float(np.sum(model.weights*dk)),
                           last_plane_z_m=base_case['planes'][-1]['z']))
        if factor==1.:
            with h5py.File(original/'reference_input/particles.h5') as h:
                old=h['particles/records'][:]
            assert data.dtype==old.dtype and data.tobytes()==old.tobytes(),'Baseline sphere inputs differ'
            continue
        name=f'q{factor:g}'.replace('.','p');d=root/name;d.mkdir();(d/'reference_input').mkdir()
        with h5py.File(d/'reference_input/particles.h5','x') as h:
            g=h.create_group('particles');g.attrs['format_version']=4
            g.attrs['scope']='Charge-specific worldline crossings; common cold rest initial distribution'
            g.create_dataset('records',data=data)
        expected={}
        for plane in base_case['planes']:
            for key,value in model.event(plane['z']).items():expected[plane['name']+'__'+key]=value
        expected.update(weights=model.weights,shell=model.shell,directions=model.directions)
        np.savez(d/'reference_input/expected_events.npz',**expected)
        geo=geometry(data,base_case['planes'][-1]['z'])
        # Fix the actual reference-grid and PIC sizes across all charge cases.
        full=base['mesh']['cells'][2];active=base['mesh']['moving_window']['cells']
        origin=-.5*full*DZ;offset=math.floor((geo['box_translation_m']-origin)/DZ)-64
        step=math.ceil((base['stop']['z']-geo['reference_initial_lab_z_m'])/(math.sqrt(G_BOOST*G_BOOST-1)*DZ))
        uz=data['proper_velocity'][:,2].astype(np.longdouble);boost=np.longdouble(G_BOOST)
        plus=(np.sqrt(1+uz*uz)+uz)/(boost+np.sqrt(boost*boost-1));held=(plus-1/plus)/2;held/=np.sqrt(1+held*held)
        mean_held=float(np.mean(held))
        lower=offset-(1+mean_held)*(step+3)-1;upper=offset+active
        assert lower>=0 and upper<=full-1,(name,lower,upper,full)
        geo.update(reference_cells=full,active_cells=active,initial_offset=offset,origin_reference_m=origin,
                   expected_stop_step=step,duration_s=base['mesh']['duration'],
                   minimum_reference_query_cells=lower,maximum_reference_query_cells=upper,
                   stop_reference_lab_z_m=base['stop']['z'])
        cfg=copy.deepcopy(base);cfg['output']['manifest']=str(d/'run-manifest.yaml')
        cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        cfg['beam']['reference']['initial_center_z']=geo['initial_center_lab_z_m']
        cfg['beam']['input'].update(file=str(d/'reference_input/particles.h5'),electrons=model.electrons)
        for key in ['detectors','trajectory','energy_ledger']:cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        case=dict(model=model.manifest(),geometry=geo,planes=base_case['planes'],macroparticles=count,
                  charge_factor=factor,weights_sum=float(np.sum(model.weights)),input_sha256=sha(d/'reference_input/particles.h5'),
                  timeout_s=1400,reference_input_note=plan['scope'],scientific_acceptance='Exploratory scaling; no production acceptance')
        save(d/'case.json',case);plan['cases'][name]=dict(directory=str(d),factor=factor,timeout_s=1400,steps=step)
    save(root/'THEORY.json',theory);save(root/'plan.json',plan)
    paths=[root/'THEORY.json',root/'plan.json',binary,Path(plan['binding_probe']),original/'reference_input/particles.h5',baseline/'analysis/RESULT.json',*scripts.iterdir()]
    for item in plan['cases'].values():
        d=Path(item['directory']);paths.extend([d/'run.yaml',d/'case.json',*list((d/'reference_input').iterdir())])
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths if p.is_file()})
    print(json.dumps(theory,indent=2));print('Prepared',root)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace,a.root)
