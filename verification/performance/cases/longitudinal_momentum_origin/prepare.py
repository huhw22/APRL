#!/usr/bin/env python3
import argparse,copy,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np,yaml
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'space_charge_physics'))
from model import DT,G_BOOST,c

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def main(w,root,build,required_source=None):
    original=w/'results/energy_ledger_test_record_20261030/physical_benchmarks/space_charge_pic_20261030_r2/uniform'
    base=yaml.safe_load((original/'run.yaml').read_text());old=json.loads((original/'case.json').read_text())
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for name in ('prepare.py','analyze.py','run.py','run.sbatch','render.py','ideal_pusher.py','first_step_boundary_probe.cpp','DESIGN.md'):
        shutil.copy2(here/name,root/'scripts'/name)
    for name in ('model.py','launcher.py'):shutil.copy2(here.parent/'space_charge_physics'/name,root/'scripts'/name)
    # launcher imports prepare.require/save; these small helpers are local below.
    cases={}
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    for name in ('base','inlet_guard_256','transverse_240'):
        d=root/name;d.mkdir();(d/'momentum').mkdir();cfg=copy.deepcopy(base)
        cfg['output']['manifest']=str(d/'run-manifest.yaml')
        cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        cfg['mesh']['duration']=float(502*DT)
        cfg['stop']['z']=float(old['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*499.75*DT)
        cfg['trajectory']['directory']=str(d/'trajectory')
        cfg['detectors']['directory']=str(d/'detectors')
        cfg['detectors']['particle_planes']=[p for p in cfg['detectors']['particle_planes'] if p['name']=='entry']
        cfg['energy_ledger']['directory']=str(d/'ledger')
        assert cfg['stop']['z']>max(cfg['beam']['reference']['initial_center_z'],max(p['z'] for p in cfg['detectors']['particle_planes']))
        if name=='inlet_guard_256':
            cfg['mesh']['moving_window']['lower_guard_cells']=256
            cfg['mesh']['moving_window']['cells']+=192
        if name=='transverse_240':cfg['mesh']['cells'][:2]=[240,240]
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        cases[name]=dict(directory=str(d),original=str(original),timeout_s=1200,
          binary=str(build/'aprl'),build=str(build),binary_sha256=sha(build/'aprl'),expected_steps=500)
    # A valid downstream stop requires >64 steps for this frame origin.
    for name,enabled in [('probe_off',False),('probe_on',True)]:
        d=root/name;d.mkdir();(d/'momentum').mkdir();cfg=copy.deepcopy(base)
        cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['duration']=float(102*DT)
        cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        cfg['stop']['z']=float(old['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*99.75*DT)
        cfg['trajectory']['directory']=str(d/'trajectory');cfg['detectors']['directory']=str(d/'detectors')
        cfg['detectors']['particle_planes']=[p for p in cfg['detectors']['particle_planes'] if p['name']=='entry']
        cfg['energy_ledger']['directory']=str(d/'ledger')
        assert cfg['stop']['z']>max(cfg['beam']['reference']['initial_center_z'],max(p['z'] for p in cfg['detectors']['particle_planes']))
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        cases[name]=dict(directory=str(d),original=str(original),timeout_s=120,binary=str(build/'aprl'),
          build=str(build),binary_sha256=sha(build/'aprl'),expected_steps=100,probe=enabled)
    save(root/'plan.json',dict(workspace=str(w),root=str(root),build=str(build),nominal_dt_s=DT,
      original=str(original),cases=cases,binding_probe=str(w/'build/cluster-tools/binding_probe'),required_source=str(required_source) if required_source else None))
    files=[root/'plan.json',*list((root/'scripts').iterdir()),*(Path(v['directory'])/'run.yaml' for v in cases.values()),build/'aprl',Path(base['beam']['input']['file'])]
    if required_source:
        prior=json.loads((required_source/'plan.json').read_text())
        assert prior['cases']['probe_off']['binary_sha256']==sha(build/'aprl')
        assert '100% tests passed, 0 tests failed out of 77' in (required_source/'preflight/required.log').read_text()
        files.extend([required_source/'preflight/required.log',required_source/'plan.json'])
    save(root/'FROZEN.json',{str(f):sha(f) for f in files})
    print(root)

def require(v,message):
    if not v:raise RuntimeError(message)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--required-source',type=Path)
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve(),a.build.resolve(),a.required_source.resolve() if a.required_source else None)
