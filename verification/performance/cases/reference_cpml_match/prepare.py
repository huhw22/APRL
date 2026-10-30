#!/usr/bin/env python3
"""Bounded same-input reference-boundary experiments; never submit jobs."""
import argparse,copy,hashlib,json,shutil
from pathlib import Path
import numpy as np,yaml

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def require(v,msg):
    if not v:raise RuntimeError(msg)

def main(w,root,build):
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent;phys=here.parent/'space_charge_physics'
    for n in ('prepare.py','run.py','run.sbatch'):shutil.copy2(here/n,root/'scripts'/n)
    for n in ('model.py','launcher.py'):shutil.copy2(phys/n,root/'scripts'/n)
    shutil.copy2(phys/'analyze.py',root/'scripts/analyze_particles.py')
    shutil.copy2(here.parent/'longitudinal_momentum_origin/analyze.py',root/'scripts/analyze_momentum.py')
    shutil.copy2(here.parent/'discrete_energy_ledger/report.py',root/'scripts/report_discrete.py')
    shutil.copy2(here.parent/'discrete_energy_ledger/review_discrete.py',root/'scripts/review_discrete.py')
    shutil.copy2(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py',root/'scripts/relative_ledger.py')
    b=w/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    original={k:b/'space_charge_pic_20261030_r2'/k for k in ('uniform','converging')}
    original['gaussian']=b/'space_charge_pic_gaussian_20261030_r4/gaussian'
    short=b/'longitudinal_momentum_origin_20261030_r3'
    cases={}
    for name,kind,mode,source in [
        ('default_100','uniform','default',short/'probe_off'),
        ('outer_500','uniform','pilot',short/'base'),
        ('uniform_full','uniform','full',original['uniform']),
        ('converging_full','converging','full',original['converging']),
        ('gaussian_full','gaussian','full',original['gaussian'])]:
        d=root/name;d.mkdir();(d/'momentum').mkdir();cfg=copy.deepcopy(yaml.safe_load((source/'run.yaml').read_text()))
        if mode!='default':cfg['mesh']['moving_window']['reference_through_cpml']=True
        cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        cfg['detectors']['directory']=str(d/'detectors');cfg['trajectory']['directory']=str(d/'trajectory');cfg['energy_ledger']['directory']=str(d/'ledger')
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        if mode=='full':
            shutil.copy2(original[kind]/'case.json',d/'case.json');(d/'reference_input').symlink_to(original[kind]/'reference_input',target_is_directory=True)
        steps=100 if mode=='default' else 500 if mode=='pilot' else json.loads((original[kind]/'case.json').read_text())['geometry']['expected_stop_step']
        cases[name]=dict(directory=str(d),kind=kind,mode=mode,original=str(original[kind]),comparison=str(source),
            steps=steps,timeout_s=1600 if mode=='full' else 900,binary=str(build/'aprl'))
    plan=dict(workspace=str(w),root=str(root),build=str(build),binary=str(build/'aprl'),binary_sha256=sha(build/'aprl'),
        binding_probe=str(w/'build/cluster-tools/binding_probe'),nominal_dt_s=1.2543262415359034e-6/299792458.,cases=cases)
    save(root/'plan.json',plan)
    files=[root/'plan.json',build/'aprl',*list((root/'scripts').iterdir())]
    for case in cases.values():
        d=Path(case['directory']);files.append(d/'run.yaml')
        cfg=yaml.safe_load((d/'run.yaml').read_text());files.append(Path(cfg['beam']['input']['file']))
        if case['mode']=='full':files.extend([d/'case.json',d/'reference_input/expected_events.npz'])
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})
    print(root)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve(),a.build.resolve())
