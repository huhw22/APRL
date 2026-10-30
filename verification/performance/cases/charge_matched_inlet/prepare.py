#!/usr/bin/env python3
import argparse,copy,hashlib,json,shutil
from pathlib import Path
import yaml
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def require(v,msg):
    if not v:raise RuntimeError(msg)
def main(w,root,build):
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir();here=Path(__file__).resolve().parent;p=here.parent
    for n in ('prepare.py','run.py','run.sbatch','DESIGN.md'):shutil.copy2(here/n,root/'scripts'/n)
    for n in ('model.py','launcher.py'):shutil.copy2(p/'space_charge_physics'/n,root/'scripts'/n)
    shutil.copy2(p/'space_charge_physics/analyze.py',root/'scripts/analyze_particles.py')
    shutil.copy2(p/'longitudinal_momentum_origin/analyze.py',root/'scripts/analyze_momentum.py')
    for src,dst in [('report.py','report_discrete.py'),('review_discrete.py','review_discrete.py')]:shutil.copy2(p/'discrete_energy_ledger'/src,root/'scripts'/dst)
    shutil.copy2(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py',root/'scripts/relative_ledger.py')
    b=w/'results/energy_ledger_test_record_20261030/physical_benchmarks';prior=b/'reference_surface_current_20261030';origin=b/'surface_residual_origin_20261030_r2'
    original={k:b/'space_charge_pic_20261030_r2'/k for k in ('uniform','converging')};original['gaussian']=b/'space_charge_pic_gaussian_20261030_r4/gaussian'
    cases={}
    for name,kind,mode,source in [('default_100','uniform','default',prior/'default_100'),('pilot_500','converging','pilot',origin/'base'),*[(k+'_full',k,'full',prior/(k+'_full')) for k in original]]:
        d=root/name;d.mkdir();(d/'momentum').mkdir();(d/'diagnostic').mkdir();cfg=copy.deepcopy(yaml.safe_load((source/'run.yaml').read_text()))
        if mode!='default':cfg['mesh']['moving_window']['inlet_charge_match']=True
        cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        for key in ('detectors','trajectory','energy_ledger'):cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        if mode=='full':shutil.copy2(original[kind]/'case.json',d/'case.json');(d/'reference_input').symlink_to(original[kind]/'reference_input',target_is_directory=True)
        steps=100 if mode=='default' else 500 if mode=='pilot' else json.loads((original[kind]/'case.json').read_text())['geometry']['expected_stop_step']
        cases[name]=dict(kind=kind,mode=mode,directory=str(d),steps=steps,original=str(original[kind]),comparison=str(source),timeout_s=1650 if mode=='full' else 900)
    plan=dict(root=str(root),workspace=str(w),build=str(build),binary=str(build/'aprl'),binary_sha256=sha(build/'aprl'),cases=cases,
        binding_probe=str(w/'build/cluster-tools/binding_probe'),nominal_dt_s=1.2543262415359034e-6/299792458.)
    save(root/'plan.json',plan);files=[root/'plan.json',build/'aprl',*list((root/'scripts').iterdir())]
    for case in cases.values():
        d=Path(case['directory']);files.append(d/'run.yaml');files.append(Path(yaml.safe_load((d/'run.yaml').read_text())['beam']['input']['file']))
        if case['mode']=='full':files.extend([d/'case.json',d/'reference_input/expected_events.npz'])
    save(root/'FROZEN.json',{str(p):sha(p) for p in files});print(root)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve(),a.build.resolve())
