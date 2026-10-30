#!/usr/bin/env python3
"""Same converging input: separate inlet and transverse-boundary influence."""
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
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def require(v,msg):
    if not v:raise RuntimeError(msg)
def main(w,root,followup=False):
    here=Path(__file__).resolve().parent;base=w/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    prior=base/'reference_surface_current_20261030';original=base/'space_charge_pic_20261030_r2/converging'
    build=w/'build/reference-surface-current-20261030'
    assert json.loads((prior/'preflight/COMPLETION.json').read_text())['checks_pass']
    assert sha(build/'aprl')==json.loads((prior/'plan.json').read_text())['binary_sha256']
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    for n in ('prepare.py','run.py','run.sbatch','DESIGN.md'):shutil.copy2(here/n,root/'scripts'/n)
    for n in ('model.py','launcher.py'):shutil.copy2(here.parent/'space_charge_physics'/n,root/'scripts'/n)
    shutil.copy2(here.parent/'longitudinal_momentum_origin/analyze.py',root/'scripts/analyze_momentum.py')
    shutil.copy2(here.parent/'discrete_energy_ledger/report.py',root/'scripts/report_discrete.py')
    shutil.copy2(here.parent/'discrete_energy_ledger/review_discrete.py',root/'scripts/review_discrete.py')
    old=json.loads((original/'case.json').read_text());basecfg=yaml.safe_load((prior/'converging_full/run.yaml').read_text())
    cases={};beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    names=('inlet_guard_128','inlet_guard_256_long') if followup else ('base','inlet_guard_256','transverse_240')
    for name in names:
        steps=1000 if name.endswith('_long') else 500
        d=root/name;d.mkdir();(d/'momentum').mkdir();cfg=copy.deepcopy(basecfg)
        cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['duration']=float((steps+2)*DT)
        cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        cfg['stop']['z']=float(old['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*(steps-.25)*DT)
        cfg['detectors']['directory']=str(d/'detectors');cfg['detectors']['particle_planes']=[p for p in cfg['detectors']['particle_planes'] if p['name']=='entry']
        cfg['trajectory']['directory']=str(d/'trajectory');cfg['energy_ledger']['directory']=str(d/'ledger')
        if name.startswith('inlet_guard_'):
            guard=int(name.split('_')[2]);cfg['mesh']['moving_window']['lower_guard_cells']=guard;cfg['mesh']['moving_window']['cells']+=guard-64
        if name=='transverse_240':cfg['mesh']['cells'][:2]=[240,240]
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        cases[name]=dict(directory=str(d),steps=steps,timeout_s=1500,kind='converging')
    plan=dict(root=str(root),workspace=str(w),build=str(build),binary=str(build/'aprl'),binary_sha256=sha(build/'aprl'),
        nominal_dt_s=DT,original=str(original),prior=str(prior),cases=cases,binding_probe=str(w/'build/cluster-tools/binding_probe'))
    save(root/'plan.json',plan)
    files=[root/'plan.json',*list((root/'scripts').iterdir()),build/'aprl',Path(basecfg['beam']['input']['file']),original/'case.json',prior/'preflight/required.log',prior/'preflight/DEFAULT_BITWISE.json']
    files.extend(Path(v['directory'])/'run.yaml' for v in cases.values())
    save(root/'FROZEN.json',{str(p):sha(p) for p in files});print(root)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--followup',action='store_true')
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve(),a.followup)
