#!/usr/bin/env python3
"""Controlled transverse-resolution check with the same physical aperture."""
import argparse,copy,json,shutil,sys
from pathlib import Path
import numpy as np,yaml
from prepare import sha,save,require

def main(workspace,root):
    b=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks';prior=b/'charge_inlet_plateau_origin_20261030_r3';parent=b/'charge_matched_inlet_20261030'
    require(json.loads((prior/'guard128_500/COMPLETION.json').read_text())['checks_pass'],'Base guard incomplete')
    old=json.loads((parent/'plan.json').read_text());require(sha(old['binary'])==old['binary_sha256'],'Binary changed')
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for n in ('prepare.py','prepare_refine.py','run.py','run.sbatch'):
        shutil.copy2(here/n,root/'scripts'/n)
    for n in ('launcher.py','analyze_momentum.py','model.py','report_discrete.py','review_discrete.py'):
        shutil.copy2(parent/'scripts'/n,root/'scripts'/n)
    sys.path.insert(0,str(root/'scripts'));from model import DT,G_BOOST,c
    d=root/'xy240_300';d.mkdir();(d/'momentum').mkdir()
    card=copy.deepcopy(yaml.safe_load((prior/'guard128_500/run.yaml').read_text()))
    card['mesh']['cells'][:2]=[240,240]
    card['mesh']['cell_size'][:2]=[float(160*5e-6/240)]*2
    card['boundary']['cells'][:2]=[15,15]
    card['mesh']['duration']=float(302*DT)
    meta=json.loads((Path(old['cases']['pilot_500']['original'])/'case.json').read_text())
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    card['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*299.75*DT)
    card['output']['manifest']=str(d/'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):
        card[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    plan=dict(root=str(root),parent=str(parent),build=old['build'],binary=old['binary'],binary_sha256=old['binary_sha256'],binding_probe=old['binding_probe'],nominal_dt_s=old['nominal_dt_s'],
        cases={'xy240_300':dict(directory=str(d),guard=128,steps=300,timeout_s=1100,original=old['cases']['pilot_500']['original'],comparison=str(prior/'guard128_500'))})
    save(root/'plan.json',plan)
    paths=[root/'plan.json',Path(plan['binary']),Path(card['beam']['input']['file']),parent/'preflight/COMPLETION.json',
           *(p for p in (root/'scripts').iterdir() if p.is_file()),d/'run.yaml']
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace.resolve(),a.root.resolve())
