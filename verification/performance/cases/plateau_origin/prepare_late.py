#!/usr/bin/env python3
"""One identical 64-cell run through step 1000 with sparse pusher-force trace."""
import argparse,copy,json,shutil,sys
from pathlib import Path
import numpy as np,yaml
from prepare import sha,save,require

def main(workspace,root):
    b=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    parent=b/'charge_matched_inlet_20261030';old=json.loads((parent/'plan.json').read_text());require(sha(old['binary'])==old['binary_sha256'],'Binary changed')
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for n in ('prepare.py','prepare_late.py','run.py','run.sbatch','early_force.py'):
        shutil.copy2(here/n,root/'scripts'/n)
    for n in ('launcher.py','analyze_momentum.py','model.py','report_discrete.py','review_discrete.py'):
        shutil.copy2(parent/'scripts'/n,root/'scripts'/n)
    d=root/'guard64_1000';d.mkdir();(d/'momentum').mkdir()
    card=copy.deepcopy(yaml.safe_load((parent/'pilot_500/run.yaml').read_text()))
    sys.path.insert(0,str(root/'scripts'))
    from model import DT,G_BOOST,c
    meta=json.loads((Path(old['cases']['pilot_500']['original'])/'case.json').read_text())
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    card['mesh']['duration']=float(1002*DT)
    card['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*999.75*DT)
    card['output']['manifest']=str(d/'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):
        card[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    plan=dict(root=str(root),parent=str(parent),build=old['build'],binary=old['binary'],binary_sha256=old['binary_sha256'],
        binding_probe=old['binding_probe'],nominal_dt_s=old['nominal_dt_s'],cases={'guard64_1000':dict(directory=str(d),guard=64,steps=1000,timeout_s=1100,original=old['cases']['pilot_500']['original'],comparison=str(parent/'pilot_500'))})
    save(root/'plan.json',plan)
    paths=[root/'plan.json',Path(plan['binary']),Path(card['beam']['input']['file']),parent/'preflight/COMPLETION.json',
           *(p for p in (root/'scripts').iterdir() if p.is_file()),d/'run.yaml']
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace.resolve(),a.root.resolve())
