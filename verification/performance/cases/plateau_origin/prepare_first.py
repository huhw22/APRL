#!/usr/bin/env python3
"""Freeze one first-push trace with distant inlet and unchanged source state."""
import argparse,copy,json,shutil,sys
from pathlib import Path
import numpy as np,yaml
from prepare import sha,save,require

def main(workspace,root):
    bench=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    prior=bench/'charge_inlet_plateau_origin_20261030_r3';parent=bench/'charge_matched_inlet_20261030'
    old=json.loads((parent/'plan.json').read_text())
    audit=json.loads((bench/'charge_matched_source_audit_20261030_r2/plan.json').read_text())
    require(sha(old['binary'])==old['binary_sha256'] and sha(audit['binary'])==audit['binary_sha256'],'Binary changed')
    require(json.loads((prior/'guard256_500/COMPLETION.json').read_text())['checks_pass'],'Prior run incomplete')
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for n in ('prepare.py','prepare_first.py','run_first.py','run_first.sbatch'):
        shutil.copy2(here/n,root/'scripts'/n)
    shutil.copy2(parent/'scripts/launcher.py',root/'scripts/launcher.py')
    shutil.copy2(parent/'scripts/model.py',root/'scripts/model.py')
    d=root/'guard256_first_140';d.mkdir();(d/'momentum').mkdir()
    card=copy.deepcopy(yaml.safe_load((prior/'guard256_500/run.yaml').read_text()))
    sys.path.insert(0,str(root/'scripts'));from model import DT,G_BOOST,c
    meta=json.loads((Path(old['cases']['pilot_500']['original'])/'case.json').read_text())
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    card['mesh']['duration']=float(142*DT)
    card['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*139.75*DT)
    card['output']['manifest']=str(d/'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):
        card[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    plan=dict(root=str(root),binary=audit['binary'],binary_sha256=audit['binary_sha256'],build=audit['build'],binding_probe=old['binding_probe'],
              nominal_dt_s=old['nominal_dt_s'],steps=140,case=str(d),comparison=str(prior/'guard256_500'),timeout_s=300)
    save(root/'plan.json',plan)
    paths=[root/'plan.json',Path(plan['binary']),Path(card['beam']['input']['file']),
           *(p for p in (root/'scripts').iterdir() if p.is_file()),d/'run.yaml']
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path)
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve())
