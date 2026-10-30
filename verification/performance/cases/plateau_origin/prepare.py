#!/usr/bin/env python3
"""Freeze guarded inlet contrasts using the already accepted physics binary."""
import argparse,copy,hashlib,json,shutil
from pathlib import Path
import numpy as np,yaml

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def require(value,message):
    if not value:raise RuntimeError(message)
def main(workspace,root):
    b=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    parent=b/'charge_matched_inlet_20261030';old=json.loads((parent/'plan.json').read_text())
    assert json.loads((parent/'preflight/COMPLETION.json').read_text())['checks_pass']
    assert sha(old['binary'])==old['binary_sha256']
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for n in ('prepare.py','run.py','run.sbatch','review.py','early_force.py'):
        shutil.copy2(here/n,root/'scripts'/n)
    for n in ('launcher.py','analyze_momentum.py','model.py','report_discrete.py','review_discrete.py'):
        shutil.copy2(parent/'scripts'/n,root/'scripts'/n)
    config=yaml.safe_load((parent/'pilot_500/run.yaml').read_text())
    plan=dict(root=str(root),parent=str(parent),binary=old['binary'],binary_sha256=old['binary_sha256'],build=old['build'],
              binding_probe=old['binding_probe'],nominal_dt_s=old['nominal_dt_s'],cases={})
    for guard in (128,256):
        name=f'guard{guard}_500';d=root/name;d.mkdir();(d/'momentum').mkdir()
        card=copy.deepcopy(config);card['mesh']['moving_window']['lower_guard_cells']=guard
        card['mesh']['moving_window']['cells']+=guard-64
        card['output']['manifest']=str(d/'run-manifest.yaml')
        card['mesh']['moving_window']['reference_directory']=str(d/'reference')
        for key in ('detectors','trajectory','energy_ledger'):
            card[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
        (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
        plan['cases'][name]=dict(directory=str(d),guard=guard,steps=500,original=old['cases']['pilot_500']['original'],comparison=str(parent/'pilot_500'))
    save(root/'plan.json',plan)
    inputs=[root/'plan.json',Path(plan['binary']),Path(config['beam']['input']['file']),
            parent/'preflight/COMPLETION.json',*(root/'scripts').iterdir(),*(root/n/'run.yaml' for n in plan['cases'])]
    save(root/'FROZEN.json',{str(p):sha(p) for p in inputs})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace.resolve(),a.root.resolve())
