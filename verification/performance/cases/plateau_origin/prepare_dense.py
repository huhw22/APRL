#!/usr/bin/env python3
"""Freeze a same-card 1000-step run with a full-step read-only force trace."""
import argparse,copy,json,shutil,sys
from pathlib import Path
import yaml
from prepare import sha,save,require

def main(workspace,root):
    b=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    prior=b/'charge_inlet_plateau_late_force_20261030'
    old=json.loads((prior/'plan.json').read_text());require(sha(old['binary'])==old['binary_sha256'],'Binary changed')
    build=workspace/'build/plateau-force-stride-r2-20261030';binary=build/'aprl'
    require(binary.exists(),'Diagnostic stride binary missing')
    require(json.loads((prior/'guard64_1000/COMPLETION.json').read_text())['checks_pass'],'Prior run incomplete')
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for n in ('prepare.py','prepare_dense.py','run_dense.py','run_dense.sbatch'):
        shutil.copy2(here/n,root/'scripts'/n)
    shutil.copy2(prior/'scripts/launcher.py',root/'scripts/launcher.py')
    d=root/'guard64_dense_1000';d.mkdir();(d/'momentum').mkdir()
    card=copy.deepcopy(yaml.safe_load((prior/'guard64_1000/run.yaml').read_text()))
    card['output']['manifest']=str(d/'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):
        card[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    plan=dict(root=str(root),binary=str(binary),binary_sha256=sha(binary),build=str(build),binding_probe=old['binding_probe'],
              nominal_dt_s=old['nominal_dt_s'],probe_time_divisor=1,steps=1000,
              case=str(d),comparison=str(prior/'guard64_1000'),timeout_s=1620)
    save(root/'plan.json',plan)
    paths=[root/'plan.json',Path(plan['binary']),Path(card['beam']['input']['file']),
           *(p for p in (root/'scripts').iterdir() if p.is_file()),d/'run.yaml']
    save(root/'FROZEN.json',{str(p):sha(p) for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path)
    a=p.parse_args();main(a.workspace.resolve(),a.root.resolve())
