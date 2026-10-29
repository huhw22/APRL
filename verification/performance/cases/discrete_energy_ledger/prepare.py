#!/usr/bin/env python3
import argparse,copy,hashlib,json,shutil
from pathlib import Path
import yaml

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def save(path,obj):Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();ws=a.workspace.resolve();build=a.build.resolve();root=a.root.resolve();root.mkdir()
    baseline=ws/'results/energy_transport_audit_20261029/campaign-r3/science/free-drift-r1-new-on'
    oldplan=json.loads((ws/'results/energy_transport_audit_20261029/campaign-r3/plan.json').read_text())
    card=yaml.safe_load((baseline/'run.yaml').read_text());assert not card['sources']['magnetic_elements'] and not card['sources']['incident_waves']
    cases={}
    for name in ('off','on'):
        directory=root/name;directory.mkdir();c=copy.deepcopy(card)
        c['energy_ledger']['discrete_operator_audit']=name=='on'
        c['output']['manifest']=str(directory/'run-manifest.yaml');c['mesh']['moving_window']['reference_directory']=str(directory/'reference')
        for key in ('detectors','trajectory','energy_ledger'):c[key]['directory']=str(directory/('ledger' if key=='energy_ledger' else key))
        (directory/'run.yaml').write_text(yaml.safe_dump(c,sort_keys=False));cases[name]=str(directory)
    scripts=root/'scripts';scripts.mkdir();here=Path(__file__).resolve().parent
    for name in ('prepare.py','run.py','run.sbatch','review_discrete.py','operator_checks.cpp'):shutil.copy2(here/name,scripts/name)
    transport=scripts/'transport';transport.mkdir()
    for name in ('prepare.py','run.py','review.py'):shutil.copy2(here.parent/'energy_transport_audit'/name,transport/name)
    plan={'workspace':str(ws),'root':str(root),'build':str(build),'binary':str(build/'aprl'),
          'binary_sha256':sha(build/'aprl'),'operator_checks':str(build/'operator_checks'),
          'binding_probe':oldplan['binding_probe'],'baseline':str(baseline),'particles':31675,'cases':cases}
    save(root/'plan.json',plan)
    paths=[root/'plan.json',build/'aprl',build/'operator_checks',Path(plan['binding_probe']),Path(card['beam']['input']['file'])]
    paths+=[Path(v)/'run.yaml' for v in cases.values()];paths+=[p for p in scripts.rglob('*') if p.is_file()]
    save(root/'frozen.json',{str(p):sha(p) for p in paths});print(root)
if __name__=='__main__':main()
