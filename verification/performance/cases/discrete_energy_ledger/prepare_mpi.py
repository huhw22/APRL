#!/usr/bin/env python3
import argparse,copy,json,shutil
from pathlib import Path
import yaml
from prepare import sha,save
def main():
    p=argparse.ArgumentParser();p.add_argument('--single',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();single=a.single.resolve();root=a.root.resolve()
    old=json.loads((single/'plan.json').read_text());gate=json.loads((single/'COMPLETION.json').read_text())
    assert gate['checks_pass'] and gate['binary_sha256']==old['binary_sha256'] and sha(old['binary'])==old['binary_sha256']
    root.mkdir();plan=copy.deepcopy(old);plan.update(root=str(root),single_node_root=str(single),cases={})
    for name in ('off','on'):
        directory=root/name;directory.mkdir();c=yaml.safe_load((Path(old['cases'][name])/'run.yaml').read_text())
        c['output']['manifest']=str(directory/'run-manifest.yaml');c['mesh']['moving_window']['reference_directory']=str(directory/'reference')
        for k in ('detectors','trajectory','energy_ledger'):c[k]['directory']=str(directory/('ledger' if k=='energy_ledger' else k))
        (directory/'run.yaml').write_text(yaml.safe_dump(c,sort_keys=False));plan['cases'][name]=str(directory)
    scripts=root/'scripts';scripts.mkdir();here=Path(__file__).resolve().parent
    for name in ('mpi.py','mpi.sbatch','review_discrete.py','prepare_mpi.py'):shutil.copy2(here/name,scripts/name)
    shutil.copytree(single/'scripts/transport',scripts/'transport')
    save(root/'plan.json',plan)
    paths=[root/'plan.json',Path(old['binary']),Path(old['operator_checks']),Path(old['binding_probe']),single/'COMPLETION.json']
    paths+=[Path(d)/'run.yaml' for d in plan['cases'].values()];paths+=[p for p in scripts.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    save(root/'frozen.json',{str(p):sha(p) for p in paths});print(root)
if __name__=='__main__':main()
