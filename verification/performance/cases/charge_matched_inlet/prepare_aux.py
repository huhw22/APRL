#!/usr/bin/env python3
import argparse,copy,json,shutil
from pathlib import Path
import numpy as np,yaml
from prepare import sha,save
def main(parent,root):
    root.mkdir(parents=True,exist_ok=False);shutil.copytree(parent/'scripts',root/'scripts',ignore=shutil.ignore_patterns('__pycache__'));here=Path(__file__).resolve().parent
    for n in ('prepare_aux.py','run_two_nodes.py','two_nodes.sbatch'):shutil.copy2(here/n,root/'scripts'/n)
    (root/'preflight').symlink_to(parent/'preflight',target_is_directory=True)
    plan=copy.deepcopy(json.loads((parent/'plan.json').read_text()));plan['root']=str(root);plan['comparison']=str(parent/'pilot_500');pilot=plan['cases']['pilot_500'];plan['cases']={}
    from model import DT,G_BOOST,c
    meta=json.loads((Path(pilot['original'])/'case.json').read_text());beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    for name in ('guard32_500','two_nodes_200'):
        d=root/name;d.mkdir();(d/'momentum').mkdir();(d/'diagnostic').mkdir();cfg=copy.deepcopy(yaml.safe_load((parent/'pilot_500/run.yaml').read_text()))
        cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
        for key in ('detectors','trajectory','energy_ledger'):cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
        if name=='guard32_500':cfg['mesh']['moving_window']['lower_guard_cells']=32;cfg['mesh']['moving_window']['cells']-=32
        else:
            cfg['mesh']['duration']=float(202*DT);cfg['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*199.75*DT)
        (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
        plan['cases'][name]=dict(pilot,directory=str(d),comparison=str(Path(pilot['original'])),steps=500 if name=='guard32_500' else 200)
    save(root/'plan.json',plan);files=[root/'plan.json',Path(plan['binary']),*[p for p in (root/'scripts').iterdir() if p.is_file()],*(root/n/'run.yaml' for n in plan['cases'])]
    frozen={str(p):sha(p) for p in files};save(root/'FROZEN.json',frozen);save(root/'two_nodes_200/FROZEN.json',{})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.parent.resolve(),a.root.resolve())
