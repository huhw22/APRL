#!/usr/bin/env python3
import argparse,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np,yaml
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'space_charge_physics'))
from model import DT,G_BOOST,c

def main(root):
    here=Path(__file__).resolve().parent;d=root/'two_nodes_200';d.mkdir(exist_ok=False)
    plan=json.loads((root/'plan.json').read_text());old=json.loads((Path(plan['cases']['surface_500']['original'])/'case.json').read_text())
    cfg=yaml.safe_load((root/'surface_500/run.yaml').read_text());cfg['mesh']['duration']=float(202*DT)
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    cfg['stop']['z']=float(old['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*199.75*DT)
    cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for k in ('trajectory','detectors','energy_ledger'):cfg[k]['directory']=str(d/('ledger' if k=='energy_ledger' else k))
    (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    for n in ('run_two_nodes.py','two_nodes.sbatch','prepare_two_nodes.py'):shutil.copy2(here/n,root/'scripts'/n)
    files=[d/'run.yaml',*(root/'scripts'/n for n in ('run_two_nodes.py','two_nodes.sbatch','prepare_two_nodes.py'))]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (d/'FROZEN.json').write_text(json.dumps({str(p):sha(p) for p in files},indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root.resolve())
