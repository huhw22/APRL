#!/usr/bin/env python3
"""Freeze short read-only first-push audits; never rewrite the original runs."""
import argparse,copy,json,shutil
from pathlib import Path
import numpy as np,yaml
from prepare import sha,save

def main(parent,root,build):
    root.mkdir(parents=True,exist_ok=False)
    shutil.copytree(parent/'scripts',root/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
    here=Path(__file__).resolve().parent
    for n in ('prepare_source_audit.py','run_source_audit.py','source_audit.sbatch'):
        shutil.copy2(here/n,root/'scripts'/n)
    plan=copy.deepcopy(json.loads((parent/'plan.json').read_text()))
    plan.update(root=str(root),parent=str(parent),build=str(build),binary=str(build/'aprl'),cases={})
    from model import DT,G_BOOST,c
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2))
    for kind,steps in [('uniform',150),('gaussian',270)]:
        original=parent/(kind+'_full');meta=json.loads((original/'case.json').read_text())
        for mode in ('off','on'):
            name=kind+'_'+mode;d=root/name;d.mkdir();(d/'momentum').mkdir()
            cfg=copy.deepcopy(yaml.safe_load((original/'run.yaml').read_text()))
            cfg['output']['manifest']=str(d/'run-manifest.yaml')
            cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
            cfg['mesh']['moving_window']['inlet_charge_match']=mode=='on'
            cfg['mesh']['duration']=float((steps+2)*DT)
            cfg['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*(steps-.25)*DT)
            # Short source audit needs only the entry plane, as in the pilot.
            cfg['detectors']['particle_planes']=cfg['detectors']['particle_planes'][:1]
            for key in ('detectors','trajectory','energy_ledger'):
                cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
            (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
            comparison=original if mode=='on' else parent.parent/'reference_surface_current_20261030'/(kind+'_full')
            plan['cases'][name]=dict(directory=str(d),comparison=str(comparison),steps=steps)
    plan['binary_sha256']=sha(build/'aprl');save(root/'plan.json',plan)
    files=[root/'plan.json',build/'aprl',*[p for p in (root/'scripts').iterdir() if p.is_file()],*(root/n/'run.yaml' for n in plan['cases'])]
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('root',type=Path);p.add_argument('build',type=Path);a=p.parse_args();main(a.parent.resolve(),a.root.resolve(),a.build.resolve())
