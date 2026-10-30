#!/usr/bin/env python3
"""Prepare the tested matched-reference / thicker-CPML control, never submit."""
import argparse,json,shutil
from pathlib import Path
import yaml
from prepare import sha,save

def main(parent,root):
    assert json.loads((parent/'preflight/COMPLETION.json').read_text())['pass_all']
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir()
    here=Path(__file__).resolve().parent
    for name in ('run.py','prepare.py','analyze.py','run.sbatch'):
        shutil.copy2(here/name,root/'scripts'/name)
    for name in ('model.py','launcher.py'):shutil.copy2(parent/'scripts'/name,root/'scripts'/name)
    (root/'preflight').symlink_to(parent/'preflight',target_is_directory=True)
    name='matched_domain_thick_cpml';d=root/name;d.mkdir();(d/'momentum').mkdir()
    cfg=yaml.safe_load((parent/'base/run.yaml').read_text())
    assert cfg['mesh']['cells'][:2]==[160,160] and cfg['boundary']['cells'][:2]==[10,10]
    cfg['mesh']['cells'][:2]=[192,192];cfg['boundary']['cells'][:2]=[26,26]
    cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):
        cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    plan=json.loads((parent/'plan.json').read_text());item=dict(plan['cases']['base']);item['directory']=str(d)
    plan.update(root=str(root),cases={name:item},matched_reference_to=str(parent/'base/reference'),required_source=str(parent))
    save(root/'plan.json',plan)
    files=[root/'plan.json',d/'run.yaml',*list((root/'scripts').iterdir()),Path(item['binary']),parent/'preflight/COMPLETION.json',Path(cfg['beam']['input']['file'])]
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})
    save(root/'DESIGN.json',dict(physical_reference_aperture_m=[-.00035,.00035],same_reference_interior_nodes=139,
        base_mesh=160,base_cpml=10,matched_mesh=192,matched_cpml=26,cpml_width_um=[50,130],same_pusher_and_field_solver=True))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();main(a.parent.resolve(),a.root.resolve())
