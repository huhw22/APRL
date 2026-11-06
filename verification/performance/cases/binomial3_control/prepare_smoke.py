"""Freeze old-mode invariance and matched-three-point full-flow controls."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import sys

import yaml


def main(workspace, output, binary, build_record):
    workspace, output, binary, build_record = map(lambda p: Path(p).resolve(),
        (workspace, output, binary, build_record))
    output.mkdir(parents=True, exist_ok=False)
    scripts = output/'scripts'; scripts.mkdir()
    cases_root = Path(__file__).resolve().parent.parent
    for folder, name, target in [
        ('binomial3_control','run_smoke.py','run_smoke.py'),
        ('energy_critical_control','prepare.py','prepare.py'),
        ('space_charge_physics','launcher.py','launcher.py'),
        ('source_coupling_integration','binding.py','binding.py'),
        ('source_coupling_integration','analyze_closed.py','analyze_closed.py'),
        ('discrete_energy_ledger','report.py','report_discrete.py'),
        ('discrete_energy_ledger','review_discrete.py','review_discrete.py')]:
        shutil.copy2(cases_root/folder/name,scripts/target)
    sys.path.insert(0,str(scripts))
    from prepare import save, sha
    original = workspace/'results/longitudinal_main_flow_20261103/flow-r2/drift-linear'
    old = workspace/'build/bf-force-probe-20261105/aprl'
    cases=[]
    for name, program, shape, nodes in [
        ('old-linear',old,'linear',1),('new-linear',binary,'linear',1),
        ('old-notch',old,'linear-notch',1),('new-notch',binary,'linear-notch',1),
        ('binomial3-single',binary,'linear-binomial3',1),
        ('binomial3-dual',binary,'linear-binomial3',2)]:
        directory=output/name;directory.mkdir()
        def rewrite(value):
            if isinstance(value,dict): return {k:rewrite(v) for k,v in value.items()}
            if isinstance(value,list): return [rewrite(v) for v in value]
            return value.replace(str(original),str(directory)) if isinstance(value,str) else value
        base=yaml.safe_load((original/'run.yaml').read_text())
        card=rewrite(base);card['beam']=copy.deepcopy(base['beam'])
        card['mesh']['longitudinal_shape']=shape
        card['detectors']['substep_particle_planes']=True
        for folder in ['detectors','trajectory','ledger']:(directory/folder).mkdir()
        # The reference generator requires a nonexistent output directory.
        assert not (directory/'reference').exists()
        (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
        cases.append(dict(name=name,directory=str(directory),binary=str(program),
            shape=shape,nodes=nodes,ranks=16*nodes,threads=8,limit_seconds=300))
    plan=dict(workspace=str(workspace),directory=str(output),binary=str(binary),
        build_record=str(build_record),expected_particles=31675,expected_steps=257,
        cases=cases,physical_acceptance=False,
        scope='Unchanged drift input and physics; full-flow charge/work checks, old paths bitwise, dual-node consistency.')
    save(output/'plan.json',plan)
    paths=[old,output/'plan.json',*scripts.iterdir()]
    for case in cases:
        path=Path(case['directory'])/'run.yaml'
        paths.extend([path,Path(yaml.safe_load(path.read_text())['beam']['input']['file'])])
    save(output/'FROZEN_INPUTS.json',{str(p):sha(p) for p in paths if p.is_file()})
    print(output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('workspace');p.add_argument('output')
    p.add_argument('--binary',required=True);p.add_argument('--build-record',required=True)
    a=p.parse_args();main(a.workspace,a.output,a.binary,a.build_record)
