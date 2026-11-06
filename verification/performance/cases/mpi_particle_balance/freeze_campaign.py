"""Install and freeze reviewed executable cards without launching a job."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil

import yaml


def sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):value.update(block)
    return value.hexdigest()


def main(root,binary,record,profile,bounded_prefix=False):
    root=Path(root).resolve();binary=Path(binary).resolve();record=Path(record).resolve()
    assert not (root/'FROZEN_INPUTS.json').exists(),'Refusing to refreeze a campaign'
    plan=json.loads((root/'plan.json').read_text());workspace=Path(plan['workspace'])
    assert bool(plan.get('bounded_duration_prefix',False))==bounded_prefix,'Select the explicit bounded-prefix freezer for a duration prefix'
    assert all(bool(case.get('bounded_duration_prefix',False))==bounded_prefix for case in plan['cases'])
    cases_root=Path(__file__).resolve().parent.parent
    scripts=root/'scripts';scripts.mkdir(exist_ok=False)
    install=[
        ('mpi_particle_balance','prefix_run.py' if bounded_prefix else 'run.py','run.py'),
        ('mpi_particle_balance','prefix_run.sbatch' if bounded_prefix else 'run.sbatch','run.sbatch'),
        ('mpi_particle_balance','analyze_closed.py','analyze_closed.py'),
        ('mpi_particle_balance','performance_report.py','performance_report.py'),
        ('energy_critical_control','prepare.py','prepare.py'),('space_charge_physics','launcher.py','launcher.py'),
        ('source_coupling_integration','binding.py','binding.py'),
        ('source_coupling_integration','analyze_closed.py','original_closed.py'),
        ('discrete_energy_ledger','report.py','report_discrete.py'),
        ('discrete_energy_ledger','review_discrete.py','review_discrete.py')]
    if bounded_prefix:
        install += [('mpi_particle_balance','prefix_closed.py','prefix_closed.py'),
                    ('full_scale_followup','run.py','duration_review.py')]
    for folder,name,target in install:
        shutil.copy2(cases_root/folder/name,scripts/target)
    shutil.copy2(workspace/'results/physical_electron_acceptance_20261106/weight-scan-r2/scripts/common.py',scripts/'common.py')
    for path in scripts.glob('*.py'):ast.parse(path.read_text())
    plan.update(binary=str(binary),build_record=str(record),accepted_for_execution=True,profile_interval=profile)
    plan['runner_protocol']='post-close-duration-prefix' if bounded_prefix else 'full-physical-stop'
    plan.pop('next_step',None)
    frozen_cards={};inputs={};paths=[]
    for case in plan['cases']:
        directory=Path(case['directory']);assert not (directory/'aprl.log').exists()
        assert not (directory/'reference').exists(),'Reference writer needs a nonexistent output'
        if case['mode']!='old':case.update(binary=str(binary),build_record=str(record))
        card=directory/'run.yaml';config=yaml.safe_load(card.read_text());assert case['threads']==1
        assert config['mesh']['particle_coupling']=='conservative-source'
        assert sha(card)==case['card_sha256'],'Prepared card changed'
        source=config['beam']['input']['file'];inputs[source]=case['input_sha256']
        assert sha(source)==inputs[source],'Prepared source changed'
        frozen_cards[str(card)]=sha(card);path=directory/'run-plan.json';path.write_text(json.dumps(case,indent=2)+'\n')
        paths.extend([card,path,Path(case['binary'])])
    plan['frozen_cards']=frozen_cards
    (root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    (root/'FROZEN_CARDS.json').write_text(json.dumps(frozen_cards,indent=2)+'\n')
    paths+=list(root.glob('*.json'))+[p for p in scripts.iterdir() if p.is_file()]+[Path(x) for x in inputs]
    paths+=[workspace/'cluster.env',workspace/'build/cluster-tools/binding_probe']
    frozen={str(p):sha(p) for p in paths if p.name!='FROZEN_INPUTS.json'}
    (root/'FROZEN_INPUTS.json').write_text(json.dumps(frozen,indent=2)+'\n')
    print(json.dumps(dict(prepared=True,cases=len(plan['cases']),nodes=plan.get('nodes'),
        ranks=plan.get('ranks'),threads=1,profile_interval=profile,frozen_files=len(frozen)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root');p.add_argument('--binary',required=True)
    p.add_argument('--build-record',required=True);p.add_argument('--profile',type=int,default=0)
    p.add_argument('--bounded-prefix',action='store_true')
    a=p.parse_args();main(a.root,a.binary,a.build_record,a.profile,a.bounded_prefix)
