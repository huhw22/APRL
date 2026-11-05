"""Run one bounded control case inside a one-node Slurm allocation."""
import argparse
import json
import os
from pathlib import Path
import re

from prepare import require, save, sha
from launcher import launch
from binding import binding
from analyze import analyze_case


def run(directory, allow_full=False):
    directory = Path(directory).resolve()
    root = directory.parent
    plan = json.loads((directory/'run-plan.json').read_text())
    campaign_plan = json.loads((root/'plan.json').read_text())
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Use a bounded Slurm allocation')
    require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == 1, 'Control case requires one node')
    if plan['stage'] == 'full':
        require(allow_full, 'Explicitly enable the full short control after build and smoke acceptance')
        smoke_record = Path(campaign_plan['smoke_record'])
        smoke = json.loads(smoke_record.read_text())
        require(smoke.get('program_integrity_pass') is True, 'Smoke regression not accepted')
        smoke_identity = json.loads((smoke_record.parent/'FROZEN.json').read_text())
    frozen = json.loads((root/'FROZEN_INPUTS.json').read_text())
    require(all(sha(path) == digest for path,digest in frozen.items()), 'Prepared cards or scripts changed')
    identity_path = Path(plan['build_record'])/'COMPLETE.json'
    identity = json.loads(identity_path.read_text())
    require(identity.get('passed') is True, 'Candidate program has not passed build/module validation')
    binary = Path(plan['binary'])
    require(Path(identity['binary']).resolve() == binary, 'Candidate binary path differs')
    require(sha(binary) == identity['binary_sha256'], 'Candidate binary identity differs')
    if plan['stage'] == 'full':
        require(smoke_identity.get(str(binary)) == identity['binary_sha256'],
                'Smoke did not validate this candidate binary')
        frozen[str(smoke_record)] = sha(smoke_record)
        frozen[str(smoke_record.parent/'FROZEN.json')] = sha(smoke_record.parent/'FROZEN.json')
    frozen[str(binary)] = identity['binary_sha256']
    frozen[str(identity_path)] = sha(identity_path)
    save(directory/'FROZEN.json', frozen)
    env = os.environ.copy()
    for key in ['APRL_PARTICLE_MOMENTUM_PROBE', 'APRL_PARTICLE_MOMENTUM_PROBE_FIRST',
                'APRL_PARTICLE_MOMENTUM_PROBE_STRIDE', 'APRL_PARTICLE_MOMENTUM_PROBE_DT',
                'APRL_INLET_CONSTRAINT_PROBE', 'APRL_GRID_WORK_AUDIT', 'APRL_FIELD_STAGE_AUDIT',
                'APRL_PARTICLE_WORK_GEOMETRY', 'APRL_SOURCE_EXIT_AUDIT',
                'APRL_DIAGNOSTIC_STEPS', 'APRL_FORCE_SPLIT_DIRECTORY']:
        env.pop(key, None)
    command = ['srun', '--mpi=pmix_v3', '--exclusive', '--exact', '--nodes=1',
        '--ntasks=16', '--ntasks-per-node=16', '--cpus-per-task=8', '--cpu-bind=cores',
        '--kill-on-bad-exit=1', '--wait=30']
    workspace = Path(campaign_plan['workspace'])
    launch(command+[str(workspace/'build/cluster-tools/binding_probe'), '16', '8'],
           directory, env, 120, 'binding.log')
    save(directory/'BINDING.json', binding(directory/'binding.log', 1, 16, 8))
    runtime = launch(command+[str(binary), str(directory/'run.yaml')], directory, env, plan['case_limit_seconds'])
    log = (directory/'aprl.log').read_text()
    require(all(s in log for s in ['Configured stop reached:', 'phase=solve-complete', 'phase=mpi-finalize-complete']), 'Incomplete normal stop')
    require(not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT', log), 'Error in control log')
    result = analyze_case(directory, plan)
    result['runtime_seconds'] = runtime['seconds']
    result['job'] = os.environ['SLURM_JOB_ID']
    result['binary_sha256'] = identity['binary_sha256']
    require(all(sha(path) == digest for path,digest in frozen.items()), 'Frozen identity changed while running')
    save(directory/'RESULT.json', result)
    save(directory/'PIC_COMPLETE.json', dict(run=runtime, analysis_complete=True, stage=plan['stage']))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    p.add_argument('--allow-full', action='store_true')
    a = p.parse_args()
    run(a.directory, a.allow_full)
