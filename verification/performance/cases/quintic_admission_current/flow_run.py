"""Run sealed old/new linear/quintic short cards inside one Slurm allocation."""
import json
import os
import re
import sys
from pathlib import Path

from flow_prepare import require, save, sha
from launcher import launch
from binding import binding
from flow_analyze import analyze_case, compare


def main(root, candidate_record=None):
    root = Path(root).resolve()
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Run inside a bounded Slurm allocation')
    plan = json.loads((root / 'plan.json').read_text())
    candidate_record = Path(candidate_record or (Path(plan['candidate_build_record']) / 'COMPLETE.json')).resolve()
    frozen = json.loads((root / 'FROZEN_INPUTS.json').read_text())
    require(all(sha(p) == h for p, h in frozen.items()), 'Prepared inputs changed')
    identity = json.loads(candidate_record.read_text())
    require(identity.get('passed') is True, 'Candidate build was not accepted')
    candidate = Path(plan['candidate_binary'])
    require(Path(identity['binary']).resolve() == candidate, 'Wrong candidate program')
    require(sha(candidate) == identity['binary_sha256'], 'Candidate binary identity differs')
    frozen[str(candidate)] = identity['binary_sha256']
    frozen[str(candidate_record)] = sha(candidate_record)
    save(root / 'FROZEN.json', frozen)
    results = {}
    for case in plan['cases']:
        directory = Path(case['directory'])
        command = ['srun', '--mpi=pmix_v3', '--exclusive', '--exact', '--nodes=1',
                   '--ntasks=16', '--ntasks-per-node=16', '--cpus-per-task=8',
                   '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']
        env = os.environ.copy()
        for key in ['APRL_PARTICLE_MOMENTUM_PROBE', 'APRL_PARTICLE_WORK_GEOMETRY',
                    'APRL_GRID_WORK_AUDIT', 'APRL_FIELD_STAGE_AUDIT', 'APRL_SOURCE_EXIT_AUDIT']:
            env.pop(key, None)
        print('START', case['name'], flush=True)
        launch(command + [plan['binding_probe'], '16', '8'], directory, env, 60, 'binding.log')
        save(directory / 'BINDING.json', binding(directory / 'binding.log', 1, 16, 8))
        runtime = launch(command + [case['binary'], str(directory / 'run.yaml')],
                         directory, env, case['limit_seconds'])
        log = (directory / 'aprl.log').read_text()
        require(all(x in log for x in ['Configured stop reached:', 'phase=solve-complete',
                                       'phase=mpi-finalize-complete']), 'Incomplete normal stop')
        require(not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT', log),
                'Error in short-flow log')
        result = analyze_case(directory, plan['particles'], plan['expected_steps'])
        result['runtime_seconds'] = runtime['seconds']
        save(directory / 'RESULT.json', result)
        results[case['name']] = result
        save(root / 'PARTIAL.json', results)
        print('COMPLETE', case['name'], flush=True)
    comparisons = compare(root)
    require(all(sha(p) == h for p, h in frozen.items()), 'Frozen program or inputs changed')
    save(root / 'COMPLETE.json', dict(program_integrity_pass=True, physical_acceptance=False,
                                     cases=results, comparison=comparisons))
    print('ALL COMPLETE', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
