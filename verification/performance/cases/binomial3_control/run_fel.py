"""Bounded two-node execution, gated on completed three-point acceptance."""
import json
import os
from pathlib import Path
import re
import sys
import traceback

from common import require,save,sha
from binding import binding
from launcher import launch


def run(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text())
    case=plan['cases'][0];directory=Path(case['directory'])
    require(plan['accepted_for_execution'],'Card not accepted for execution')
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ['SLURM_JOB_NUM_NODES'])==case['nodes'],
        'Wrong/missing bounded allocation')
    frozen=json.loads((root/'FROZEN_INPUTS.json').read_text())
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen inputs changed')
    identity=json.loads((Path(plan['build_record'])/'COMPLETE.json').read_text())
    acceptance=json.loads((Path(plan['acceptance_record'])/'COMPLETE.json').read_text())
    review=json.loads((Path(plan['acceptance_record'])/'INDEPENDENT_REVIEW.json').read_text())
    require(identity.get('passed') and acceptance.get('program_integrity_pass'),
        'Module or full-flow acceptance missing')
    require(review.get('passed') is True,'Independent convergence/output review missing')
    require(identity['binary_sha256']==acceptance['binary_sha256']==sha(plan['binary']),
        'Binary differs from acceptance')
    require(not (directory/'reference').exists(),'Reference output directory already exists')
    command=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes='+str(case['nodes']),
        '--ntasks='+str(case['ranks']),'--ntasks-per-node=16','--cpus-per-task=8',
        '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
    env['APRL_PROFILE_INTERVAL']='0'
    probe=str(Path(plan['workspace'])/'build/cluster-tools/binding_probe')
    launch(command+[probe,str(case['ranks']),str(case['threads'])],directory,env,120,'binding.log')
    save(directory/'BINDING.json',binding(directory/'binding.log',case['nodes'],case['ranks'],case['threads']))
    runtime=launch(command+[plan['binary'],str(directory/'run.yaml')],directory,env,
        case['case_limit_seconds'],'aprl.log')
    log=(directory/'aprl.log').read_text()
    require(all(token in log for token in ['Configured stop reached:','phase=solve-complete',
        'phase=mpi-finalize-complete','linear-binomial3']), 'No complete three-point normal stop')
    require(not re.search(r'Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log),
        'Error in simulation log')
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen inputs changed during run')
    save(directory/'PIC_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],
        binary_sha256=identity['binary_sha256'],runtime=runtime,physical_acceptance=False,
        analysis_pending=True))
    print('PIC COMPLETE; closed-file analysis runs in its separate dependent job.',flush=True)


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve()
    try:run(root)
    except Exception as error:
        save(root/'one4one/FAILED.json',dict(complete=False,error=str(error),
            traceback=traceback.format_exc(),job=os.environ.get('SLURM_JOB_ID')))
        raise
