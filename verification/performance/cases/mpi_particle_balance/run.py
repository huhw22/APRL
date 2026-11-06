"""Bounded single-thread MPI acceptance; inputs and current production stay frozen."""
import json
import os
from pathlib import Path
import re
import sys
import traceback

from common import require,save,sha
from launcher import launch
from binding import binding
from analyze_closed import verify_case,compare_cases


def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text())
    require(plan.get('accepted_for_execution') is True,'Campaign preflight not accepted')
    require(os.environ.get('SLURM_JOB_ID'),'Use a bounded allocation')
    frozen=json.loads((root/'FROZEN_INPUTS.json').read_text())
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen cards/scripts changed')
    build=json.loads((Path(plan['build_record'])/'COMPLETE.json').read_text())
    require(build.get('passed') is True and sha(plan['binary'])==build['binary_sha256'],'Unaccepted candidate')
    env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
    env['APRL_PROFILE_INTERVAL']=str(plan.get('profile_interval',0))
    require(env['OMP_NUM_THREADS']=='1','Only single-thread MPI is accepted in this campaign')
    results={}
    for case in plan['cases']:
        d=Path(case['directory']);nodes=case['nodes'];ranks=case['ranks']
        require(case['threads']==1 and ranks%nodes==0,'Wrong MPI-only layout')
        command=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes='+str(nodes),
            '--ntasks='+str(ranks),'--ntasks-per-node='+str(ranks//nodes),'--cpus-per-task=1',
            '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
        print('START',case['name'],flush=True)
        launch(command+[str(Path(plan['workspace'])/'build/cluster-tools/binding_probe'),str(ranks),'1'],
            d,env,120,'binding.log')
        save(d/'BINDING.json',binding(d/'binding.log',nodes,ranks,1))
        runtime=launch(command+[case['binary'],str(d/'run.yaml')],d,env,case['limit_seconds'])
        log=(d/'aprl.log').read_text()
        require(all(x in log for x in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),
            'Simulation did not close normally')
        require(not re.search(r'Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log),'Simulation error')
        save(d/'PIC_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],runtime=runtime,
            binary_sha256=sha(case['binary']),physical_acceptance=False))
        result=verify_case(case);result['runtime_seconds']=runtime['seconds']
        save(d/'RESULT.json',result);results[case['name']]=result
        print('COMPLETE',case['name'],flush=True)
    comparisons=compare_cases(plan)
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen inputs changed during execution')
    save(root/'COMPLETE.json',dict(passed=True,program_integrity_pass=True,physical_acceptance=False,
        job=os.environ['SLURM_JOB_ID'],binary_sha256=build['binary_sha256'],cases=results,comparisons=comparisons))
    print('ALL COMPLETE',flush=True)


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve()
    try:main(root)
    except Exception as error:
        save(root/'FAILED.json',dict(passed=False,error=str(error),traceback=traceback.format_exc(),
            job=os.environ.get('SLURM_JOB_ID')))
        raise
