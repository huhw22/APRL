#!/usr/bin/env python3
"""Run one bounded, frozen charge case and analyze only closed outputs."""
import json,os,sys
from pathlib import Path
from prepare_scaling import sha,save
from launcher import launch,binding
from analyze import analyze

root=Path(sys.argv[1]).resolve();name=sys.argv[2];plan=json.loads((root/'plan.json').read_text());item=plan['cases'][name];d=Path(item['directory'])
frozen=json.loads((root/'FROZEN.json').read_text());assert all(sha(p)==h for p,h in frozen.items())
assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
for key in ['APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_MOMENTUM_PROBE_FIRST','APRL_PARTICLE_MOMENTUM_PROBE_STRIDE','APRL_INLET_CONSTRAINT_PROBE']:env.pop(key,None)
sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
try:
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    launch(sr+[plan['binary'],str(d/'run.yaml')],d,env,item['timeout_s'])
    result=analyze(d);assert result['execution_integrity_pass']
    assert abs(result['steps']-item['steps'])<=1
    assert all(sha(p)==h for p,h in frozen.items())
    save(d/'COMPLETION.json',dict(checks_pass=True,physics_accepted=False,job_id=os.environ['SLURM_JOB_ID'],steps=result['steps']))
except BaseException as err:
    save(d/'FAILED.json',dict(error=str(err)));raise
