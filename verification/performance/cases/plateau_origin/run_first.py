#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
import h5py
from prepare import sha,save
from launcher import launch,binding

root=Path(sys.argv[1]);plan=json.loads((root/'plan.json').read_text());d=Path(plan['case'])
frozen=json.loads((root/'FROZEN.json').read_text());assert all(sha(p)==h for p,h in frozen.items())
assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
env['APRL_PARTICLE_MOMENTUM_PROBE']=str(d/'momentum')
env['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s'])
env['APRL_PARTICLE_MOMENTUM_PROBE_FIRST']='1'
env.pop('APRL_INLET_CONSTRAINT_PROBE',None)
sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
try:
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    launch(sr+[plan['binary'],str(d/'run.yaml')],d,env,plan['timeout_s'])
    with h5py.File(d/'ledger/energy-ledger.h5') as h:
        g=h['energy_ledger'];assert int(g['complete'][()])==1
        records=g['records'][:];assert len(records)==141 and int(records['step'][-1])==140
    assert all(sha(p)==h for p,h in frozen.items())
    save(d/'COMPLETION.json',dict(checks_pass=True,physics_accepted=False,job_id=os.environ['SLURM_JOB_ID'],
        steps=140,scope='Read-only first-push trace; same PIC binary, input, geometry, inlet as guard256_500'))
except BaseException as error:
    save(d/'FAILED.json',dict(error=str(error)));raise
