#!/usr/bin/env python3
"""Same 16x8 decomposition across two hosts; compare closed one-host prefix."""
import json,os,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch,binding

root=Path(sys.argv[1]);d=root/'two_nodes_200';plan=json.loads((root/'plan.json').read_text())
frozen=json.loads((root/'FROZEN.json').read_text())
extra=json.loads((d/'FROZEN.json').read_text());frozen.update(extra)
require(all(sha(p)==h for p,h in frozen.items()),'Frozen input changed')
require(json.loads((root/'surface_500/COMPLETION.json').read_text())['checks_pass'],'Pilot incomplete')
require(os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8' and os.environ['SLURM_JOB_NUM_NODES']=='2','Wrong allocation')
env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
for k in ('APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_MOMENTUM_PROBE_DT'):env.pop(k,None)
sr=['srun','--mpi=pmix_v3','--exact','--nodes=2','--ntasks=16','--ntasks-per-node=8','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
try:
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log')
    save(d/'BINDING.json',binding(d/'binding.log',2))
    launch(sr+[plan['binary'],str(d/'run.yaml')],d,env,600)
    checked={}
    for sub in ('ledger','trajectory','detectors'):
        for p in (d/sub).glob('*.h5'):
            with h5py.File(p) as h,h5py.File(root/'surface_500'/sub/p.name) as old:
                g=next(iter(h));x=h[g+'/records'][:];y=old[g+'/records'][:len(x)]
                require(int(h[g+'/complete'][()])==1 and int(old[g+'/complete'][()])==1,'Incomplete H5')
                require(x.shape==y.shape and x.dtype==y.dtype,'Prefix shape/schema mismatch')
                require(all(np.ascontiguousarray(x[n]).tobytes()==np.ascontiguousarray(y[n]).tobytes() for n in x.dtype.names if 'seconds' not in n),'Scientific prefix mismatch: '+str(p))
                checked[sub+'/'+p.name]=len(x)
    require(all(sha(p)==h for p,h in frozen.items()),'Frozen file changed')
    save(d/'COMPLETION.json',dict(checks_pass=True,job_id=os.environ['SLURM_JOB_ID'],steps=200,files_bitwise=checked,scope='same 16 MPI x8; two hosts vs one host, 200-step scientific prefix'))
except BaseException as error:
    save(d/'FAILED.json',dict(error=str(error)));raise
