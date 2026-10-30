#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save
from launcher import launch,binding
from analyze_momentum import analyze
import report_discrete

root=Path(sys.argv[1]);name=sys.argv[2];plan=json.loads((root/'plan.json').read_text());case=plan['cases'][name];d=Path(case['directory'])
frozen=json.loads((root/'FROZEN.json').read_text());assert all(sha(p)==h for p,h in frozen.items())
assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
env['APRL_PARTICLE_MOMENTUM_PROBE']=str(d/'momentum');env['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s']);env.pop('APRL_PARTICLE_MOMENTUM_PROBE_FIRST',None)
env.pop('APRL_INLET_CONSTRAINT_PROBE',None)
sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
try:
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    launch(sr+[plan['binary'],str(d/'run.yaml')],d,env,case.get('timeout_s',900))
    a=analyze(d,Path(case['original']),d/'momentum-analysis','converging');assert a['steps']==case['steps'] and a['all_released_step']==123 and a['first_exit_step'] is None
    with h5py.File(d/'detectors/entry.h5') as h,h5py.File(Path(case['comparison'])/'detectors/entry.h5') as g:
        x=h['particle_plane'];y=g['particle_plane'];assert int(x['complete'][()])==int(y['complete'][()])==1
        p=x['records'][:];q=y['records'][:];assert len(p)==len(q)==32768 and p.dtype==q.dtype
        p=p[np.argsort(p['source_id'])];q=q[np.argsort(q['source_id'])]
        changes={n:dict(components=int(np.count_nonzero(p[n]!=q[n])),maximum_absolute_difference=float(np.max(np.abs(p[n]-q[n])))) for n in p.dtype.names if p[n].tobytes()!=q[n].tobytes()}
        assert set(changes)<=set(['proper_velocity'])
        save(d/'ENTRY_DIAGNOSTIC.json',dict(records=32768,changed_fields=changes,all_other_fields_bitwise=True))
    r,s,_,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');_,_,summary=report_discrete.calculate(r,s,steps)
    save(d/'LEDGER_CHECK.json',summary)
    assert all(sha(p)==h for p,h in frozen.items());save(d/'COMPLETION.json',dict(checks_pass=True,physics_accepted=False,job_id=os.environ['SLURM_JOB_ID'],steps=case['steps'],scope='Same physical source and core; boundary placement or stop only; no production acceptance'))
except BaseException as error:save(d/'FAILED.json',dict(error=str(error)));raise
