#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch,binding
from analyze_momentum import analyze
import report_discrete

root=Path(sys.argv[1]);name=sys.argv[2];plan=json.loads((root/'plan.json').read_text());case=plan['cases'][name];d=Path(case['directory'])
frozen=json.loads((root/'FROZEN.json').read_text())
assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
try:
    assert all(sha(p)==h for p,h in frozen.items())
    env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    env['APRL_PARTICLE_MOMENTUM_PROBE']=str(d/'momentum');env['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s'])
    sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    launch(sr+[plan['binary'],str(d/'run.yaml')],d,env,case['timeout_s'])
    mo=analyze(d,Path(plan['original']),d/'analysis',case['kind']);assert mo['steps']==case['steps'] and mo['first_exit_step'] is None
    r,s,_,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');_,_,summary=report_discrete.calculate(r,s,steps);save(d/'LEDGER_CHECK.json',summary)
    if name=='base':
        with h5py.File(d/'ledger/energy-ledger.h5') as h,h5py.File(Path(plan['prior'])/'converging_full/ledger/energy-ledger.h5') as g:
            a=h['energy_ledger/records'][:];b=g['energy_ledger/records'][:len(a)]
            assert a.dtype==b.dtype and a.shape==b.shape
            assert all(np.ascontiguousarray(a[n]).tobytes()==np.ascontiguousarray(b[n]).tobytes() for n in a.dtype.names if 'seconds' not in n)
        save(d/'PREFIX_BITWISE.json',dict(pass_all=True,rows=len(a),scope='all non-timing ledger fields vs prior full converging prefix'))
    assert all(sha(p)==h for p,h in frozen.items())
    save(d/'COMPLETION.json',dict(checks_pass=True,job_id=os.environ['SLURM_JOB_ID'],physics_accepted=False))
except BaseException as error:
    save(d/'FAILED.json',dict(error=str(error)));raise
