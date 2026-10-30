#!/usr/bin/env python3
import os,json,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch,binding
from analyze_momentum import analyze as momentum
from analyze_particles import analyze as particles
import report_discrete

root=Path(sys.argv[1]);name=sys.argv[2];plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text());d=root/name;d.mkdir(exist_ok=True)
try:
    assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8';assert all(sha(p)==h for p,h in frozen.items())
    env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    for k in ('APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_MOMENTUM_PROBE_DT','APRL_INLET_CONSTRAINT_PROBE'):env.pop(k,None)
    sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    if name=='preflight':
        te=env.copy();te['OMP_NUM_THREADS']='1';launch(['ctest','--test-dir',plan['build'],'-L','required','--output-on-failure','-j','1'],d,te,320,'required.log')
        case=plan['cases']['default_100'];p=Path(case['directory']);launch(sr+[plan['binary'],str(p/'run.yaml')],p,env,300)
        checked={}
        for sub in ('ledger','trajectory','detectors'):
            for f in (Path(case['comparison'])/sub).glob('*.h5'):
                with h5py.File(f) as h,h5py.File(p/sub/f.name) as g:
                    group=next(iter(h));x=h[group+'/records'][:];y=g[group+'/records'][:];assert int(g[group+'/complete'][()])==1 and x.dtype==y.dtype and x.shape==y.shape
                    assert all(np.ascontiguousarray(x[n]).tobytes()==np.ascontiguousarray(y[n]).tobytes() for n in x.dtype.names if 'seconds' not in n)
                    checked[sub+'/'+f.name]=len(x)
        save(d/'DEFAULT_BITWISE.json',checked)
    else:
        assert json.loads((root/'preflight/COMPLETION.json').read_text())['checks_pass'];case=plan['cases'][name];p=Path(case['directory']);ev=env.copy()
        if case['mode']=='pilot':ev['APRL_PARTICLE_MOMENTUM_PROBE']=str(p/'momentum');ev['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s']);ev['APRL_INLET_CONSTRAINT_PROBE']=str(p/'diagnostic')
        launch(sr+[plan['binary'],str(p/'run.yaml')],p,ev,case['timeout_s'])
        mo=momentum(p,Path(case['original']),p/'momentum-analysis',case['kind']);assert mo['steps']==case['steps']
        if case['mode']=='full':assert particles(p)['execution_integrity_pass']
        # Trusted plane is an immutable phase-space boundary condition.
        with h5py.File(p/'detectors/entry.h5') as h,h5py.File(Path(case['comparison'])/'detectors/entry.h5') as g:
            a=h['particle_plane/records'][:];b=g['particle_plane/records'][:];assert len(a)==len(b)==32768
            a=a[np.argsort(a['source_id'])];b=b[np.argsort(b['source_id'])]
            assert all(np.ascontiguousarray(a[n]).tobytes()==np.ascontiguousarray(b[n]).tobytes() for n in a.dtype.names)
        save(p/'ENTRY_PHASE_SPACE_BITWISE.json',dict(pass_all=True,particles=len(a),all_record_fields=True))
        r,s,_,steps=report_discrete.load(p/'ledger',p/'ledger/energy-ledger.h5');_,_,summary=report_discrete.calculate(r,s,steps);save(p/'LEDGER_CHECK.json',summary)
    assert all(sha(p)==h for p,h in frozen.items());save(d/'COMPLETION.json',dict(checks_pass=True,physics_accepted=False,job_id=os.environ['SLURM_JOB_ID']))
except BaseException as error:save(d/'FAILED.json',dict(error=str(error)));raise
