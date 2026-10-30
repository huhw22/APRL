#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch,binding

root=Path(sys.argv[1]);kind=sys.argv[2];plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text())
require(all(sha(p)==h for p,h in frozen.items()),'Frozen file changed')
env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
env.pop('APRL_INLET_CONSTRAINT_PROBE',None)
sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
try:
    d=root/(kind+'_off');launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
    if kind=='uniform':
        te=env.copy();te['OMP_NUM_THREADS']='1';launch(['ctest','--test-dir',plan['build'],'-L','required','--output-on-failure','-j','1'],d,te,320,'required.log')
    first={};prefix={}
    dtype=np.dtype([('step','<u8'),('source_id','<u8'),('time','<f8'),('dt','<f8'),('charge','<f8'),('position','<f8',(3,)),('before','<f8',(3,)),('after','<f8',(3,)),('electric','<f8',(3,)),('magnetic','<f8',(3,))]);require(dtype.itemsize==160,'Wrong trace schema')
    for mode in ('off','on'):
        name=kind+'_'+mode;case=plan['cases'][name];d=root/name;ev=env.copy()
        ev.update(APRL_PARTICLE_MOMENTUM_PROBE=str(d/'momentum'),APRL_PARTICLE_MOMENTUM_PROBE_DT=str(plan['nominal_dt_s']),APRL_PARTICLE_MOMENTUM_PROBE_FIRST='1')
        launch(sr+[plan['binary'],str(d/'run.yaml')],d,ev,500)
        checked={}
        for sub in ('ledger','trajectory','detectors'):
            for p in (d/sub).glob('*.h5'):
                with h5py.File(p) as h,h5py.File(Path(case['comparison'])/sub/p.name) as g:
                    group=next(iter(h));a=h[group+'/records'][:];b=g[group+'/records'][:len(a)]
                    require(int(h[group+'/complete'][()])==int(g[group+'/complete'][()])==1,'Unclosed output')
                    require(a.dtype==b.dtype and a.shape==b.shape,'Prefix shape mismatch')
                    require(all(np.ascontiguousarray(a[n]).tobytes()==np.ascontiguousarray(b[n]).tobytes() for n in a.dtype.names if 'seconds' not in n),'Read-only diagnostic changed scientific prefix: '+str(p))
                    checked[sub+'/'+p.name]=len(a)
        prefix[mode]=checked
        trace=np.concatenate([np.fromfile(p,dtype=dtype) for p in sorted((d/'momentum').glob('*.bin'))])
        trace=trace[np.lexsort((trace['time'],trace['source_id']))];_,indices=np.unique(trace['source_id'],return_index=True);a=trace[indices]
        require(len(a)==32768 and np.array_equal(a['source_id'],np.arange(1,32769)),'Incomplete first-push source coverage')
        require(all(np.isfinite(a[n]).all() for n in a.dtype.names),'Non-finite trace')
        require(int(a['step'].max())<case['steps'],'First push after stop')
        first[mode]=a;np.save(d/'FIRST_SOURCE_STATES.npy',a,allow_pickle=False)
    fields=('source_id','step','time','dt','charge','position','before')
    equality={n:np.ascontiguousarray(first['off'][n]).tobytes()==np.ascontiguousarray(first['on'][n]).tobytes() for n in fields}
    require(all(equality.values()),'Actual first-push input differs')
    result=dict(pass_all=True,particles=32768,fields_bitwise=equality,first_step_range=[int(first['on']['step'].min()),int(first['on']['step'].max())],scientific_prefix_bitwise=prefix,job_id=os.environ['SLURM_JOB_ID'],scope='Earliest actual pre-force state per source ID; after/E/B intentionally excluded')
    require(all(sha(p)==h for p,h in frozen.items()),'Frozen file changed');save(root/(kind+'_SOURCE_CHECK.json'),result)
except BaseException as error:save(root/(kind+'_FAILED.json'),dict(error=str(error)));raise
