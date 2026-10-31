#!/usr/bin/env python3
import json,os,re,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch
from binding import binding
from analyze_closed import verify
from analyze_particles import analyze
import report_discrete
root=Path(sys.argv[1]).resolve();group=sys.argv[2];plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text());assert os.environ.get('SLURM_JOB_ID');assert all(sha(p)==h for p,h in frozen.items());env=os.environ.copy()
for key in ['APRL_PARTICLE_WORK_GEOMETRY','APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT','APRL_INLET_CONSTRAINT_PROBE']:env.pop(key,None)
dtype=np.dtype([('step','<u8'),('source_id','<u8'),('time','<f8'),('dt','<f8'),('charge','<f8'),('position','<f8',(3,)),('before','<f8',(3,)),('after','<f8',(3,)),('electric','<f8',(3,)),('magnetic','<f8',(3,))]);assert dtype.itemsize==160
for name in plan['jobs'][group]:
 case=plan['cases'][name];d=Path(case['directory']);sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];launch([*sr,plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1));ev=env.copy();ev.update(APRL_PARTICLE_MOMENTUM_PROBE=str(d/'momentum'),APRL_PARTICLE_MOMENTUM_PROBE_DT=str(plan['nominal_dt_s']),APRL_PARTICLE_MOMENTUM_PROBE_FIRST='1')
 launch([*sr,plan['binary'],str(d/'run.yaml')],d,ev,case['timeout_s']);log=(d/'aprl.log').read_text();assert all(x in log for x in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']);assert not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log)
 files=verify(d);physical=analyze(d);assert physical['execution_integrity_pass'];rows,sums,ranks,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');assert steps==case['steps'];_,_,ledger=report_discrete.calculate(rows,sums,steps)
 paths=list((d/'reference').glob('source-coupling-rank-*.csv'));assert len(paths)==16
 for p in paths:
  assert p.with_suffix('.complete').read_text().split()==['1',str(steps)];x=np.genfromtxt(p,delimiter=',',names=True);assert np.array_equal(x['step'],np.arange(1,steps+1));assert np.all(x['self_force_defect']<1e-13) and np.all(x['relative_fixed_point_work_defect']<1e-8)
 traces=[np.fromfile(p,dtype=dtype) for p in sorted((d/'momentum').glob('*.bin'))];assert traces;first=np.concatenate(traces);first=first[np.lexsort((first['time'],first['source_id']))];_,index=np.unique(first['source_id'],return_index=True);first=first[index];assert len(first)==32768;np.save(d/'FIRST_SOURCE_STATES.npy',first,allow_pickle=False)
 source_comparison={}
 if case['first_source_reference']:
  prior=np.load(case['first_source_reference']);source_comparison={n:np.ascontiguousarray(first[n]).tobytes()==np.ascontiguousarray(prior[n]).tobytes() for n in ['step','source_id','time','dt','charge','position','before']};assert all(source_comparison.values()),source_comparison
 save(d/'RESULT.json',{'execution_integrity_pass':True,'scientific_files':files,'original_first_force_state_bitwise':source_comparison,'strict_ledger':ledger,'physics_summary':physical,'full_physics_accepted':False,'job':os.environ['SLURM_JOB_ID']});print('COMPLETE',name,flush=True)
assert all(sha(p)==h for p,h in frozen.items());save(root/(group+'_COMPLETE.json'),{'execution_integrity_pass':True,'cases':plan['jobs'][group],'job':os.environ['SLURM_JOB_ID'],'physics_accepted':False})
