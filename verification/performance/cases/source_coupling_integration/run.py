#!/usr/bin/env python3
import json,os,re,sys,subprocess,time
from pathlib import Path
import numpy as np
import h5py
from prepare import sha,save,require
from launcher import launch
from binding import binding
from analyze_closed import verify,bitwise
import report_discrete
root=Path(sys.argv[1]).resolve();plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text());assert os.environ.get('SLURM_JOB_ID');assert all(sha(p)==h for p,h in frozen.items())
results={};env=os.environ.copy()
for key in ['APRL_PARTICLE_WORK_GEOMETRY','APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT','APRL_INLET_CONSTRAINT_PROBE','APRL_PARTICLE_MOMENTUM_PROBE']:env.pop(key,None)
for name,case in plan['cases'].items():
 if case.get('reuse_closed'):
  result=json.loads((Path(case['reuse_closed'])/'RESULT.json').read_text());assert result['program_integrity_pass'];verify(case['reuse_closed']);results[name]=dict(result,reused_closed=case['reuse_closed']);print('REUSED',name,flush=True);continue
 d=Path(case['directory']);ranks=case['ranks'];nodes=plan.get('nodes',1);sr=['srun','--mpi=pmix_v3','--exact','--nodes='+str(nodes),'--ntasks='+str(ranks),'--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
 if nodes==2:sr+=['--ntasks-per-node=8']
 launch([*sr,plan['binding_probe'],str(ranks),'8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',nodes,ranks,8))
 launch([*sr,plan['binary'],str(d/'run.yaml')],d,env,case['limit_seconds'])
 log=(d/'aprl.log').read_text();assert all(q in log for q in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']);assert not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log)
 closed=verify(d);rows,sums,count,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');_,_,ledger=report_discrete.calculate(rows,sums,steps)
 assert np.all(rows['active_macroparticles']+rows['removed_macroparticles']==plan['particles'])
 for file in (d/'detectors').glob('*.h5'):
  with h5py.File(file) as h:
   if 'particle_plane' in h:
    data=h['particle_plane/records'][:];assert len(data)==len(np.unique(data['particle_id']))==plan['particles'],str(file)
 if case['mode']=='off':comparison=bitwise(Path(case['reference']),d)
 else:
  comparison={};paths=list((d/'reference').glob('source-coupling-rank-*.csv'));assert len(paths)==ranks
  for p in paths:
   assert p.with_suffix('.complete').read_text().split()==['1',str(steps)];x=np.genfromtxt(p,delimiter=',',names=True);assert np.array_equal(x['step'],np.arange(1,steps+1));assert np.all(x['self_force_defect']<1e-13) and np.all(x['relative_fixed_point_work_defect']<1e-8);assert np.all(x['iterations']<=8)
 if case.get('same_method_reference'):
  reference=Path(case['same_method_reference'])
  if (reference/'RESULT.json').exists():comparison=bitwise(reference,d)
  else:comparison={'pending_closed_reference':str(reference)}
 result={'program_integrity_pass':True,'mode':case['mode'],'ranks':ranks,'steps':steps,'closed_scientific_files':closed,'default_or_same_method_bitwise_reference':comparison,'strict_ledger':ledger,'full_sphere_physics_accepted':False};save(d/'RESULT.json',result);results[name]=result;save(root/'PARTIAL.json',results);print('COMPLETE',name,steps,flush=True)
if plan.get('rank_mode_negative_configs'):
 nodes=plan.get('nodes',1);command=['srun','--mpi=pmix_v3','--exact','--nodes='+str(nodes),'--ntasks=2','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=10']
 if nodes==2:command+=['--ntasks-per-node=1']
 command += [sys.executable,str(root/'scripts/select_rank_mode.py'),plan['binary'],*plan['rank_mode_negative_configs']]
 started=time.monotonic()
 with (root/'inconsistent-mode/aprl.log').open('x') as f:q=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=40)
 log=(root/'inconsistent-mode/aprl.log').read_text();assert q.returncode!=0 and 'particle_coupling must agree on every MPI rank' in log
 save(root/'inconsistent-mode/RESULT.json',{'passed':True,'returncode':q.returncode,'seconds':time.monotonic()-started,'rejected_before_initialization':True})
assert all(sha(p)==h for p,h in frozen.items());save(root/'COMPLETE.json',{'program_preflight_pass':True,'job':os.environ['SLURM_JOB_ID'],'cases':results,'full_sphere_physics_accepted':False});print('ALL COMPLETE',flush=True)
