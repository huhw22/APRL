import os,json,re,sys,time
from pathlib import Path
import numpy as np,h5py
from prepare_flow import sha,save
from launcher import launch
from binding import binding
from analyze_closed import verify,bitwise
import report_discrete
r=Path(sys.argv[1]).resolve();p=json.loads((r/'plan.json').read_text());assert os.environ.get('SLURM_JOB_ID');frozen=json.loads((r/'FROZEN.json').read_text());assert all(sha(q)==h for q,h in frozen.items());results={}
for c in p['cases']:
 d=Path(c['directory']);nodes=c['nodes'];ranks=c['ranks'];cmd=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes='+str(nodes),'--ntasks='+str(ranks),'--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];env=os.environ.copy()
 for k in ['APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_WORK_GEOMETRY','APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT','APRL_SOURCE_EXIT_AUDIT']:env.pop(k,None)
 print('START',c['name'],time.time(),flush=True);launch([*cmd,p['binding_probe'],str(ranks),'8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',nodes,ranks,8));launch([*cmd,c['binary'],str(d/'run.yaml')],d,env,c['limit_seconds']);log=(d/'aprl.log').read_text();assert all(s in log for s in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']);assert not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log)
 closed=verify(d);rows,sums,_,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5');_,_,ledger=report_discrete.calculate(rows,sums,steps);assert np.all(rows['active_macroparticles']+rows['removed_macroparticles']==p['particles'])
 for h in (d/'detectors').glob('*.h5'):
  with h5py.File(h) as f:
   if 'particle_plane' in f:
    data=f['particle_plane/records'][:];assert len(data)==len(np.unique(data['particle_id']))==p['particles'],str(h)
 for q in (d/'reference').glob('source-coupling-rank-*.csv'):
  assert q.with_suffix('.complete').read_text().split()==['1',str(steps)];a=np.genfromtxt(q,delimiter=',',names=True);assert np.all(a['self_force_defect']<1e-13) and np.all(a['relative_fixed_point_work_defect']<1e-8)
 comparison={}
 if c['method']=='linear':comparison=bitwise(r/(c['kind']+'-baseline'),d)
 if c['method']=='quintic':
  reference=r/(c['kind']+'-linear')
  for f in (d/'detectors').glob('*.h5'):
   with h5py.File(f) as h:
    if 'particle_plane' not in h:continue
    if h['particle_plane'].attrs.get('z_lab_m',1.)!=-0.05 and 'release' not in f.name:continue
    with h5py.File(reference/'detectors'/f.name) as b:
     aa=h['particle_plane/records'][:];bb=b['particle_plane/records'][:];key='particle_id';aa=np.sort(aa,order=key);bb=np.sort(bb,order=key);assert aa.dtype==bb.dtype and aa.tobytes()==bb.tobytes(),'Release state changed '+str(f)
 result=dict(program_integrity_pass=True,kind=c['kind'],method=c['method'],ranks=ranks,steps=steps,closed_scientific_files=closed,default_bitwise_reference=comparison,strict_ledger=ledger,physical_acceptance=False);save(d/'RESULT.json',result);results[c['name']]=result;save(r/'PARTIAL.json',results);print('COMPLETE',c['name'],flush=True)
assert all(sha(q)==h for q,h in frozen.items());save(r/'COMPLETE.json',dict(program_preflight_pass=True,cases=results,physical_acceptance=False));print('ALL COMPLETE',flush=True)
