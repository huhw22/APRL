"""Frozen inputs, bounded MPI launches, and closed-output checks."""
import hashlib,json,os,re,subprocess,time
from pathlib import Path
import h5py,yaml
import numpy as np

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def save(path,value):
 p=Path(path);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n');q.replace(p)
def duration_steps(card):
 assert card['mesh']['field_solver']=='cowan-z'
 dt=card['mesh']['cell_size'][2]/299792458.
 assert dt>0 and card['mesh']['duration']>0
 t=0.;steps=0
 while t<card['mesh']['duration']:
  t+=dt;steps+=1
  assert steps<10000000,'Unbounded diagnostic duration'
 return steps,t,dt
def completed_steps(case,log,physical):
 card=yaml.safe_load((Path(case['directory'])/'run.yaml').read_text())
 budget,final_time,dt=duration_steps(card)
 assert budget==case['duration_step_budget'],(budget,case['duration_step_budget'])
 matches=re.findall(r'^Moving window completed: steps=(\d+)\b',log,re.M)
 assert len(matches)==1,'Missing or ambiguous completed-step report'
 steps=int(matches[0]);assert 0<steps<=budget,(steps,budget)
 if physical:
  assert 'Configured stop reached:' in log
  assert 'mesh.duration was exhausted before the configured physical stop' not in log
 else:
  assert steps==budget,(steps,budget)
  assert 'Configured stop reached:' not in log
 return steps,dt
def frozen(root):
 assert os.environ.get('SLURM_JOB_ID'),'Run on Slurm'
 for path,digest in json.loads((root/'FROZEN.json').read_text()).items():assert sha(path)==digest,path
 return json.loads((root/'PLAN.json').read_text())
def environment(directory,observer=True):
 env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
 env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PLACES='cores',OMP_PROC_BIND='close',
  OMPI_MCA_btl='self,vader,tcp',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',HDF5_USE_FILE_LOCKING='FALSE')
 if observer:env.update(APRL_LAB_WORLDTUBE=str(directory/'worldtube'),APRL_LAB_FIELD_AUDIT=str(directory/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(directory/'lab-audit'))
 return env
def srun():
 return ['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
def binding(root,plan,filename):
 p=root/filename
 with p.open('x') as f:subprocess.run(srun()+[plan['binding_probe'],'16','8'],env=environment(root,False),stdout=f,stderr=subprocess.STDOUT,check=True,timeout=90)
 text=p.read_text();rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)',text)
 assert '[binding] result=PASS failures=0' in text and len(rows)==128
 assert {(int(r[0]),int(r[1])) for r in rows}=={(r,t) for r in range(16) for t in range(8)}
 assert all(r[-1]=='PASS' for r in rows) and len({(r[2],r[4],r[5]) for r in rows})==128
 return dict(physical_cores=128,mpi=16,threads=8,hosts=sorted({r[2] for r in rows}))
def launch(plan,case,observer=True):
 d=Path(case['directory']);assert not (d/'aprl.log').exists(),'Refuse replacing attempted output'
 limit=case['timeout_seconds'];start=time.monotonic()
 command=['timeout','--signal=TERM','--kill-after=30',str(limit)]+srun()+[plan['binary'],str(d/'run.yaml')]
 save(d/'STARTED.json',dict(job_id=os.environ['SLURM_JOB_ID'],started_unix=time.time(),command=command,observer=observer,self_field_feedback=True))
 with (d/'aprl.log').open('x') as f:
  result=subprocess.run(command,cwd=d,env=environment(d,observer),stdout=f,stderr=subprocess.STDOUT,timeout=limit+90)
 record=dict(exit_code=result.returncode,seconds=time.monotonic()-start,job_id=os.environ['SLURM_JOB_ID'])
 save(d/'LAUNCH.json',record);return record
def validate(case,physical=True,observer=True):
 d=Path(case['directory']);log=(d/'aprl.log').read_text();steps,dt=completed_steps(case,log,physical)
 assert not re.search('Segmentation fault|OpenFabrics device|DUE TO TIME LIMIT',log)
 if physical:
  assert 'Configured stop reached:' in log and 'phase=solve-complete' in log and 'phase=mpi-finalize-complete' in log
  assert 'Simulation failure' not in log
 else:
  assert 'mesh.duration was exhausted before the configured physical stop' in log
  assert 'phase=duration-exhaustion-readiness-sync-complete' in log
 assert 'self_field_feedback is disabled' not in log and '[one-way-source]' not in log
 with h5py.File(d/'ledger/energy-ledger.h5') as f:
  g=f['energy_ledger'];assert int(g['complete'][()])==int(physical)
  assert int(g['committed_records'][()])==len(g['records'])==steps+1
  assert int(g['records'][-1]['step'])==steps
  assert np.array_equal(g['records']['step'],np.arange(steps+1))
  assert np.all(np.diff(g['records']['time_box_s'])>0)
  assert abs(float(g['records'][-1]['time_box_s'])-steps*dt)<1e-10*dt*max(1,steps)
 samples={};particle_records={}
 for path in sorted((d/'detectors').glob('*.h5')):
  with h5py.File(path) as f:
   g=f[next(iter(f))];assert int(g['complete'][()])==int(physical),path
   if 'committed_samples' in g:
    n=int(g['committed_samples'][()]);assert n==len(g['time_s'])==len(g['fields']);samples[path.name]=n
   elif 'records' in g:
    n=int(g['committed_records'][()]);assert n==len(g['records']),path
    if 'particle_plane' in f or 'particle_substep_plane' in f:
     particle_records[path.name]=n
     if physical:assert n==case['particle_count'],(path,n)
    elif physical:
     a=g['records'][-1];assert int(a['particle_downstream_count'])==case['particle_count'] and int(a['particle_upstream_count'])==0,path
 if observer:
  for marker in ['[lab-worldtube]','[lab-particle-sides]','[lab-field-audit]','[discrete-energy-ledger]']:assert marker in log,marker
  assert (d/'worldtube/COMPLETE').read_text().split()==[str(int(physical)),str(steps)]
  meta=json.loads((d/'worldtube/METADATA.json').read_text())
  assert (d/'worldtube/field_lines.bin').stat().st_size==steps*meta['nz']*6*8
  for rank in range(16):
   assert (d/'worldtube'/f'particle-sides-rank-{rank:05d}.csv').is_file()
   for prefix in ['lab-field','paired-momentum','whole-cowan']:
    assert (d/'lab-audit'/f'{prefix}-rank-{rank:05d}.complete').read_text().split()==[str(int(physical)),str(steps)]
   assert (d/'ledger'/f'discrete-field-rank-{rank:05d}.complete').read_text().split()==[str(int(physical)),str(steps)]
 return dict(complete=physical,actual_steps=steps,duration_step_budget=case['duration_step_budget'],
  legacy_expected_steps=case.get('legacy_expected_steps'),step_gate='actual physical stop plus mutually consistent closed records',
  all_outputs_closed=True,field_samples=samples,particle_plane_records=particle_records,physical_precision_accepted=False)
