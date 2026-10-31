#!/usr/bin/env python3
"""Read-only postcheck of closed cases; separate launcher guard and algebra toy.

Never changes a frozen simulation, its failed batch status or its outputs.
"""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np,yaml
import transport_review

assert os.environ.get('SLURM_JOB_ID')
out=Path(sys.argv[1]).resolve(); campaign=Path(sys.argv[2]).resolve()
plan=json.loads((campaign/'plan.json').read_text())
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def frozen_check(path):
 for p,h in json.loads(path.read_text()).items():
  assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
frozen_check(campaign/'FROZEN.json');frozen_check(out/'FROZEN.json')
results={c['name']:json.loads((Path(c['directory'])/'RESULT.json').read_text()) for c in plan['cases']}
assert all(v['integrity_passed'] for v in results.values())
comparisons={}
for layout in ['native','mpi-small']:
 for kind in ['undulator','free-drift']:
  off=campaign/f'{layout}-{kind}-off';on=campaign/f'{layout}-{kind}-on'
  pairs=[(off,on)]
  if layout=='native':pairs.insert(0,(campaign/f'{layout}-{kind}-baseline',off))
  for first,last in pairs:
   item=transport_review.compare(first,last);assert item['checks_pass'],item
   with h5py.File(first/'ledger/energy-ledger.h5') as f,h5py.File(last/'ledger/energy-ledger.h5') as g:
    x=f['energy_ledger/records'][:];y=g['energy_ledger/records'][:]
   bits={n:x[n].tobytes()==y[n].tobytes() for n in x.dtype.names if 'seconds' not in n}
   item['ledger_fields_bitwise']=bits
   for n,exact in bits.items():
    if not exact:
     assert np.all(abs(x[n]-y[n])<=512*np.finfo(float).eps*np.maximum(abs(x[n])+abs(y[n]),1e-100)),n
   comparisons[first.name+'__'+last.name]=item
save(out/'COMPARISONS.json',comparisons)
negative=out/'inconsistent-rank-setting';negative.mkdir()
card=yaml.safe_load((campaign/'native-undulator-off/run.yaml').read_text())
card['output']['manifest']=str(negative/'run-manifest.yaml')
card['mesh']['moving_window']['reference_directory']=str(negative/'reference')
for key,folder in [('detectors','detectors'),('trajectory','trajectory'),('energy_ledger','ledger')]:card[key]['directory']=str(negative/folder)
(negative/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
wrapper=negative/'rank_env.py'
wrapper.write_text('import os,sys\nrank=int(os.environ["SLURM_PROCID"])\nif rank==0:os.environ["APRL_PARTICLE_WORK_GEOMETRY"]=sys.argv[3]\nelse:os.environ.pop("APRL_PARTICLE_WORK_GEOMETRY",None)\nos.execv(sys.argv[1],[sys.argv[1],sys.argv[2]])\n')
candidate=next(c['binary'] for c in plan['cases'] if c['mode']=='on');start=time.monotonic()
with (negative/'aprl.log').open('x') as f:
 r=subprocess.run(['timeout','--signal=TERM','--kill-after=5','30','srun','--mpi=pmix_v3','--exact','--ntasks=2','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1',sys.executable,str(wrapper),candidate,str(negative/'run.yaml'),str(negative/'geometry')],stdout=f,stderr=subprocess.STDOUT,timeout=40)
assert r.returncode not in [0,124,137]
assert 'identical valid settings on all ranks' in (negative/'aprl.log').read_text()
save(negative/'RESULT.json',{'negative_control_passed':True,'seconds':time.monotonic()-start,'returncode':r.returncode,'rank_source':'SLURM_PROCID available before MPI_Init','physical_simulation_started':False})
with (out/'linear-response.log').open('x') as f:
 subprocess.run([sys.executable,str(out/'scripts/linear_response_probe.py'),str(out/'linear-response')],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=90)
frozen_check(campaign/'FROZEN.json');frozen_check(out/'FROZEN.json')
save(out/'COMPLETE.json',{'checks_passed':True,'job':os.environ['SLURM_JOB_ID'],'simulation_campaign':str(campaign),'closed_scientific_cases':len(results),'physical_state_bitwise_pass':True,'original_campaign_batch_status_unchanged':'FAILED at launcher-wrapper negative control after all science cases; no simulation rerun','physics_method_changes':False,'comparisons':comparisons,'cases':results})
print('READ-ONLY POSTCHECK, LAUNCHER GUARD AND CLASSICAL TOY COMPLETE',flush=True)
