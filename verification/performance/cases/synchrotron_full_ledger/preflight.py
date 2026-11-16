"""Short full-input diagnostic invariance gate, deliberately not a physical run."""
import json,os,sys
from pathlib import Path
import h5py,numpy as np
from common import frozen,binding,launch,validate,save,sha

R=Path(sys.argv[1]).resolve();P=frozen(R);bound=binding(R,P,'preflight-binding.log');results=[]
for case in P['preflight']:
 run=launch(P,case,case['observer']);assert run['exit_code'] not in [0,124,137]
 result=validate(case,False,case['observer']);results.append(dict(case=case['name'],diagnostic_preflight_passed=True,
  physical_run_complete=False,actual_steps=result['actual_steps'],run=run))
 print(case['name'],'SHORT_GATE_PASSED',flush=True)
a=Path(P['preflight'][0]['directory']);b=Path(P['preflight'][1]['directory']);same=[]
for p in sorted((a/'detectors').glob('*.h5')):
 with h5py.File(p) as f,h5py.File(b/'detectors'/p.name) as g:
  left=f[next(iter(f))];right=g[next(iter(g))]
  for key in ['records','fields','time_s']:
   if key in left:assert np.array_equal(left[key][:],right[key][:]),(p,key)
 same.append(p.name)
with h5py.File(a/'ledger/energy-ledger.h5') as f,h5py.File(b/'ledger/energy-ledger.h5') as g:
 x=f['energy_ledger/records'][:];y=g['energy_ledger/records'][:]
 for key in x.dtype.names:
  if not key.endswith('_rank_seconds'):assert np.array_equal(x[key],y[key]),key
save(R/'PREFLIGHT_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],binary_sha256=sha(P['binary']),binding=bound,cases=results,
 off_on_scientific_h5_bitwise=same,physical_results_produced=False,reused_observer_validation=P['reused_observer_validation']))
print('PREFLIGHT_COMPLETE',flush=True)
