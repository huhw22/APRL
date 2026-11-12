"""Four matched sample/control cases per method, guarded by the short gate."""
import json,os,sys
from pathlib import Path
from common import frozen,binding,launch,validate,save,sha

R=Path(sys.argv[1]).resolve();method=sys.argv[2];P=frozen(R)
gate=json.loads((R/'PREFLIGHT_COMPLETE.json').read_text());assert gate['complete'] and gate['binary_sha256']==sha(P['binary'])
assert method in ['linear-binomial3','quintic-compensated']
bound=binding(R,P,'binding-'+method+'.log');completed=[]
for case in P['cases']:
 if case['shape']!=method:continue
 d=Path(case['directory']);assert sha(case['input'])==case['input_sha256']
 print('START',case['name'],flush=True)
 run=launch(P,case);assert run['exit_code']==0
 result=validate(case);result.update(name=case['name'],shape=method,seed=case['seed'],kind=case['kind'],run=run,analysis_pending=True)
 save(d/'COMPLETE.json',result);completed.append(result)
 print('COMPLETE',case['name'],flush=True)
save(R/('COMPLETE_'+method+'.json'),dict(complete=True,cases=completed,binding=bound,job_id=os.environ['SLURM_JOB_ID'],postprocessing_pending=True))
