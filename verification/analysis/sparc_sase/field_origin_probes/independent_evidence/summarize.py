from pathlib import Path
import json,sys,os
os.environ['MPLCONFIGDIR']='/tmp/aprl-field-evidence-mpl'
import numpy as np
sys.path.insert(0,'source/aprl/verification/performance/cases/discrete_energy_ledger')
from report import load,calculate
LD=np.longdouble
root=Path('results/sparc_sase_comparison_20261103')
out=root/'momentum_field_origin_20261105/independent_evidence'
result={}
for case in ['energy','momentum']:
 d=root/'compact-production-r1'/case/'ledger'
 a,s,ranks,steps=load(d,d/'energy-ledger.h5')
 ix,q,summary=calculate(a,s,steps)
 cumulative={n:np.r_[LD(0),np.cumsum(x,dtype=LD)] for n,x in s.items() if not n.startswith('max_') and '_prediction_scale_' not in n and n not in ['native_begin_J','native_end_J']}
 native=np.r_[s['native_begin_J'][0],s['native_end_J']]
 sample=[]
 for step in [1,10,50,100,200,500,1000,2000,4000,8000,12000,16000,20000,steps]:
  sample.append(dict(step=step,native_field_J=float(native[step]),terms={n:float(x[step]) for n,x in cumulative.items()}))
 result[case]=dict(frame='boost',ranks=ranks,steps=steps,summary=summary,final_native_field_J=float(native[-1]),initial_native_field_J=float(native[0]),field_terms_final={n:float(x[-1]) for n,x in cumulative.items()},samples=sample)
 np.savez_compressed(out/(case+'-strict-traces.npz'),native_field_J=np.asarray(native,float),**{n:np.asarray(x,float) for n,x in cumulative.items()})
 print(case,json.dumps(result[case]['field_terms_final']),flush=True)
(out/'STRICT_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
