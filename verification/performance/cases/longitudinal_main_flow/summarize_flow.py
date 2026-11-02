"""Describe closed short-run costs without extrapolating to production."""
import json,re,sys
from pathlib import Path
import numpy as np
r=Path(sys.argv[1]).resolve();plan=json.loads((r/'plan.json').read_text());assert json.loads((r/'COMPLETE.json').read_text())['program_preflight_pass'];cases={}
for c in plan['cases']:
 d=Path(c['directory']);points={int(s):float(t) for s,t in re.findall(r'phase=progress step=(\d+) elapsed_s=([^ ]+)',(d/'aprl.log').read_text())};assert 80 in points and 160 in points
 iterations=[]
 for p in (d/'reference').glob('source-coupling-rank-*.csv'):
  a=np.genfromtxt(p,delimiter=',',names=True);iterations.append(float(np.max(a['iterations'])))
 cases[c['name']]=dict(nodes=c['nodes'],ranks=c['ranks'],threads=8,seconds_per_step_80_to_160=(points[160]-points[80])/80,max_iterations=max(iterations),strict_ledger=True,profile_enabled=False)
comparison={}
for kind in ['drift','undulator']:
 base=cases[kind+'-baseline']['seconds_per_step_80_to_160'];linear=cases[kind+'-linear']['seconds_per_step_80_to_160'];q=cases[kind+'-quintic']['seconds_per_step_80_to_160'];comparison[kind]=dict(old=base,new_linear=linear,new_quintic=q,quintic_increment_over_new_linear_percent=100*(q/linear-1),samples=1,production_cost_accepted=False)
out=dict(cases=cases,comparisons=comparison,scope='One same-allocation short pair per geometry; no production cost extrapolation; fixed progression interval 80-160',physical_acceptance=False);(r/'COST.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
