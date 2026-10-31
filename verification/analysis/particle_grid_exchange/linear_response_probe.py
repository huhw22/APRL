#!/usr/bin/env python3
"""Classical coupling-only experiment: explicit, one correction, linear implicit.

Frozen periodic CIC weights isolate the velocity/electric-current exchange.
There is no orbit advance, charge-continuity test, Maxwell curl, CPML, boost or
production speed claim. This is a small algebra/phase experiment, not ECSIM's
full algorithm and not a replacement for the APRL small physical cases.
"""
import csv,json,os,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve,expm

assert os.environ.get('SLURM_JOB_ID'),'Run the exploration in a compute allocation'
out=Path(sys.argv[1]);out.mkdir(exist_ok=False)
ng,np_=16,256;dx=1/ng;x=(np.arange(np_)+.37)/np_;z=x/dx;i=np.floor(z).astype(int);s=z-i
W=np.zeros((np_,ng));W[np.arange(np_),i]=1-s;W[np.arange(np_),(i+1)%ng]=s
weight=np.full(np_,1/np_);D=(W.T*weight)/dx;M=D@W
initial_v=.02*np.cos(2*np.pi*x)+.003*np.sin(6*np.pi*x);initial_e=.01*np.sin(2*np.pi*np.arange(ng)/ng)
A=np.block([[np.zeros((np_,np_)),-W],[D,np.zeros((ng,ng))]])
initial=np.r_[initial_v,initial_e]
def energy(v,e):return float((np.dot(weight*v,v)+dx*np.dot(e,e))/2)
H0=energy(initial_v,initial_e);results=[]
for dt in [.1,.5,1.]:
 steps=int(round(30/dt));exact=expm(A*(steps*dt))@initial
 for mode in ['explicit_old_field','one_predictor_correction','linear_implicit_response']:
  v=initial_v.copy();e=initial_e.copy();history=[]
  factor=cho_factor(np.eye(ng)+dt*dt*M/4)
  start=time.monotonic()
  for n in range(steps):
   v0=v.copy();e0=e.copy()
   if mode=='explicit_old_field':v=v0-dt*(W@e0);e=e0+dt*(D@v)
   elif mode=='one_predictor_correction':
    vp=v0-dt*(W@e0);ep=e0+dt*(D@vp)
    v=v0-dt*(W@((e0+ep)/2));e=e0+dt*(D@((v0+v)/2))
   else:
    emid=cho_solve(factor,e0+dt*(D@v0)/2)
    v=v0-dt*(W@emid);e=2*emid-e0
   history.append((n+1,energy(v,e)/H0-1))
  elapsed=time.monotonic()-start
  scale=np.sqrt(np.dot(weight*exact[:np_],exact[:np_])+dx*np.dot(exact[np_:],exact[np_:]))
  error=np.sqrt(np.dot(weight*(v-exact[:np_]),v-exact[:np_])+dx*np.dot(e-exact[np_:],e-exact[np_:]))/scale
  row={'mode':mode,'dt_normalized':dt,'steps':steps,'final_relative_energy_error':history[-1][1],
       'maximum_absolute_relative_energy_error':max(abs(h[1]) for h in history),'final_state_relative_error_vs_exact_linear_ODE':float(error),
       'nonlinear_iterations':0,'linear_field_solves_per_step':1 if mode=='linear_implicit_response' else 0,
       'toy_seconds_not_a_production_speed_estimate':elapsed}
  results.append(row)
  with (out/f'{mode}_dt{dt:g}.csv').open('w') as f:
   writer=csv.writer(f);writer.writerow(['step','relative_energy_error']);writer.writerows(history)
  print(json.dumps(row),flush=True)
implicit=[v for v in results if v['mode']=='linear_implicit_response']
assert all(v['maximum_absolute_relative_energy_error']<1e-11 for v in implicit)
assert all(implicit[k]['final_state_relative_error_vs_exact_linear_ODE']<implicit[k+1]['final_state_relative_error_vs_exact_linear_ODE'] for k in range(2))
(out/'RESULT.json').write_text(json.dumps({'algebra_checks_pass':True,'model':'Classical fixed-weight electric/particle coupling subsystem, not full PIC',
    'particles':np_,'nodes':ng,'energy_definition':'0.5 sum(w v^2) + 0.5 dx sum(E^2), q/m=-1, epsilon0=1',
    'linear_system':'[I + dt^2/(4 dx) W^T diag(w) W] E_mid = E_old + dt/(2 dx) W^T diag(w) v_old',
    'matrix_fixed_and_factorization_reused_only_in_this_frozen-weight_toy':True,
    'charge_conservation_full_Maxwell_and_production_cost_accepted':False,'results':results},indent=2)+'\n')
