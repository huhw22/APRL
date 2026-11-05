"""Read-only phase-harmonic counterfactuals; not a force or mode attribution."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,hashlib
import numpy as np
import h5py
from scipy.constants import c
TOP=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
OUT=Path(__file__).resolve().parent
BASE=TOP/'energy_startup_origin_20261105'
cache=BASE/'independent/independent_phase_arrays.npz'
with np.load(cache) as f: dp=f['phase_residual'];ids=f['source_id']
with h5py.File(TOP/'compact-preflight-r1/particles.h5','r') as f:
 a=f['particles/records'].fields(['source_id','macro_weight'])[:]
assert np.array_equal(a['source_id'],ids)
w=a['macro_weight'];del a
W=w.sum();group=np.arange(len(w))//10240;wg=np.bincount(group,weights=w)
lam=json.loads((TOP/'compact-production-r1/plan.json').read_text())['physics']['wavelength_m']
def timeplane(name):
 p=TOP/'compact-production-r1/energy/detectors'/(name+'-particles.h5')
 with h5py.File(p,'r') as f:
  g=f['particle_plane'];assert g['complete'][()]==1
  a=g['records'].fields(['source_id','time_s'])[:]
  assert len(a)==g['committed_records'][()]
 a=a[np.argsort(a['source_id'])];assert np.array_equal(a['source_id'],ids)
 return a['time_s']
tund=timeplane('und-entry');tp=timeplane('p1')
def phase(t):
 t=t.astype(np.longdouble)
 return np.exp(1j*np.asarray(np.remainder((t-np.longdouble(tund[0]))*np.longdouble(c/lam),1)*np.longdouble(2*np.pi),float))
pin=phase(tund);pactual=phase(tp);pref=pactual*np.exp(-1j*dp)
within=dp-dp.reshape(-1,16).mean(axis=1).repeat(16)
def summary(ph):
 b=np.sum(w*ph)/W
 z=(np.bincount(group,weights=w*ph.real)+1j*np.bincount(group,weights=w*ph.imag))/wg
 rms=np.sqrt(np.sum(wg*abs(z)**2)/W)
 return dict(global_BF=float(abs(b)),local_RMS=float(rms))
rows={'reference':summary(pref),'actual':summary(pactual),'within_actual':summary(pref*np.exp(1j*within))}
coeff={}
for name,ns in [('n10',[10]),('n10_20_30',[10,20,30]),('n1',[1]),('n1_to_3',[1,2,3]),('n1_to_8',list(range(1,9)))]:
 columns=[]
 for n in ns:
  z=pin**n
  for v in [z.real,z.imag]:columns.append(v-v.reshape(-1,16).mean(axis=1).repeat(16))
 X=np.column_stack(columns);beta=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w*within));fit=X@beta
 row=summary(pref*np.exp(1j*fit));row.update(fitted_phase_rms_rad=float(np.sqrt(np.sum(w*fit**2)/W)),within_phase_variance_fraction=float(1-np.sum(w*(within-fit)**2)/np.sum(w*within**2)))
 rows[name]=row;coeff[name]=beta.tolist()
 del X,columns,fit
assert abs(rows['reference']['global_BF']-2.421389838119538e-5)<1e-9
assert abs(rows['actual']['global_BF']-.004450123462378471)<1e-9
result=dict(case='energy',new_pic=False,read_only=True,particle_count=len(w),definition='Weighted within-family phase regression on harmonics of recorded undulator-entrance optical phase, then exact finite phase-map summation on the same full particle cohort.',counterfactuals=rows,coefficients=coeff,limitation='Descriptive phase rearrangement only. Not source-force attribution, not a direct projection of Cowan modes. BF and phase residual retain raw whole-Maxwell-step crossing-time reconstruction. Input Fawley shot noise breaks ideal exact pi-pair symmetry.',symmetry='For perfectly uniform 16-clone pi pairs, any perturbation delta_phi(phi+pi)=delta_phi(phi), including only n10 and multiples, leaves fundamental BF zero exactly; this exact symmetry is not assumed for the actual noisy input.')
(OUT/'PHASE_PROJECTION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
