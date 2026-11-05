"""Strict energy-weighted field-step overlap from native ledger square terms."""
from pathlib import Path
import json
import numpy as np
LD=np.longdouble
ROOT=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
OUT=Path(__file__).resolve().parent
with np.load(ROOT/'momentum_field_origin_20261105/independent_evidence/energy-strict-traces.npz') as f:a={k:f[k] for k in f.files}
with np.load(OUT/'rank_time_metrics.npz') as f:r=f['rank_records']
old=r[:,:,1].sum(axis=0);end=r[:,:,2].sum(axis=0)
delta=lambda key:np.diff(a[key])
te=delta('time_electric_J');tb=delta('time_magnetic_J');D=te-tb
# TE compares actual E after electric CPML and before inlet overwrite to old E.
# TB compares B just after raw magnetic curl (before magnetic CPML) to old B.
# Therefore match exactly these end-state norms; do not subtract electric CPML.
new=end-delta('remap_J')-delta('inlet_overwrite_J')-delta('magnetic_cpml_J')
assert np.all(new>0) and np.all(old>0)
S=old+new
cosine=(S-D)/(2*np.sqrt(old*new))
nearby_new=new-delta('electric_cpml_J')
approx=(old+nearby_new-D)/(2*np.sqrt(old*nearby_new))
assert np.all(cosine>=-1-LD('1e-10')) and np.all(cosine<=1+LD('1e-10'))
rows=[]
for n in [1,10,50,100,200,237,391,500,1000,1314,1500,1537,2000,4000,5324,6936,8215,8469,11000,14000,23984]:
 i=n-1
 rows.append(dict(step=n,native_old_J=float(old[i]),matched_native_new_J=float(new[i]),delta_square_D_J=float(D[i]),cosine_energy_weighted_old_vs_new=float(cosine[i]),opposite_fraction_D_over_2S=float(D[i]/(2*S[i])),electric_cpml_alternative_difference=float(approx[i]-cosine[i])))
threshold=[]
for t in [.9,.5,0.,-.25,-.5,-.8,-.9,-.95,-.99]:
 matches=np.flatnonzero(cosine<t)
 i=int(matches[0]) if len(matches) else None
 threshold.append(dict(cosine_below=t,first_step=i+1 if i is not None else None,value=float(cosine[i]) if i is not None else None))
# Projected energies of the step-average and step-difference fields are exact.
# mean-field energy = (2S-D)/4, difference-half field energy = D/4.
result=dict(frame='boost',new_simulation=False,actual_recorded_square_terms=True,definition='Old/new E/B vectors in the positive native-energy inner product; magnetic endpoint before magnetic CPML, electric endpoint after electric CPML before inlet overwrite.',
 formula='cos=(Uold+Unew-(TE-TB))/(2 sqrt(Uold Unew))',electric_cpml_not_subtracted_reason='TE already uses after-electric-CPML values.',max_electric_cpml_alternative_cosine_difference=float(np.max(np.abs(approx-cosine))),
 min_cosine=float(cosine.min()),max_cosine=float(cosine.max()),thresholds=threshold,selected_steps=rows,
 causal_limit='A near-minus-one overlap establishes actual E/B change dominated by sign reversal at the fixed boost grid during the Maxwell field step. It does not independently measure kz, prove exact Jordan dominance, or identify which part of the fields does work on particles.')
(OUT/'FLIP_OVERLAP.json').write_text(json.dumps(result,indent=2)+'\n')
np.savez_compressed(OUT/'flip_overlap_traces.npz',step=np.arange(1,len(old)+1),old_J=old,new_matched_J=new,D_J=D,cosine=cosine,D_over_2S=D/(2*S))
print(json.dumps(result,indent=2))
