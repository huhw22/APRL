"""Read-only timing/rank support of energy-run time-grid-scale field activity."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import h5py
from scipy.constants import c
LD=np.longdouble
ROOT=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
OUT=Path(__file__).resolve().parent
PROD=ROOT/'compact-production-r1/energy'
assert json.loads((PROD/'PIC_COMPLETE.json').read_text())['run']['exit_code']==0
with np.load(ROOT/'momentum_field_origin_20261105/independent_evidence/energy-strict-traces.npz') as f:
 a={k:f[k] for k in f.files}
with np.load(OUT/'rank_time_metrics.npz') as f:ranks=f['rank_records']
# rank columns are actual CSV step,native_begin,native_end,timeE,timeB.
assert np.array_equal(ranks[0,:,0],np.arange(1,23985))
assert np.max(np.abs(ranks[:,:,2].sum(axis=0)/a['native_field_J'][1:]-1))<1e-12
assert np.max(np.abs(ranks[:,:,3].sum(axis=0)/np.diff(a['time_electric_J'])-1))<1e-9
U=a['native_field_J'];te=np.r_[LD(0),np.diff(a['time_electric_J'])];tb=np.r_[LD(0),np.diff(a['time_magnetic_J'])]
ratio=te/U
with h5py.File(PROD/'detectors/entry-transport-energy.h5') as f:
 g=f['energy_plane'];assert g['complete'][()]==1
 gamma=LD(g.attrs['boost_gamma']);beta=LD(g.attrs['boost_beta']);t0=LD(g.attrs['time_lab_origin_s'])
 samples=g['records'][:]
 dt=LD(samples['time_box_s'][-1])/23984
 # Check time transform against all online samples, rather than assume origins.
 assert np.max(np.abs(gamma*(samples['time_lab_s'].astype(LD)-t0)-samples['time_box_s']))<LD('3e-25')
planes={};event_steps={};event_weights={};ids={}
for name in ['entry','und-entry','p1']:
 with h5py.File(PROD/'detectors'/(name+'-particles.h5')) as f:
  g=f['particle_plane'];assert g['complete'][()]==1
  x=g['records'].fields(['source_id','time_s','weight'])[:]
  assert len(x)==int(g['committed_records'][()])==1419264
  z=LD(g.attrs['plane_z_m'])
 n=gamma*(x['time_s'].astype(LD)-beta*z/LD(c)-t0)/dt
 w=x['weight'].astype(LD);order=np.argsort(n);cdf=np.cumsum(w[order])/w.sum()
 event_steps[name]=n;event_weights[name]=w;ids[name]=x['source_id']
 quantiles={str(q):float(n[order[min(int(np.searchsorted(cdf,q)),len(n)-1)]]) for q in [.001,.01,.1,.5,.9,.99,.999]}
 quantiles['min']=float(n.min());quantiles['max']=float(n.max())
 planes[name]=dict(z_lab_m=float(z),particle_count=len(x),weighted_crossing_step_quantiles=quantiles)
thresholds=[]
for threshold in [.01,.1,.25,.5,1.,1.5,1.8,1.9,1.99]:
 n=int(np.flatnonzero(ratio>threshold)[0]);fractions={name:float(event_weights[name][event_steps[name]<=n].sum()/event_weights[name].sum()) for name in planes}
 thresholds.append(dict(TE_U_threshold=threshold,first_step=n,time_box_ps=float(n*dt*LD('1e12')),native_field_J=float(U[n]),TE_U=float(ratio[n]),crossed_charge_fraction=fractions))
rows=[]
for n in [1,10,50,100,200,500,1000,1500,2000,4000,5324,6929,8215,8500,11000,14000,23984]:
 u=ranks[:,n-1,2];et=ranks[:,n-1,3]
 rows.append(dict(step=n,time_box_ps=float(n*dt*LD('1e12')),native_field_J=float(U[n]),TE_U=float(ratio[n]),TB_TE=float(tb[n]/te[n]) if te[n] else None,
   rank_energy_fractions=[float(v) for v in u/u.sum()],rank_TE_U=[float(v) if den>0 else None for v,den in zip(et/np.maximum(u,LD('1e-200')),u)]))
flight={}
for name in ['und-entry','p1']:
 i=np.argsort(ids['entry']);j=np.argsort(ids[name]);assert np.array_equal(ids['entry'][i],ids[name][j])
 dn=event_steps[name][j]-event_steps['entry'][i]
 flight[name]=dict(min_steps=float(dn.min()),max_steps=float(dn.max()),mean_steps=float(np.average(dn,weights=event_weights['entry'][i])),p10_steps=float(np.quantile(dn,.1)),p50_steps=float(np.quantile(dn,.5)),p90_steps=float(np.quantile(dn,.9)))
# Stored rank0 inlet positions anchor plane locations within the scanning mesh.
window=np.genfromtxt(PROD/'reference/window-rank-0.csv',delimiter=',',names=True)
dz=LD('8.293564128124092e-6');origin0=LD(window['origin_box_z_m'][0])+LD(window['step'][0])*dz
# Header frame spatial offset can be derived from a fixed z=0 plane and t0.
# boost.cpp makes box-plane z = zlab/gamma - beta*c*tbox + spatial translation.
# The actual production source audit independently reconstructed this offset.
source_events=json.loads((ROOT/'momentum_field_origin_20261105/independent_source/ACTUAL_ADMISSION_EVENTS.json').read_text())
shift=LD(source_events['box_translation_m'])
for name,plane in planes.items():
 xi=(LD(plane['z_lab_m'])/gamma+shift-origin0)/dz
 plane['window_cell_at_step0']=float(xi)
 plane['window_cell_at_step23984']=float(xi+(1-beta)*23984)
plane_fields=[]
for name in ['entry','und-entry']:
 with h5py.File(PROD/'detectors'/(name+'-transport-energy.h5')) as f:x=f['energy_plane/records'][:]
 for step in [100,200,500,1000,2000,4000,7000,8500,14000,23984]:
  k=int(np.argmin(abs(x['time_box_s']-float(step*dt))));y=x[k]
  plane_fields.append(dict(plane=name,step=step,field_power_lab_W=float(y['physical_field_power_lab_W']),cumulative_field_J=float(y['physical_field_energy_lab_J']),particle_downstream_count=int(y['physical_particle_downstream_count'])))
result=dict(case='energy',new_simulation=False,dt_box_s=float(dt),boost_gamma=float(gamma),boost_beta=float(beta),
  metrics_scope='TE/U is an independently recorded time-grid-scale indicator, not a measured kz=pi mode fraction.',
  first_thresholds=thresholds,plane_event_support=planes,per_particle_flight=flight,selected_step_metrics=rows,online_plane_records=plane_fields,
  causal_limit='The high-time-grid activity precedes most charge entering the undulator and overlaps its first-segment passage, so timing does not exclude a contribution. Existing three particle planes do not resolve BF onset within z=.01..2.166m, and global ledger/rank bins do not identify force spectra or causal contribution.')
(OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'thresholds.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['ratio_threshold','first_step','time_box_ps','native_J','entry_crossed_charge_fraction','und_entry_crossed_charge_fraction','p1_crossed_charge_fraction'])
 for x in thresholds:w.writerow([x['TE_U_threshold'],x['first_step'],x['time_box_ps'],x['native_field_J'],*[x['crossed_charge_fraction'][k] for k in ['entry','und-entry','p1']]])
print(json.dumps(dict(thresholds=thresholds,planes=planes,flight=flight),indent=2))
