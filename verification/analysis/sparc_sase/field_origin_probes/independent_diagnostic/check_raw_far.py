from pathlib import Path
import h5py,numpy as np,json,hashlib
root=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
out=root/'momentum_field_origin_20261105/independent_diagnostic'
c=np.longdouble(299792458.)
# Values defined in frozen src/physical_constants.h / constants.cpp.
mu0=np.longdouble(1.25663706212e-6)
eps0=np.longdouble(8.8541878128e-12)
rows=[]
for case in ['energy','momentum']:
 d=root/'compact-production-r1'/case/'detectors'
 with h5py.File(d/'exit-transport-energy.h5') as f:
  g=f['energy_plane'];assert int(g['complete'][()])==1;a=g['records'][:];G=np.longdouble(g.attrs['boost_gamma']);beta=np.longdouble(g.attrs['boost_beta']);area=np.longdouble(g.attrs['dx_m'])*np.longdouble(g.attrs['dy_m'])
 with h5py.File(d/'far.h5') as f:
  g=f['field_plane'];assert int(g['complete'][()])==1;t=g['time_s'][:];v=g['fields'];assert len(t)==len(v)==int(g['committed_samples'][()])
  choose=np.unique(np.rint(np.linspace(0,len(t)-1,65)).astype(int));rec=[]
  for n in choose:
   j=int(np.searchsorted(a['time_lab_s'],t[n]));assert a['time_lab_s'][j]==t[n]
   b=v[n];E=b['electric_V_per_m'].astype(np.longdouble);B=b['magnetic_T'].astype(np.longdouble)
   power=area/mu0*(E[...,0]*B[...,1]-E[...,1]*B[...,0]);density=area*.5*(eps0*np.sum(E*E,axis=-1)+np.sum(B*B,axis=-1)/mu0)
   Eb=E.copy();Bb=B.copy();Eb[...,0]=G*(E[...,0]-beta*c*B[...,1]);Eb[...,1]=G*(E[...,1]+beta*c*B[...,0]);Bb[...,0]=G*(B[...,0]+beta*E[...,1]/c);Bb[...,1]=G*(B[...,1]-beta*E[...,0]/c)
   Plab=float(np.sum(power[10:-10,10:-10]));Ulab=float(np.sum(density[10:-10,10:-10]));Ub=float(np.sum(.5*area*(eps0*np.sum(Eb[10:-10,10:-10]**2,axis=-1)+np.sum(Bb[10:-10,10:-10]**2,axis=-1)/mu0)))
   # Forward/backward transverse characteristic energies, sum is transverse lab density.
   Fp=area*eps0/4*np.sum((E[...,0]+c*B[...,1])**2+(E[...,1]-c*B[...,0])**2)
   Fm=area*eps0/4*np.sum((E[...,0]-c*B[...,1])**2+(E[...,1]+c*B[...,0])**2)
   rec.append(dict(index=int(n),time_s=float(t[n]),field_power_physical_W=Plab,ledger_power_physical_W=float(a['physical_field_power_lab_W'][j]),field_density_physical_J_per_m=Ulab,ledger_density_physical_J_per_m=float(a['physical_field_energy_per_length_lab_J_per_m'][j]),inverse_boost_density_physical_J_per_m=Ub,max_E_lab=float(np.max(np.abs(E))),max_cB_lab=float(np.max(np.abs(c*B))),max_E_boost=float(np.max(np.abs(Eb))),max_cB_boost=float(np.max(np.abs(c*Bb))),forward_characteristic_full_J_per_m=float(Fp),backward_characteristic_full_J_per_m=float(Fm)))
  scale=max(abs(q['ledger_power_physical_W']) for q in rec)
  rows.append(dict(case=case,matched_samples=len(rec),power_max_abs_error_W=max(abs(q['field_power_physical_W']-q['ledger_power_physical_W']) for q in rec),power_error_relative_to_peak=max(abs(q['field_power_physical_W']-q['ledger_power_physical_W']) for q in rec)/scale,density_error_relative_to_peak=max(abs(q['field_density_physical_J_per_m']-q['ledger_density_physical_J_per_m']) for q in rec)/max(abs(q['ledger_density_physical_J_per_m']) for q in rec),records=rec))
(out/'RAW_FAR_CHECK.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps([{k:v for k,v in q.items() if k!='records'} for q in rows],indent=2))
