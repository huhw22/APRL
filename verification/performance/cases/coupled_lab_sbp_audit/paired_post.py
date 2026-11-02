import json,os
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import work_replay_post as replay
R=Path(__file__).resolve().parents[1];replay.R=R
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def read(p):return np.genfromtxt(p,delimiter=',',names=True,dtype=None,encoding='utf8')
def merge(d,prefix,staged=True):
 paths=[p for p in sorted((d/'lab-audit').glob(prefix+'-rank-*.csv')) if p.stem.rsplit('-',1)[-1].isdigit()];assert len(paths)==16;out={}
 for p in paths:
  assert p.with_suffix('.complete').read_text().split()[0]=='1';a=read(p)
  groups=np.unique(a['stage']) if staged else ['records']
  for label in groups:
   rows=a[a['stage']==label] if staged else a;assert np.array_equal(rows['step'],np.arange(1,len(rows)+1))
   value={k:rows[k].astype(np.longdouble) for k in a.dtype.names if k not in ['step','stage']}
   if label not in out:out[label]=value
   else:
    for k in value:assert value[k].shape==out[label][k].shape;out[label][k]+=value[k]
 return out
def identity(a,b):
 for name in ['release','aux','exit']:
  x,_=replay.rows(a/'detectors'/(name+'.h5'),'particle_plane');y,_=replay.rows(b/'detectors'/(name+'.h5'),'particle_plane');x.sort(order='particle_id');y.sort(order='particle_id');assert np.array_equal(x,y),name
 with h5py.File(a/'detectors/exit-field.h5') as f,h5py.File(b/'detectors/exit-field.h5') as g:
  x=f['field_plane'];y=g['field_plane'];assert x['complete'][()]==y['complete'][()]==1;assert np.array_equal(x['time_s'][:],y['time_s'][:]) and x['fields'].shape==y['fields'].shape
  for n in range(0,len(x['fields']),16):assert np.array_equal(x['fields'][n:n+16],y['fields'][n:n+16])
 return True
def analyze(name):
 d=R/name;old=R.parent/'r1'/('undulator-off' if name=='undulator-on' else 'drift-on');assert identity(old,d)
 work=replay.run_case(name);Q=work['charge_abs_C'];s=merge(d,'paired-momentum');c=merge(d,'lab-constraints',False)['records'];a,_=replay.rows(d/'ledger/energy-ledger.h5','energy_ledger');force=read(d/'particle-work-step.csv');n=len(a)-1
 f={k:np.zeros(n,dtype=np.longdouble) for k in force.dtype.names if k!='step'}
 for k in f:f[k][force['step']-1]=force[k]
 G=np.longdouble(a['boost_gamma'][0]);beta=np.longdouble(a['boost_beta'][0]);assert abs(float(np.sum(f['delta_Klab_J']))-a['delta_particle_kinetic_lab_J'][-1])/Q<1e-6
 def sums(stage):return sum(s[stage]['stress_'+face+'_J'] for face in ['xm','xp','ym','yp','zm','zp'])
 QE=sum(s['magnetic_identity']['port_stress_'+face+'_J'] for face in ['xm','xp','ym','yp','zm','zp']);QB=sums('electric_identity');FE=s['magnetic_identity']['interior_gauss_electric_impulse_J'];FB=s['electric_identity']['gauss_magnetic_impulse_J'];PY=s['magnetic_identity']['yee_magnetic_delta_J'];PC=s['magnetic_identity']['cowan_magnetic_extra_J'];PA=s['electric_identity']['full_ampere_curl_delta_J'];PJ=s['electric_identity']['full_current_delta_J'];PH=s['electric_identity']['held_node_delta_J']
 maxE=max(abs(PY+QE+FE))/Q;maxB=max(abs(PA+QB+FB))/Q;assert maxE<1e-7 and maxB<1e-7,(name,maxE,maxB)
 parts={}
 for k in ['native_Uprime_J','paired_cPzprime_J','lab_charge_J']:
  v={tag:s[tag][k] for tag in s};parts[k]=dict(delta=v['remap']-v['begin_native'],cpml=v['magnetic_cpml']-v['magnetic_curl']+v['electric_cpml']-v['electric_curl_local'],inlet=v['inlet_overwrite']-v['electric_cpml'],remap=v['remap']-v['electric_mpi'])
 def op(k):return parts[k]['cpml']+parts[k]['inlet']+parts[k]['remap']
 impulseExternal=f['external_Ez_cImpulse_J']+f['external_B_cImpulse_J'];impulseConstraint=c['frozen_cImpulse_z_J']+c['carrier_cImpulse_z_J']
 RP=parts['paired_cPzprime_J']['delta']+f['delta_cPzprime_J']+QE+QB-op('paired_cPzprime_J')-PH-impulseExternal-impulseConstraint
 forceCompatibility=f['delta_cPzprime_J']-impulseExternal-impulseConstraint+PJ-FE
 pred=forceCompatibility-FB+PC;predictionError=RP-pred;assert max(abs(np.cumsum(predictionError)))/Q<1e-6,(name,max(abs(np.cumsum(predictionError)))/Q)
 frozenWork=np.diff(a['frozen_constraint_work_J']).astype(np.longdouble);carrierWork=np.diff(a['carrier_constraint_work_J']).astype(np.longdouble);boundaryWork=np.diff(a['boundary_source_work_J']).astype(np.longdouble)
 assert max(abs(np.cumsum(frozenWork)-np.cumsum(c['frozen_workprime_J'])))/Q<1e-7
 assert max(abs(np.cumsum(carrierWork)-np.cumsum(c['carrier_workprime_J']+c['handoff_workprime_reclassified_J'])))/Q<1e-7
 # Exact native operator energy ledger; Q is its independently paired curl
 # exchange, not the physically aligned plane Poynting diagnostic.
 de={}
 for p in sorted((d/'ledger').glob('discrete-field-rank-*.csv')):
  assert p.with_suffix('.complete').read_text().split()[0]=='1';v=read(p)
  for k in v.dtype.names:
   if k=='step':continue
   if k not in de:de[k]=v[k].astype(np.longdouble)
   else:de[k]+=v[k]
 extK=f['external_lab_work_J']/G-beta*impulseExternal
 RE=parts['native_Uprime_J']['delta']+de['discrete_curl_exchange_J']+f['delta_Kprime_J']-extK-frozenWork-carrierWork-boundaryWork-op('native_Uprime_J')
 numericE=de['time_electric_J']+de['time_magnetic_J']+de['current_time_pairing_J']+de['cowan_operator_J']
 energyIdentityError=RE-numericE;assert max(abs(np.cumsum(energyIdentityError)))/Q<1e-6,(name,max(abs(np.cumsum(energyIdentityError)))/Q)
 # An interval with no frozen sources, carriers or handoff is the clearest
 # direct audit of active physical coupling. Source boundary terms remain.
 clean=(c['frozen_count']==0)&(c['carrier_count']==0)&(c['handoff_workprime_reclassified_J']==0)&(abs(a['released_charge_C'][1:])>Q*(1-1e-10))&(a['removed_macroparticles'][1:]==0)
 assert clean.sum()>20
 Ecompat=f['grid_Ez_cImpulse_J']-FE;Bcompat=f['grid_B_cImpulse_J']+PJ
 lab=G*(RE+beta*RP);frozenLab=G*(c['frozen_workprime_J']+beta*c['frozen_cImpulse_z_J']);carrierLab=G*(c['carrier_workprime_J']+beta*c['carrier_cImpulse_z_J'])
 def integrated(v,mask=None):return float(np.sum(v if mask is None else v[mask])/Q)
 result=dict(original_physics_bitwise=True,fixed_plane_particle_gain_eV=work['plane_gain_eV'],whole_cohort_particle_gain_eV=integrated(f['delta_Klab_J']),frozen_constraint_lab_input_eV=integrated(frozenLab),carrier_constraint_lab_input_raw_eV=integrated(carrierLab),handoff_energy_reclassified_boost_eV=integrated(c['handoff_workprime_reclassified_J']),native_lab_discrepancy_eV=integrated(lab),Gamma_energy_numerical_eV=integrated(G*RE),Gamma_beta_momentum_discrepancy_eV=integrated(G*beta*RP),Gamma_beta_force_compatibility_eV=integrated(G*beta*forceCompatibility),Gamma_beta_Cowan_momentum_eV=integrated(G*beta*PC),Gamma_beta_divB_eV=integrated(-G*beta*FB),clean_interval=dict(first_step=int(np.flatnonzero(clean)[0]+1),last_step=int(np.flatnonzero(clean)[-1]+1),steps=int(clean.sum()),particle_gain_eV=integrated(f['delta_Klab_J'],clean),lab_discrepancy_eV=integrated(lab,clean),Gamma_energy_numerical_eV=integrated(G*RE,clean),Gamma_beta_momentum_discrepancy_eV=integrated(G*beta*RP,clean),Gamma_beta_electric_force_compatibility_eV=integrated(G*beta*Ecompat,clean),Gamma_beta_magnetic_force_compatibility_eV=integrated(G*beta*Bcompat,clean),Gamma_beta_Cowan_momentum_eV=integrated(G*beta*PC,clean),Gamma_beta_divB_eV=integrated(-G*beta*FB,clean)),checks=dict(max_electric_SBP_error_eV=float(maxE),max_magnetic_SBP_error_eV=float(maxB),max_cumulative_momentum_identity_error_eV=float(max(abs(np.cumsum(predictionError)))/Q),max_cumulative_energy_identity_error_eV=float(max(abs(np.cumsum(energyIdentityError)))/Q)),physical_lab_conservation_accepted=False,limitations=['Native energy/momentum bookkeeping includes explicit numerical time/Cowan channels; it is not a physical conservation acceptance.','Raw carrier impulse is reported separately; energy handoff reclassification has not been assigned a fabricated momentum correction.','Clean interval excludes all preparation sources and carrier/handoff/removal events.'])
 data=np.column_stack([np.arange(1,n+1),np.cumsum(f['delta_Klab_J'])/Q,np.cumsum(frozenLab)/Q,np.cumsum(carrierLab)/Q,np.cumsum(lab)/Q,np.cumsum(G*RE)/Q,np.cumsum(G*beta*RP)/Q,np.cumsum(G*beta*forceCompatibility)/Q,np.cumsum(-G*beta*FB)/Q,np.cumsum(G*beta*PC)/Q]);np.savetxt(d/'lab_balance_paired.csv',data,delimiter=',',header='step,particle_delta_lab_eV,frozen_lab_input_eV,carrier_lab_input_raw_eV,lab_discrepancy_eV,Gamma_energy_numerical_eV,Gamma_beta_momentum_discrepancy_eV,Gamma_beta_force_compatibility_eV,Gamma_beta_divB_eV,Gamma_beta_Cowan_momentum_eV',comments='');save(d/'PAIRED_RESULT.json',result);print(name,json.dumps(result),flush=True);return result
def main():
 assert (R/'PIC_COMPLETE.json').exists() and (R/'EXIT_CODE').read_text().strip()=='0'
 result={name:analyze(name) for name in ['undulator-on','drift-on']};save(R/'RESULT.json',dict(cases=result,full_19174_root_cause_resolved=False,production_physics_changed=False))
 fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
 for ax,(name,res) in zip(axes,result.items()):
  x=read(R/name/'lab_balance_paired.csv');
  for key,label in [('Gamma_energy_numerical_eV','Energy numerical channels'),('Gamma_beta_momentum_discrepancy_eV','Momentum discrepancy'),('lab_discrepancy_eV','Sum')]:ax.plot(x['step'],x[key],label=label)
  ax.set(title=name,xlabel='Field step',ylabel='Cumulative eV / electron');ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.savefig(R/'paired_lab_balance.png',dpi=180);fig.savefig(R/'paired_lab_balance.pdf');plt.close(fig)
if __name__=='__main__':main()
