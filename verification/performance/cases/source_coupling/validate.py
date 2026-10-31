#!/usr/bin/env python3
"""Independent checks of closed prototype evidence and rendered-data units."""
import argparse, hashlib, json, re, subprocess
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);p.add_argument('analysis',type=Path);a=p.parse_args()
campaign=a.campaign.resolve();analysis=a.analysis.resolve()
plan=json.loads((campaign/'plan.json').read_text());ws=Path(plan['workspace']);source=ws/'source/aprl'
complete=json.loads((campaign/'COMPLETE.json').read_text());assert complete['checks_passed']
frozen=json.loads((campaign/'FROZEN.json').read_text())
for path,digest in frozen.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
for name in ['source_coupling.h','prototype.cpp','checks.cpp','external.cpp']:
 assert (source/'verification/performance/cases/source_coupling'/name).read_bytes()==(campaign/'sources'/name).read_bytes(),name
subprocess.run(['git','diff','--exit-code','2118ce8','--','src','CMakeLists.txt'],cwd=source,check=True,stdout=subprocess.PIPE)
accepted=ws/'build/charge-matched-inlet-20261030/aprl'
assert hashlib.sha256(accepted.read_bytes()).hexdigest()=='79d643bc8ade4db4d9e4488a367cb1bf5ef74e7c602dc3a2094462709aed9435'
for name,digest in complete['binaries_sha256'].items():assert hashlib.sha256((Path(plan['build'])/name).read_bytes()).hexdigest()==digest
normalization=json.loads((analysis/'NORMALIZATION.json').read_text());ep=np.finfo(np.longdouble).eps
roundtrip=0.;identity=0.;continuity=0.;prediction=0.;rho_cases=0;closed_cases=0
for case in plan['cases']:
 name=case['name'];directory=campaign/name;result=complete['cases'][name];assert result['complete'];closed_cases+=1
 assert result['magnetic_CPML_calls']==result['electric_CPML_calls']==case['steps']
 if not case['audit']:continue
 assert (directory/'discrete-field-rank-00000.complete').read_text().split()==['1',str(case['steps'])]
 x=np.genfromtxt(directory/'history.csv',delimiter=',',names=True,dtype=np.longdouble)
 ledger=np.genfromtxt(directory/'discrete-field-rank-00000.csv',delimiter=',',names=True,dtype=np.longdouble)
 assert np.array_equal(x['step'],np.arange(1,case['steps']+1))
 assert all(np.isfinite(x[k]).all() for k in x.dtype.names)
 actual=x['delta_particle_box_J']-x['prescribed_work_J']-x['mesh_mid_J']
 reported=x['coupling_gap_J']+x['pusher_residual_J']
 scale=np.max(abs(x['delta_particle_box_J'])+abs(x['prescribed_work_J'])+abs(x['mesh_mid_J']))
 error=float(np.max(abs(actual-reported)));identity=max(identity,error)
 assert error<=512*ep*max(scale,np.longdouble('1e-300')),(name,error,scale)
 assert np.max(abs(np.cumsum(ledger['field_operator_residual_J'])-x['field_operator_residual_J']))<=512*ep*max(scale,np.max(abs(ledger['native_begin_J'])),np.longdouble('1e-300'))
 rho_cases+=1;continuity=max(continuity,result['relative_continuity_residual_max']);assert continuity<5e-12
 for field,unit in [('electric','V_per_m'),('magnetic','T')]:
  ratio=abs(ledger['max_'+field+'_prediction_error_'+unit])/np.maximum(abs(ledger[field+'_prediction_scale_'+unit]),np.longdouble('1e-300'))
  prediction=max(prediction,float(np.max(ratio)));assert np.all(ratio<=1024*np.finfo(float).eps)
 if name in [k+'-'+m for k in ['sphere','drift','undulator'] for m in ['baseline','coupled']]:
  relative=np.genfromtxt(analysis/(name+'-relative.csv'),delimiter=',',names=True,dtype=np.longdouble)
  D=np.longdouble(normalization[case['kind']]['primary_fixed_D_J'])
  columns=[('actual_exchange_gap_percent',reported),('grid_exchange_gap_percent',x['coupling_gap_J']),('pusher_residual_percent',x['pusher_residual_J']),('field_operator_residual_percent',x['field_operator_residual_J']),('strict_total_residual_percent',reported+x['field_operator_residual_J'])]
  for key,original in columns:
   restored=relative[key]*D/100
   e=float(np.max(abs(restored-original))/max(np.max(abs(original)),np.longdouble('1e-300')))
   roundtrip=max(roundtrip,e);assert e<5e-15,(name,key,e)
vac={mode:complete['cases']['vacuum-'+mode] for mode in ['baseline','coupled']}
assert vac['baseline']['field_hash']==vac['coupled']['field_hash']
v=np.genfromtxt(campaign/'vacuum-baseline/discrete-field-rank-00000.csv',delimiter=',',names=True)
assert v['native_begin_J'][0]>0 and v['native_end_J'][-1]>0
# Reconstruct the cold-sphere kinetic gain directly from final momenta, without grid work.
kinetic={};mass=np.longdouble(float('9.1093837139e-31'))*np.longdouble(float(1e8/4096));c=np.longdouble(299792458)
for name in ['sphere-baseline','sphere-coupled','sphere-halfdt-baseline','sphere-halfdt-coupled']:
 states=np.genfromtxt(campaign/name/'particles-final.csv',delimiter=',',names=True,dtype=np.longdouble)
 u2=sum(states[key]**2 for key in ['ux','uy','uz']);K=np.sum(mass*c*c*u2/(np.sqrt(1+u2)+1))
 x=np.genfromtxt(campaign/name/'history.csv',delimiter=',',names=True,dtype=np.longdouble)
 err=float(abs(K-x['delta_particle_box_J'][-1])/abs(K));assert err<5e-14,(name,err)
 kinetic[name]={'independent_final_momenta_delta_K_J':float(K),'relative_difference':err}
# The two stopping rules measure different quantities; report the state comparison as a diagnostic.
agreement={}
for kind in ['sphere','drift','undulator']:
 reference=ws/'results/source_coupling_prototype_20261031/r5'/(kind+'-coupled')/'particles-final.csv'
 old=np.genfromtxt(reference,delimiter=',',names=True);new=np.genfromtxt(campaign/(kind+'-coupled')/'particles-final.csv',delimiter=',',names=True)
 assert np.array_equal(old['id'],new['id'])
 du=max(float(np.max(abs(old[k]-new[k]))) for k in ['ux','uy','uz']);dx=max(float(np.max(abs(old[k]-new[k]))) for k in ['x','y','z'])
 agreement[kind]={'max_normalized_momentum_component_difference':du,'max_position_component_difference_m':dx}
result={'passed':True,'frozen_files':len(frozen),'closed_cases':closed_cases,'strict_audit_cases':rho_cases,'maximum_relative_continuity_residual':continuity,'maximum_relative_field_operator_prediction_error':prediction,'maximum_actual_exchange_identity_roundoff_J':identity,'maximum_relative_percent_to_J_roundtrip_error':roundtrip,'cold_sphere_kinetic_energy_from_final_momenta':kinetic,'optimized_vs_earlier_state_change_iteration_final_particle_state':agreement,'comparison_uses_distinct_stopping_rules_not_a_physics_acceptance_gate':True,'vacuum_nonzero_and_full_field_hash_equal':True,'production_core_equal_to_parent_2118ce8':True,'accepted_production_binary_hash_unchanged':True,'independent_prototype_sources_equal_frozen_build':True,'figures_require_visual_review':True,'full_production_physics_accepted':False}
(analysis/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
