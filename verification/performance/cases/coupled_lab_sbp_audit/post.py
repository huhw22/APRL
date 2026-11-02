import json
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import require,save
import paired_post as old
R=Path(__file__).resolve().parents[1];old.R=R;old.replay.R=R
SOURCES=['magnetic_supply','magnetic_cpml','magnetic_sync','electric_hold','electric_cpml','electric_inlet','electric_sync','reference_current']
def read(p):return np.genfromtxt(p,delimiter=',',names=True,dtype=None,encoding='utf8')
def analyze(name):
    before=old.analyze(name);d=R/name;Q=json.loads((d/'RESULT.json').read_text())['charge_abs_C'];n=before['clean_interval']['last_step']
    stages={}
    paths=sorted((d/'lab-audit').glob('whole-cowan-rank-*.csv'));require(len(paths)==16,'Missing ranks')
    for path in paths:
        require(path.with_suffix('.complete').read_text().split()[0]=='1','Incomplete whole ledger');a=read(path)
        for stage in ['begin','fixed_end','remap']:
            z=a[a['stage']==stage];require(np.array_equal(z['step'],np.arange(1,len(z)+1)),'Step coverage')
            if stage not in stages:stages[stage]={k:z[k].astype(np.longdouble) for k in z.dtype.names if k not in ['stage','step']}
            else:
                for k in stages[stage]:stages[stage][k]+=z[k]
    start,end,remap=[stages[k] for k in ['begin','fixed_end','remap']];n=len(start['H_ref_J'])
    with h5py.File(d/'ledger/energy-ledger.h5') as f:a=f['energy_ledger/records'][:]
    require(len(a)==n+1,'Ledger row count');G=np.longdouble(a['boost_gamma'][0]);b=np.longdouble(a['boost_beta'][0]);force=read(d/'particle-work-step.csv');v={k:np.zeros(n,dtype=np.longdouble) for k in force.dtype.names if k!='step'}
    for k in v:v[k][force['step']-1]=force[k]
    sourceU=sum(end[k+'_U_J'] for k in SOURCES)+remap['remap_H_J'];sourceP=sum(end[k+'_P_J'] for k in SOURCES)+remap['remap_P_J'];deltaH=remap['H_ref_J']-start['H_ref_J'];deltaP=remap['cP_ref_J']-start['cP_ref_J']
    fieldE=deltaH+end['Q_H_J']+end['WJ_particle_J']-sourceU-end['C_H_J'];fieldP=deltaP+end['Q_P_J']-end['PJ_particle_J']-sourceP-end['C_P_J']+end['GaussE_J']+end['GaussB_J']
    error=lambda x:float(np.max(abs(np.cumsum(x)))*G/Q)
    continuityH=start['H_ref_J'][1:]-remap['H_ref_J'][:-1];continuityP=start['cP_ref_J'][1:]-remap['cP_ref_J'][:-1]
    checks=dict(field_energy_cumulative_error_eV=error(fieldE),field_momentum_cumulative_error_eV=error(fieldP),interstep_H_jump_eV=float(np.max(abs(continuityH))*G/Q),interstep_P_jump_eV=float(np.max(abs(continuityP))*G/Q),direct_end_H_residual_agrees_eV=float(np.max(abs(fieldE-end['field_R_H_J']))*G/Q),direct_end_P_residual_agrees_eV=float(np.max(abs(fieldP-end['field_R_P_J']))*G/Q))
    save(d/'FIELD_CHECKS.json',checks);require(max(checks.values())<1e-6,'Independent whole-field ledger fails')
    constraints=old.merge(d,'lab-constraints',False)['records'];clean=(constraints['frozen_count']==0)&(constraints['carrier_count']==0)&(constraints['handoff_workprime_reclassified_J']==0)&(abs(a['released_charge_C'][1:])>Q*(1-1e-10))&(a['removed_macroparticles'][1:]==0);require(np.sum(clean)>20,'No clean interval')
    iext=v['external_Ez_cImpulse_J']+v['external_B_cImpulse_J'];wext=v['external_lab_work_J']/G-b*iext
    energy_gap=G*(v['delta_Kprime_J']-wext-end['WJ_particle_J']);electric_gap=G*b*(v['grid_Ez_cImpulse_J']-end['GaussE_J']);magnetic_gap=G*b*(v['grid_B_cImpulse_J']+end['PJ_particle_J'])
    CE=G*end['C_H_J'];CP=G*b*end['C_P_J'];GB=-G*b*end['GaussB_J'];total=G*(deltaH+end['Q_H_J']-sourceU+v['delta_Kprime_J']-wext+b*(deltaP+end['Q_P_J']-sourceP+v['delta_cPzprime_J']-iext));predicted=CE+CP+GB+energy_gap+electric_gap+magnetic_gap
    require(float(np.max(abs(np.cumsum(total-predicted)))/Q)<1e-6,'Particle-field identity fails');values=dict(particle_change=v['delta_Klab_J'],field_reference_change=G*(deltaH+b*deltaP),native_propagation_out=G*(end['Q_H_J']+b*end['Q_P_J']),remap_exchange=G*(remap['remap_H_J']+b*remap['remap_P_J']),Cowan_energy_operator=CE,Cowan_momentum_operator=CP,Gauss_B_momentum=GB,energy_coupling_gap=energy_gap,electric_impulse_gap=electric_gap,magnetic_impulse_gap=magnetic_gap,total_numerical_discrepancy=total,unexplained=total-predicted)
    for key in SOURCES:values[key+'_source']=G*(end[key+'_U_J']+b*end[key+'_P_J'])
    mean=lambda x:float(np.sum(x[clean])/Q);answer=dict(complete=True,original_science_bitwise_identity=True,field_checks=checks,clean_interval=dict(first_step=int(np.flatnonzero(clean)[0]+1),last_step=int(np.flatnonzero(clean)[-1]+1),steps=int(np.sum(clean)),eV_per_electron={k:mean(x) for k,x in values.items()}),physical_continuum_accuracy_accepted=False,whole_record_particle_conservation_accepted=False,full_19174_root_cause_resolved=False)
    save(d/'WHOLE_RESULT.json',answer)
    np.savetxt(d/'whole_steps_J.csv',np.column_stack([np.arange(1,n+1),clean.astype(int)]+list(values.values())),delimiter=',',header='step,clean,'+','.join(values),comments='')
    return answer
def main():
    require((R/'PIC_COMPLETE.json').exists(),'Calculation incomplete');answers={name:analyze(name) for name in ['undulator-on','drift-on']};save(R/'RESULT.json',dict(cases=answers,production_physics_changed=False,physical_lab_conservation_accepted=False,full_19174_root_cause_resolved=False));print(json.dumps(answers,indent=2),flush=True)
if __name__=='__main__':main()
