#!/usr/bin/env python3
"""Independent reconstruction; numerical corrections are never residual-derived."""
import csv,json
from pathlib import Path
import h5py
import numpy as np
LD=np.longdouble
EPS=LD(np.finfo(float).eps)

def require(b,message):
    if not b:raise RuntimeError(message)
def save(path,obj):Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def validate_columns(v):
    energy=['native_begin_J','native_end_J','discrete_curl_exchange_J','applied_current_mid_J',
      'time_electric_J','time_magnetic_J','current_time_pairing_J','cowan_operator_J',
      'magnetic_cpml_J','electric_cpml_J','inlet_overwrite_J','remap_J']
    scale=sum(abs(v[n]) for n in energy);bound=512*EPS*np.maximum(scale,LD('1e-100'))
    expected=v['native_end_J']-v['native_begin_J']+v['discrete_curl_exchange_J']+v['applied_current_mid_J']
    for n in energy[4:]:expected-=v[n]
    require(np.all(abs(expected-v['field_operator_residual_J'])<=bound),'Serialized operator ledger identity differs')
    require(np.all(abs(v['current_time_pairing_J']-(v['applied_current_mid_J']-v['applied_current_old_J']))<=bound),'Current time pairing differs')
    require(np.all(abs(v['field_operator_residual_J'])<=bound),'Field-operator residual above arithmetic bound')
    require(np.all(abs(v['field_operator_residual_J']-v['electric_prediction_error_work_J']-v['magnetic_prediction_error_work_J'])<=bound),'Observed residual differs from independent stencil defects')
    for kind,unit in [('electric','V_per_m'),('magnetic','T')]:
        require(np.all(v[f'max_{kind}_prediction_error_{unit}']<=512*EPS*np.maximum(v[f'{kind}_prediction_scale_{unit}'],LD('1e-100'))),'Pointwise '+kind+' prediction mismatch')
    return bound

def inspect(root,plan):
    root=Path(root);case=Path(plan['cases']['on']);off=Path(plan['cases']['off'])
    files=sorted((case/'ledger').glob('discrete-field-rank-*.csv'))
    require([p.name for p in files]==[f'discrete-field-rank-{r:05d}.csv' for r in range(16)],'Missing rank ledger')
    require(not list((off/'ledger').glob('discrete-field-*')),'Disabled flag wrote discrete output')
    sums={};reports={};first=None
    for path in files:
        require(path.with_suffix('.complete').read_text().split()==['1','257'],'Incomplete rank ledger')
        with path.open() as f:reader=csv.DictReader(f);rows=list(reader);names=reader.fieldnames
        require([int(r['step']) for r in rows]==list(range(1,258)),'Wrong step coverage')
        values={n:np.array([LD(r[n]) for r in rows]) for n in names if n!='step'}
        require(all(np.isfinite(v).all() for v in values.values()),'Nonfinite numerical ledger')
        bound=validate_columns(values)
        if first is None:first={k:v.copy() for k,v in values.items()}
        reports[path.name]={'checks_pass':True,'max_field_operator_residual_J':float(np.max(abs(values['field_operator_residual_J']))),
                            'max_arithmetic_bound_J':float(np.max(bound))}
        for n,v in values.items():
            if n not in sums:sums[n]=np.zeros_like(v)
            if n.startswith('max_') or '_prediction_scale_' in n:sums[n]=np.maximum(sums[n],v)
            else:sums[n]+=v
    mutated={k:v.copy() for k,v in first.items()};mutated['time_electric_J'][100]+=LD('1e-16')
    caught=False
    try:validate_columns(mutated)
    except RuntimeError:caught=True
    require(caught,'Reader accepted a corrupted numerical term')
    def ledger(path):
        with h5py.File(path) as f:return f['energy_ledger/records'][:]
    rows=ledger(case/'ledger/energy-ledger.h5');old=ledger(off/'ledger/energy-ledger.h5')
    require(rows.dtype==old.dtype and rows.shape==old.shape,'v5 schema changed')
    fields={n:np.ascontiguousarray(rows[n]).tobytes()==np.ascontiguousarray(old[n]).tobytes() for n in rows.dtype.names if 'seconds' not in n}
    require(all(fields.values()),'Flag changed a scientific v5 field')
    save(root/'V5_LEDGER_BITWISE_REVIEW.json',{'checks_pass':True,'fields':fields,'timers_excluded':True})
    def a(n):return np.array(rows[n],dtype=LD)
    def cum(n):return np.r_[LD(0),np.cumsum(sums[n],dtype=LD)]
    delta_native=np.r_[LD(0),np.cumsum(sums['native_end_J']-sums['native_begin_J'],dtype=LD)]
    rf=cum('field_operator_residual_J');mid=cum('applied_current_mid_J')
    time=cum('time_electric_J')+cum('time_magnetic_J')+cum('current_time_pairing_J')
    cpml=cum('magnetic_cpml_J')+cum('electric_cpml_J')
    numerical=time+cum('cowan_operator_J')+cpml
    total=a('delta_particle_kinetic_box_J')+delta_native+cum('discrete_curl_exchange_J')-numerical-cum('inlet_overwrite_J')-cum('remap_J')-a('prescribed_boris_work_J')-a('frozen_constraint_work_J')-a('carrier_constraint_work_J')-a('boundary_source_work_J')
    coupling=a('particle_grid_work_J')-a('frozen_constraint_work_J')-a('carrier_constraint_work_J')-a('boundary_source_work_J')-mid
    reconstructed=rf+coupling+a('particle_work_balance_residual_J')+a('magnetic_rotation_delta_K_J')
    require(np.max(abs(total-reconstructed))<LD('1e-28'),'Independent total balance differs')
    spatial=(a('delta_field_energy_J')-a('inlet_overwrite_work_J')-a('moving_remap_field_change_J'))-(delta_native-cum('inlet_overwrite_J')-cum('remap_J'))
    flux=a('outward_field_energy_J')-cum('discrete_curl_exchange_J')
    bridge=spatial+flux+numerical
    require(np.max(abs(a('balance_residual_boris_J')-total-bridge))<LD('1e-28'),'Old/new bridge fails')
    phases={}
    quantities={'field_operator_residual':rf,'particle_grid_coupling_gap':coupling,'total_discrete_balance':total,
                'old_total_balance':a('balance_residual_boris_J'),'spatial_definition':spatial,'flux_definition':flux,
                'paired_time_term':time,'Cowan_operator_term':cum('cowan_operator_J'),'physical_cpml_term':cpml,
                'unapplied_current_mid_work':cum('unapplied_current_mid_J')}
    for label,i,j in [('whole',0,257),('released_no_exit',61,180),('zero_current_tail',244,257)]:
        phases[label]={'start_step':i,'end_step':j,**{n+'_J':float(v[j]-v[i]) for n,v in quantities.items()}}
    offreport=json.loads((off/'CASE_REVIEW.json').read_text());onreport=json.loads((case/'CASE_REVIEW.json').read_text())
    before=offreport['steps_40_240_seconds_per_step'];after=onreport['steps_40_240_seconds_per_step']
    result={'checks_pass':True,'rank_checks':reports,'corrupted_term_negative_control_rejected':True,'phases':phases,
      'max_abs_cumulative_field_operator_residual_J':float(np.max(abs(rf))),
      'max_abs_cumulative_total_discrete_balance_J':float(np.max(abs(total))),
      'performance':{'off_seconds_per_step':before,'on_seconds_per_step':after,'increment_percent':100*(after/before-1),
                     'scope':'One same-allocation pair, no profiling, full small-card output. Not a production or disabled-overhead estimate.'},
      'scope':'Operator identity + independent particle-grid coupling; explicit time/Cowan terms are numerical channels, not radiation or external work.',
      'physical_conservation_accepted':False}
    save(root/'DISCRETE_REVIEW.json',result)
    with (root/'discrete_values.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['step']+list(quantities));w.writerows([[i]+[str(v[i]) for v in quantities.values()] for i in range(258)])
    return result
