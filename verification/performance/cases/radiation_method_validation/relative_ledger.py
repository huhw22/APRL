"""Relative display of the current discrete ledger; no physical data changes.

One positive, fixed denominator is used throughout each completed case. Prefer
max(abs(delta_particle_kinetic_box_J)); zero particle change falls back to the
initial, then peak, native field energy. An all-zero case has no relative scale
and is rejected explicitly. Original signed joule quantities remain exportable.
"""
import csv,json
from pathlib import Path
import numpy as np

LD=np.longdouble

def construct(rows,sums,index,q):
    steps=np.asarray(rows['step'],dtype=np.int64)
    assert np.array_equal(steps,index) and index[0]==0
    field_steps=len(sums['native_begin_J'])
    def cumulative(name):return np.r_[LD(0),np.cumsum(sums[name],dtype=LD)][index]
    def a(name):return np.asarray(rows[name],dtype=LD)
    delta_k=a('delta_particle_kinetic_box_J')
    native=np.r_[LD(0),np.cumsum(sums['native_end_J']-sums['native_begin_J'],dtype=LD)][index]
    curl=cumulative('discrete_curl_exchange_J')
    time=cumulative('time_electric_J')+cumulative('time_magnetic_J')+cumulative('current_time_pairing_J')
    cowan=cumulative('cowan_operator_J');cpml=cumulative('magnetic_cpml_J')+cumulative('electric_cpml_J')
    inlet=cumulative('inlet_overwrite_J');remap=cumulative('remap_J')
    source=a('prescribed_boris_work_J')+a('boundary_source_work_J')+a('frozen_constraint_work_J')+a('carrier_constraint_work_J')
    # Derived directly from independent field/boundary/source channels, not R.
    reconstructed=-native-curl+time+cowan+cpml+inlet+remap+source
    values={
      'particle_change_J':delta_k,
      'field_boundary_source_prediction_J':reconstructed,
      'native_field_change_J':native,'outward_curl_exchange_J':curl,
      'inlet_energy_added_J':inlet,'remap_energy_added_J':remap,
      'prescribed_and_constraint_source_work_J':source,
      'time_pairing_channel_J':time,'Cowan_channel_J':cowan,'physical_CPML_channel_J':cpml,
      'total_residual_J':np.asarray(q['discrete_total_residual_J'],dtype=LD),
      'particle_grid_exchange_difference_J':np.asarray(q['particle_grid_coupling_gap_J'],dtype=LD),
      'field_operator_residual_J':np.asarray(q['field_operator_residual_J'],dtype=LD)}
    assert all(np.isfinite(v).all() for v in values.values())
    bound=2048*LD(np.finfo(float).eps)*np.maximum(sum(abs(v) for v in values.values()),LD('1e-100'))
    assert np.all(abs(delta_k-reconstructed-values['total_residual_J'])<=bound)
    initial_field=LD(sums['native_begin_J'][0])
    peak_field=max(np.max(abs(sums['native_begin_J'])),np.max(abs(sums['native_end_J'])))
    peak_k=np.max(abs(delta_k));scale=peak_k;policy='peak_absolute_particle_kinetic_change'
    if scale==0:
        scale=abs(initial_field);policy='initial_native_field_energy_zero_particle_change'
    if scale==0:
        scale=LD(peak_field);policy='peak_native_field_energy_zero_initial_field_and_particle_change'
    if not np.isfinite(scale) or scale<=0:
        raise ValueError('All energy scales are zero; a relative budget is undefined')
    pct={k.removesuffix('_J')+'_pct':100*v/scale for k,v in values.items()}
    metadata={'frame':'boost','unit':'percent','definition':'100 * signed_energy_quantity / fixed_denominator_J',
      'denominator_policy':policy,'fixed_denominator_J':float(scale),'fixed_denominator_decimal_J':str(scale),
      'particle_change_peak_absolute_J':float(peak_k),'particle_change_final_J':float(delta_k[-1]),
      'time_interval_steps':[int(index[0]),int(index[-1])],'all_field_steps':field_steps,
      'same_denominator_for_all_series':True,'signed_quantities_preserved':True,
      'normalization_is_not_a_physical_acceptance_criterion':True,
      'final_total_residual_pct_fixed_scale':float(pct['total_residual_pct'][-1]),
      'maximum_absolute_total_residual_pct_fixed_scale':float(np.max(abs(pct['total_residual_pct']))),
      'maximum_absolute_field_residual_pct_fixed_scale':float(np.max(abs(pct['field_operator_residual_pct']))),
      'final_total_residual_pct_absolute_final_particle_change':float(100*values['total_residual_J'][-1]/abs(delta_k[-1])) if delta_k[-1]!=0 else None,
      'maximum_absolute_total_residual_step':int(index[np.argmax(abs(values['total_residual_J']))]),
      'warning':'Particle energy change includes retained and removed members of the same cohort. Numerical channels are retained, not interpreted as radiation.'}
    return values,pct,metadata


def write_tables(out,index,absolute,pct,metadata):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for name,data in [('energy_budget.csv',pct),('energy_budget_absolute.csv',absolute)]:
        with (out/name).open('w') as f:
            writer=csv.writer(f);writer.writerow(['step']+list(data))
            writer.writerows([[int(step)]+[str(v[i]) for v in data.values()] for i,step in enumerate(index)])
    (out/'normalization.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def plot(out,index,pct,metadata,title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out)
    fig,axes=plt.subplots(3,1,figsize=(11,10),sharex=True,layout='constrained')
    axes[0].plot(index,pct['particle_change_pct'],color='#4b78a4',label='Particle kinetic-energy change')
    axes[0].plot(index,pct['field_boundary_source_prediction_pct'],'--',color='#d38b39',label='Independent sum of field, boundary and source channels')
    axes[0].set_title('Energy change and independent accounting')
    axes[1].plot(index,pct['total_residual_pct'],color='#4b78a4',label='Total accounting difference')
    axes[1].plot(index,pct['particle_grid_exchange_difference_pct'],'--',color='#d38b39',label='Particle-grid exchange difference')
    axes[1].set_title('Remaining difference relative to the same energy-change scale')
    axes[2].plot(index,pct['field_operator_residual_pct'],color='#714b86',label='Field operator accounting difference')
    axes[2].set_title('Field operator arithmetic consistency')
    for ax in axes:
        ax.axhline(0,color='#999999',linewidth=.7);ax.grid(alpha=.2);ax.set_ylabel('% of fixed D');ax.legend(fontsize=9)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useOffset=False)
    axes[-1].set_xlim(int(index[0]),int(index[-1]));axes[-1].set_xlabel('Maxwell step')
    ticks=list(axes[-1].get_xticks());ticks=[int(t) for t in ticks if index[0]<=t<=index[-1]]
    if ticks[-1]!=index[-1]:ticks[-1]=int(index[-1])
    axes[-1].set_xticks(ticks)
    label='D = max |particle kinetic-energy change|' if metadata['denominator_policy']=='peak_absolute_particle_kinetic_change' else 'D = native field energy (zero particle change)'
    fig.suptitle(title+'\n'+label+'; boost frame; one denominator for the complete displayed interval',fontsize=11)
    fig.savefig(out/'energy_budget.png',dpi=180);fig.savefig(out/'energy_budget.pdf');plt.close(fig)
