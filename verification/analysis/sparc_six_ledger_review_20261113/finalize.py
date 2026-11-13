"""Deliver completed six-segment ledgers without altering production results."""
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD=np.longdouble
CASES=['energy-binomial3','energy-momentum']
TITLES=['19743: energy matching + three-point kernel','19744: energy matching + momentum matching']

def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def table(p,values):
    with Path(p).open('w') as f:
        w=csv.writer(f);w.writerow(list(values))
        for row in zip(*values.values()):w.writerow([str(x) for x in row])

def main(root):
    root=Path(root);parent=root.parent;answers={};checks={}
    assert read(root/'lab-v2/LAB_SUMMARY.json')['analysis_complete']
    assert read(root/'whole-v2/WHOLE_SUMMARY.json')['complete']
    fig,axes=plt.subplots(4,2,figsize=(14,13),sharex=True,layout='constrained')
    fig2,ax2=plt.subplots(2,2,figsize=(13,8),sharex=True,layout='constrained')
    for col,case in enumerate(CASES):
        b=read(root/'boost'/case/'BOOST_RESULT.json');w=read(root/'whole-v2'/case/'WHOLE_RESULT.json');l=read(root/'lab-v2'/case/'LAB_RESULT.json')
        with np.load(root/'boost'/case/'boost_cache.npz') as a:
            records=a['records'];step=records['step'];n=len(step)
            v={k:a[k].copy() for k in a.files if k.startswith('sum_')}
            cum=lambda k:np.r_[LD(0),np.cumsum(v['sum_'+k],dtype=LD)]
            channels={k:records[k].astype(LD) for k in ['step','time_box_s','delta_particle_kinetic_box_J','delta_particle_kinetic_lab_J',
                'delta_particle_cPz_box_J','prescribed_boris_work_J','boundary_source_work_J','frozen_constraint_work_J','carrier_constraint_work_J',
                'particle_grid_work_J','particle_work_balance_residual_J','magnetic_rotation_delta_K_J','released_charge_C',
                'active_macroparticles','removed_macroparticles']}
            for k in ['native_end_J','native_begin_J','discrete_curl_exchange_J','time_electric_J','time_magnetic_J',
                'current_time_pairing_J','cowan_operator_J','electric_cpml_J','magnetic_cpml_J','inlet_overwrite_J',
                'remap_J','applied_current_mid_J','unapplied_current_mid_J','field_operator_residual_J']:
                channels['cumulative_'+k]=cum(k)
            channels['native_field_change_J']=np.r_[LD(0),np.cumsum(v['sum_native_end_J']-v['sum_native_begin_J'],dtype=LD)]
            # Absolute stores are not cumulative sums; replace the two helper
            # stock columns with their correctly sampled begin/end values.
            channels.pop('cumulative_native_end_J');channels.pop('cumulative_native_begin_J')
            table(root/'boost'/case/'all_recorded_components.csv',channels)
            absolute=np.genfromtxt(root/'boost'/case/'energy_budget_absolute.csv',delimiter=',',names=True,dtype=LD)
            relative=np.genfromtxt(root/'boost'/case/'energy_budget.csv',delimiter=',',names=True,dtype=LD)
            D=LD(b['normalization']['fixed_denominator_decimal_J'])
            pct_error=max(float(np.max(abs(relative[k.removesuffix('_J')+'_pct']*D/100-absolute[k]))) for k in absolute.dtype.names if k!='step')
            assert pct_error<1e-17
            kprime=records['delta_particle_kinetic_box_J'].astype(LD);pz=records['delta_particle_cPz_box_J'].astype(LD)
            klab=records['delta_particle_kinetic_lab_J'].astype(LD);G=LD(records['boost_gamma'][0]);beta=LD(records['boost_beta'][0])
            for key,label,style in [('particle_change_J','Particle kinetic change','-'),('field_boundary_source_prediction_J','Independent field/source prediction','--')]:
                axes[0,col].plot(step,absolute[key]*1e6,style,label=label)
            for key,label,style in [('total_residual_J','Total difference','-'),('particle_grid_exchange_difference_J','Particle/current gap','--'),('field_operator_residual_J','Field arithmetic',':')]:
                axes[1,col].plot(step,absolute[key]*1e12,style,label=label)
        with np.load(root/'whole-v2'/case/'whole_cache.npz') as a:
            cc=lambda key:np.r_[LD(0),np.cumsum(a['fixed_end_'+key],dtype=LD)]
            egap=G*(kprime-cc('WJ_particle_J'));pgap=G*beta*(pz+cc('PJ_particle_J'))
            gap=a['particle_vs_deposited_current_Lorentz_gap_J'];split_error=float(np.max(abs(egap+pgap-gap)))
            assert split_error<1e-14
            table(root/'whole-v2'/case/'particle_energy_momentum_split.csv',dict(step=step,
                Gamma_particle_energy_minus_deposited_work_J=egap,
                Gamma_beta_particle_momentum_plus_deposited_impulse_J=pgap,
                combined_lab_gap_J=gap))
            for key,label,style in [('particle_kinetic_change_lab_J','Lab particle kinetic change','-'),
                ('field_boundary_source_prediction_lab_J','Independent field/source prediction','--')]:
                axes[2,col].plot(step,a[key]*1e6,style,label=label)
            for key,label,style in [('full_window_Lorentz_accounting_difference_J','Lab numerical difference','-'),
                ('particle_vs_deposited_current_Lorentz_gap_J','Particle/current Lorentz gap','--'),
                ('field_Lorentz_identity_residual_J','Field arithmetic',':')]:
                axes[3,col].plot(step,a[key]*1e6,style,label=label)
            energy_momentum=dict(energy_part_J=float(egap[-1]),momentum_part_J=float(pgap[-1]),
                combined_gap_J=float(gap[-1]),contains_external_and_constraint_contributions=True,
                finer_force_split_measured=False)
        labrows=list(csv.DictReader((root/'lab-v2'/case/'lab_energy_absolute.csv').open()))
        assert len(labrows)==7
        known=[r for r in labrows if r['physical_lab_sum_J']!=''];assert len(known)==6
        assert labrows[-1]['particle_time_caps_J']=='' and labrows[-1]['physical_lab_sum_J']==''
        z=[float(r['z_lab_m']) for r in known]
        for key,label,style in [('particle_net_transport_J','Particle net transport incl. sides','-'),
            ('independent_field_gain_J','Field gain: ports + sides + caps','--'),('physical_lab_sum_J','Physical sum','-')]:
            ax2[0,col].plot(z,[float(r[key])*1e6 for r in known],style,marker='o',label=label)
        for key,label in [('field_ports_net_J','Field exit - entrance'),('field_transverse_net_J','Four field sides'),
            ('field_time_caps_net_J','Late - early field cap')]:
            ax2[1,col].plot([float(r['z_lab_m']) for r in labrows],[float(r[key])*1e6 for r in labrows],marker='o',label=label)
        ax2[0,col].set_title(TITLES[col]);axes[0,col].set_title(TITLES[col])
        for row,unit in enumerate(['µJ','pJ','µJ','µJ']):
            axes[row,col].set_ylabel(['Boost energy change','Boost accounting difference','Lab energy on boost slices','Lab accounting difference'][row]+' ('+unit+')')
            axes[row,col].grid(alpha=.2);axes[row,col].legend(fontsize=7);axes[row,col].axhline(0,color='grey',lw=.6)
        axes[3,col].set_xlabel('Maxwell step')
        for ax in ax2[:,col]:
            ax.set(xlabel='Laboratory exit plane z (m)',ylabel='Signed energy (µJ)',xlim=(0,15.3));ax.grid(alpha=.2);ax.legend(fontsize=8);ax.axhline(0,color='grey',lw=.6)
        ax2[0,col].axvspan(12.8,15.3,color='grey',alpha=.12)
        ax2[0,col].text(.91,.12,'p6: unknown\nparticle cap',ha='center',transform=ax2[0,col].transAxes,fontsize=8)
        wc=w['final_components_J'];last=known[-1]
        answers[case]=dict(job=[19743,19744][col],steps=26923,
            boost_particle_change_J=b['final_boost_particle_change_J'],boost_residual_J=b['final']['discrete_total_residual_J'],
            boost_residual_percent_of_abs_final_particle_change=b['normalization']['final_total_residual_pct_absolute_final_particle_change'],
            lab_particle_change_full_cohort_J=wc['particle_kinetic_change_lab_J'],lab_slice_difference_J=wc['full_window_Lorentz_accounting_difference_J'],
            lab_slice_difference_percent_of_abs_final_particle_change=100*wc['full_window_Lorentz_accounting_difference_J']/abs(wc['particle_kinetic_change_lab_J']),
            p5_particle_net_J=float(last['particle_net_transport_J']),p5_independent_field_gain_J=float(last['independent_field_gain_J']),
            p5_complete_physical_sum_J=l['last_complete_residual_J'],p6_missing_particle_cap_count=l['coverage'][-1]['unmatched_particle_ID_counts'],
            p6_physical_sum_J=None,lab_energy_momentum_gap_split=energy_momentum,
            physical_conservation_accepted=False)
        checks[case]=dict(complete_steps=n-1,boost_percent_back_conversion_max_error_J=pct_error,
            energy_momentum_gap_split_max_error_J=split_error,
            particle_Lorentz_identity_max_error_J=b['max_particle_Lorentz_identity_error_J'],
            lab_field_identity_max_cumulative_error_J=w['max_cumulative_Lorentz_field_identity_J'],
            physical_worldtube_surface_coverage_max_error=l['maximum_surface_coverage_error'],
            p0_p5_ID_boundary_counts_closed=True,p6_missing_cap_not_replaced_by_zero=True)
    fig.suptitle('Six-segment full-duration energy accounting | both coordinate frames\nLab panel uses energy/momentum charges on boost-time slices; fixed-lab tube shown separately')
    fig.savefig(root/'both_frames.png',dpi=175);fig.savefig(root/'both_frames.pdf');plt.close(fig)
    fig2.suptitle('Independent fixed-laboratory worldtube | full common retarded-time gate, all frequencies\nPhysical sum = particle transport + field gain; p6 has unmeasured late particle-cap energy')
    fig2.savefig(root/'lab_physical_absolute.png',dpi=175);fig2.savefig(root/'lab_physical_absolute.pdf');plt.close(fig2)
    frozen=read(parent/'FROZEN.json');mismatch=[p for p,h in frozen.items() if digest(p)!=h]
    assert not mismatch,mismatch
    save(root/'VALIDATION.json',dict(passed=True,scope='Analysis reconstruction and coverage checks, not physical acceptance',
        cases=checks,production_frozen_files_sha256_checked=len(frozen),production_frozen_mismatches=mismatch,
        geometry=read(root/'lab-v2/GEOMETRY_CHECK.json'),raw_field_streams_finite_scanned=True,
        original_raw_outputs_not_rewritten=True,PDF_exported_from_same_figure=True,PDF_separately_rendered=False))
    save(root/'RESULT.json',dict(analysis_complete=True,physical_conservation_accepted=False,
        no_new_simulation=True,original_loading='Historical correlated phase-clone input, not the new independent-random source.',cases=answers))
    shutil.copy2(parent.parent/'lab_worldtube_balance_20261112_r1/THEORY_LAB_COMPLETE.md',root/'THEORY_LAB_COMPLETE_previous.md')
    print(json.dumps(answers,indent=2),flush=True)

if __name__=='__main__':main(sys.argv[1])
