"""Closed-result report for an independent longitudinal coupling investigation."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import sha,save,require
ME=9.1093837015e-31;C=299792458.;N=192e5;DZ=1.2543262415359034e-6

def load(p):return json.loads(p.read_text())
def rows(p):return np.genfromtxt(p,delimiter=',',names=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();r=w/'results/longitudinal_coupling_diagnosis_20261102';out=r/'analysis';out.mkdir(exist_ok=True);checks={}
 for stage in ['r3','force_r2','spline_r3','spline_r4','fluid_r1','remainder_r1']:
  d=r/stage;require((d/'EXIT_CODE').read_text().strip()=='0','Not closed '+stage)
  for f,h in load(d/'FROZEN.json').items():require(sha(f)==h,'Frozen mismatch '+f)
  checks[stage]='closed and all frozen inputs unchanged'
 doc=w/'results/lab_energy_balance_20261102/yee_sbp_r1/PREEXISTING_DOC_HASHES.json'
 protected=load(doc)
 # Historical manifest has repository-relative names.
 for f,h in protected.items():
  target=Path(f) if Path(f).is_absolute() else w/'source/aprl'/f
  require(sha(target)==h,'Preexisting document changed '+str(target))
 binary=w/'build/source-coupling-exit-r3-20261031/aprl';require(sha(binary)=='865ebef5b79b6b8670d40c81080f15f916ab1938033afcabd050bcd3c8198752','Production binary changed')
 checks['production_binary']='unchanged';checks['preexisting_dirty_documents']='unchanged'
 fref=load(r/'fluid_r1/analysis/RESULT.json');require(fref['reference_converged'],'Reference not converged')
 base=rows(r/'force_r2/longitudinal_raw_current.csv');scalar=rows(r/'spline_r3/linear_baseline.csv')
 replay={}
 for x,y in [('particle_delta_Kprime_J','particle_delta_Kprime_J'),('particle_delta_cPzprime_J','particle_delta_cPzprime_J')]:
  # Production control CSV uses the same first physical columns.
  bx=x if x in base.dtype.names else base.dtype.names[1 if 'Kprime' in x else 2]
  error=float(np.max(np.abs(base[bx]-scalar[y])));replay[y]=error;require(error<1e-20,'Scalar/3D replay mismatch')
 summary=dict(reference=fref,scalar_original_3D_replay_max_absolute_J=replay,checks=checks,production_modified=False,full_3D_physical_acceptance=False)
 table=[]
 for stage in ['spline_r3','spline_r4','remainder_r1']:
  for name,item in load(r/stage/'RESULT.json')['cases'].items():
   d=rows(r/stage/(name+'.csv'));conf=load(r/stage/(name+'.json'))['configuration'];dt=DZ/C/conf['args'][2];D=float(np.max(np.abs(d['particle_delta_Kprime_J'])));P=float(np.max(np.abs(d['particle_delta_cPzprime_J'])));row=dict(stage=stage,case=name,degree=item['degree'],compensation=item['compensation_order'],cells=item['cells'],sheet_pairs=item['sheet_pairs'],duration_s=float(d['step'][-1]*dt),max_cPz_residual_J=P,max_mean_uz_error=P/(N*ME*C*C),max_energy_residual_over_fixed_particle_change=float(np.max(np.abs(d['energy_residual_J']))/D),fixed_particle_change_J=D,max_energy_residual_J=float(np.max(np.abs(d['energy_residual_J']))),max_gauss=item['max_gauss_L1_relative'])
   if 'spatial_impulse_gap_J' in d.dtype.names and stage=='remainder_r1':
    err=d['particle_delta_cPzprime_J']-sum(d[k] for k in ['spatial_impulse_gap_J','time_impulse_gap_J','endpoint_charge_impulse_J','iteration_impulse_gap_J']);require(np.max(abs(err))<1e-23,'Impulse split mismatch')
    for k in ['spatial_impulse_gap_J','time_impulse_gap_J','endpoint_charge_impulse_J','iteration_impulse_gap_J']:row['max_'+k]=float(np.max(abs(d[k])));row['final_'+k]=float(d[k][-1])
   table.append(row)
 summary['cases']=table;save(out/'RESULT.json',summary)
 keys=list(dict.fromkeys(k for row in table for k in row))
 with (out/'comparison.csv').open('w') as f:writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(table)
 plt.rcParams.update({'font.size':10,'figure.dpi':140,'savefig.dpi':180})
 def finish(fig,name):
  fig.text(.5,.01,'Independent 1-D model; original Cowan propagator unchanged; not full-program acceptance',ha='center',fontsize=8,color='#555555');fig.tight_layout(rect=(0,.035,1,1));fig.savefig(out/(name+'.png'));fig.savefig(out/(name+'.pdf'));plt.close(fig)
 methods=[('linear_stationary','Linear','#d97706'),('cubic_compensated_stationary','Cubic + compensation','#2563eb'),('quintic_compensated_stationary','Quintic + compensation','#159059')]
 fig,axes=plt.subplots(1,2,figsize=(11,4.3));D=float(np.max(abs(rows(r/'spline_r3/linear_stationary.csv')['particle_delta_Kprime_J'])))
 save(out/'normalization.json',dict(energy_plot_fixed_denominator_J=D,definition='max absolute same-cohort particle kinetic-energy change in linear stationary baseline; same fixed denominator for all curves',original_J_retained=True))
 for name,label,color in methods:
  d=rows(r/'spline_r3'/(name+'.csv'));t=d['step']*DZ/C*1e12;axes[0].plot(t,d['mean_delta_uz'],label=label,color=color);axes[1].plot(t,100*d['energy_residual_J']/D,color=color)
 axes[0].axhline(0,color='black',lw=.8);axes[0].set_yscale('symlog',linthresh=1e-10);axes[0].set_ylabel(r'$\Delta\langle u_z\rangle$ (physical reference = 0)');axes[0].legend(fontsize=8);axes[0].set_title('Longitudinal spurious impulse');axes[1].set_ylabel('Energy residual / fixed baseline change (%)');axes[1].set_title('Energy exchange remains paired')
 for ax in axes:ax.set_xlabel('Time (ps)');ax.grid(alpha=.2)
 finish(fig,'01_longitudinal_impulse')
 fig,axes=plt.subplots(1,3,figsize=(13,4.3));labels=[('analytic_linear','Linear','#d97706'),('analytic_cubic_continuum','Cubic + continuum compensation','#2563eb'),('analytic_quintic_continuum','Quintic + continuum compensation','#159059')]
 for name,label,color in labels:
  d=rows(r/'fluid_r1/analysis'/(name+'.csv'));mask=(d['id']%2==1);d=d[mask];x=d['initial_z_m']*1e6;scale=fref['velocity_scale']
  if name=='analytic_linear':axes[0].plot(x,d['reference_uz'],color='black',lw=1.8,label='Independent cold fluid')
  axes[0].plot(x,d['PIC_uz'],color=color,lw=1,label=label);target=axes[1] if name=='analytic_linear' else axes[2];target.plot(x,100*(d['PIC_uz']-d['reference_uz'])/scale,color=color,label=label)
 axes[0].set_ylabel(r'Final $u_z$');axes[0].set_title('Individual longitudinal momentum');axes[0].legend(fontsize=7)
 for ax in axes[1:]:ax.set_ylabel('Error / reference peak change (%)');ax.axhline(0,color='black',lw=.7)
 axes[1].set_title('Linear assignment error');axes[2].set_title('Higher longitudinal order');axes[2].legend(fontsize=7)
 for ax in axes:ax.set_xlabel('Initial longitudinal position (um)');ax.grid(alpha=.2)
 finish(fig,'02_individual_reference')
 fig,axes=plt.subplots(1,2,figsize=(11,4.3))
 for name,label,color in labels:
  d=rows(r/'spline_r4'/(name+'.csv'));j=load(r/'spline_r4'/(name+'.json'))['result'];A=j['continuum_linear_initial_E_amplitude_V_per_m'];t=d['step']*DZ/C;err=np.abs((d['fundamental_E_real_V_per_m']+1j*d['fundamental_E_imag_V_per_m'])/A-1j*np.cos(j['linear_plasma_omega_rad_s']*t));axes[0].semilogy(t*1e12,np.maximum(err,1e-18),color=color,label=label)
 for name,label,color in [('linear_stationary','Linear','#d97706'),('quintic_stationary','Quintic + baseline compensation','#159059'),('quintic_stationary_continuum','Quintic + continuum compensation','#7c3aed')]:
  d=rows(r/'remainder_r1'/(name+'.csv'));axes[1].plot(d['step']*DZ/C*1e12,d['mean_delta_uz'],color=color,label=label)
 axes[0].set_ylabel('Complex field-mode error / theoretical amplitude');axes[0].set_title('Weak-mode formula');axes[1].set_ylabel(r'$\Delta\langle u_z\rangle$ (reference = 0)');axes[1].set_yscale('symlog',linthresh=1e-10);axes[1].set_title('Longer, stronger control: residual remains')
 for ax in axes:ax.set_xlabel('Time (ps)');ax.grid(alpha=.2);ax.legend(fontsize=7)
 finish(fig,'03_formula_and_longer_control')
 print(json.dumps(dict(report=str(out),reference=fref['comparisons'],replay=replay),indent=2))
if __name__=='__main__':main()
