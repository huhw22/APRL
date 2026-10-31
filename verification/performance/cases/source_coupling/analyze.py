#!/usr/bin/env python3
"""Closed small-PIC results, independent cold-sphere reference, paired cost."""
import argparse,hashlib,json,re,shutil
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-source-coupling-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
campaign=a.campaign.resolve();root=a.output.resolve();root.mkdir(exist_ok=False);(root/'figures').mkdir()
complete=json.loads((campaign/'COMPLETE.json').read_text());assert complete['checks_passed']
plan=json.loads((campaign/'plan.json').read_text());results=complete['cases'];rows={};e=1.602176634e-19;m=9.1093837139e-31;c=299792458.;eps=8.8541878128e-12
def save(name,x):(root/name).write_text(json.dumps(x,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z))+'\n')
for case in plan['cases']:
 name=case['name'];assert results[name]['complete']
 if not case['audit']:continue
 x=np.genfromtxt(campaign/name/'history.csv',delimiter=',',names=True,dtype=np.longdouble)
 assert np.array_equal(x['step'],np.arange(1,case['steps']+1));assert all(np.isfinite(x[k]).all() for k in x.dtype.names)
 rows[name]=x
 assert results[name]['relative_continuity_residual_max']<5e-12
assert results['vacuum-baseline']['field_hash']==results['vacuum-coupled']['field_hash']
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'axes.grid':True,'grid.alpha':.18})
blue='#28639A';orange='#CB7723';pink='#AF5377';gray='#66717B'
exchange={};normalization={};fig,axes=plt.subplots(3,2,figsize=(11.8,10.1),layout='constrained')
for row,kind in enumerate(['sphere','drift','undulator']):
 pair={mode:rows[kind+'-'+mode] for mode in ['baseline','coupled']}
 D=max(float(np.max(abs(x['delta_particle_box_J']))) for x in pair.values());assert D>0
 Dterminal=max(float(abs(x['delta_particle_box_J'][-1])) for x in pair.values());assert Dterminal>0
 Dgrid=max(float(abs(x['delta_particle_box_J'][-1]-x['prescribed_work_J'][-1])) for x in pair.values())
 Dgrid=max(Dgrid,max(float(abs(x['mesh_mid_J'][-1])) for x in pair.values()));assert Dgrid>0
 normalization[kind]={'primary_fixed_D_J':D,'definition':'maximum absolute stable-cohort delta K boost over both completed paired histories','secondary_fixed_terminal_D_J':Dterminal,'secondary_particle_grid_work_D_J':Dgrid,'secondary_particle_grid_work_definition':'max of absolute terminal (actual particle delta K minus prescribed work) and mesh midpoint work across the pair; common fixed scale, never time-dependent','shared_between_baseline_and_candidate':True}
 exchange[kind]={}
 for mode,color,style in [('baseline',blue,'--'),('coupled',orange,'-')]:
  x=pair[mode];actual=x['coupling_gap_J']+x['pusher_residual_J'];operator=x['field_operator_residual_J'];strict=actual+operator
  for col,denominator in enumerate([D,Dgrid]):
   axes[row,col].plot(x['time_s']*1e15,100*actual/denominator,color=color,ls=style,label=mode+' · actual ΔK − Wext − J·E',lw=1.8)
  exchange[kind][mode]={'terminal_actual_exchange_gap_J':float(actual[-1]),'maximum_actual_exchange_gap_J':float(np.max(abs(actual))),
   'terminal_grid_work_minus_mesh_work_J':float(x['coupling_gap_J'][-1]),'terminal_local_pusher_roundoff_J':float(x['pusher_residual_J'][-1]),
   'maximum_local_pusher_roundoff_J':float(np.max(abs(x['pusher_residual_J']))),'terminal_field_operator_residual_J':float(operator[-1]),
   'maximum_field_operator_residual_J':float(np.max(abs(operator))),'terminal_strict_total_residual_J':float(strict[-1]),
   'terminal_gap_percent_of_common_peak_D':100*float(actual[-1])/D,'terminal_gap_percent_of_common_terminal_D':100*float(actual[-1])/Dterminal,
   'terminal_gap_percent_of_common_particle_grid_D':100*float(actual[-1])/Dgrid,
   'terminal_actual_particle_grid_work_J':float(x['delta_particle_box_J'][-1]-x['prescribed_work_J'][-1]),'terminal_mesh_midpoint_work_J':float(x['mesh_mid_J'][-1]),
   'final_particle_delta_box_J':float(x['delta_particle_box_J'][-1]),'mean_lab_loss_eV':float(x['mean_lab_loss_eV'][-1]),'std_lab_loss_eV':float(x['std_lab_loss_eV'][-1]),
   'mean_uz_final':float(x['mean_uz'][-1]),'max_outer_iterations':results[kind+'-'+mode]['maximum_outer_iterations'],'mean_outer_iterations':results[kind+'-'+mode]['mean_outer_iterations']}
  shutil.copy2(campaign/(kind+'-'+mode)/'history.csv',root/(kind+'-'+mode+'-absolute.csv'))
  with (root/(kind+'-'+mode+'-relative.csv')).open('w') as f:
   f.write('step,time_s,actual_exchange_gap_percent,grid_exchange_gap_percent,pusher_residual_percent,field_operator_residual_percent,strict_total_residual_percent\n')
   for i in range(len(x)):f.write(','.join([str(int(x['step'][i])),format(x['time_s'][i],'.18e'),*[format(100*v[i]/D,'.18e') for v in [actual,x['coupling_gap_J'],x['pusher_residual_J'],operator,strict]]])+'\n')
 old=exchange[kind]['baseline'];new=exchange[kind]['coupled'];exchange[kind]['maximum_gap_reduction_factor']=old['maximum_actual_exchange_gap_J']/max(new['maximum_actual_exchange_gap_J'],1e-300)
 for col in [0,1]:
  ax=axes[row,col];ax.axhline(0,color=gray,lw=.7);ax.set_xlabel('Computational time (fs)');ax.set_ylabel('Exchange defect (% of fixed scale)');ax.set_title(kind+' · '+('common peak |ΔK| scale' if col==0 else 'common terminal particle-grid work scale'));ax.legend(frameon=False,fontsize=8)
fig.suptitle('Dynamic 3D PIC · original Cowan/CPML and charge-conserving deposition · actual particle energy',fontsize=12)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('exchange.'+ext),dpi=180)
plt.close(fig);save('EXCHANGE.json',exchange);save('NORMALIZATION.json',normalization)

# Independent nonrelativistic uniform-sphere expansion in dimensionless time.
R=60e-6;N=1e8;kappa=N*e*e/(4*np.pi*eps*m*R**3)
times=np.unique(np.concatenate([np.asarray(rows[n]['time_s'],float) for n in ['sphere-baseline','sphere-coupled','sphere-halfdt-baseline','sphere-halfdt-coupled']]))
reference=solve_ivp(lambda tau,y:[y[1],1/y[0]**2],[0,float(times[-1]*np.sqrt(kappa))],[1.,0.],rtol=1e-12,atol=1e-14,dense_output=True)
assert reference.success
def theory(t):return reference.sol(np.asarray(t,float)*np.sqrt(kappa))
nbase=results['sphere-baseline']['particles']//8;rms0=R*np.sqrt(np.mean(((np.arange(nbase)+.5)/nbase)**(2/3)))
sphere={};fig,axes=plt.subplots(1,2,figsize=(11.3,4.6),layout='constrained');aa,bb=theory(times)
axes[0].plot(times*1e15,100*(aa-1),color='#25282C',lw=2,label='Independent uniform-sphere ODE')
axes[1].plot(times*1e15,.5*N*m*rms0*rms0*kappa*bb*bb*1e12,color='#25282C',lw=2,label='Independent uniform-sphere ODE')
for name,color,style,marker in [('sphere-baseline',blue,'--','o'),('sphere-coupled',orange,'-','s'),('sphere-halfdt-baseline',blue,':','^'),('sphere-halfdt-coupled',orange,'-.','D')]:
 x=rows[name];a1,b1=theory(x['time_s']);expansion=np.asarray(x['rms_radius_m'],float)/rms0-1;referenceK=.5*N*m*rms0*rms0*kappa*b1*b1
 sphere[name]={'time_s':float(x['time_s'][-1]),'final_radius_growth_percent':100*float(expansion[-1]),'theory_final_radius_growth_percent':100*float(a1[-1]-1),
  'final_relative_error_in_radius_growth_percent':100*(float(expansion[-1])/(a1[-1]-1)-1),'final_delta_K_J':float(x['delta_particle_box_J'][-1]),'theory_final_delta_K_J':float(referenceK[-1]),
  'final_relative_error_in_delta_K_percent':100*(float(x['delta_particle_box_J'][-1])/referenceK[-1]-1)}
 label=name.replace('halfdt','dt,dz/2')
 axes[0].plot(x['time_s']*1e15,100*expansion,color=color,ls=style,lw=1.6,marker=marker,markevery=max(1,len(x)//6),ms=3,label=label)
 axes[1].plot(x['time_s']*1e15,x['delta_particle_box_J']*1e12,color=color,ls=style,lw=1.6,marker=marker,markevery=max(1,len(x)//6),ms=3,label=label)
for ax in axes:ax.set_xlabel('Computational/rest-frame time (fs)');ax.legend(frameon=False,fontsize=8)
axes[0].set_ylabel('RMS radius growth (%)');axes[1].set_ylabel('Particle kinetic-energy gain (pJ)')
fig.suptitle('Cold uniform sphere · matched dt,dz refinement · same domain and physical end time',fontsize=12)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('sphere.'+ext),dpi=180)
plt.close(fig)
save('SPHERE_REFERENCE.json',{'ODE':'a_double_dot=kappa/a^2, a(0)=1, a_dot(0)=0; kappa=N e^2/(4 pi eps0 m R^3)','kappa_s_minus2':kappa,'quiet_start_rms_radius_m':rms0,'maximum_theory_surface_speed_over_c':float(R*np.sqrt(kappa)*np.max(bb)/c),'classical_reference_not_a_radiation_model':True,'comparison':sphere})

external=np.genfromtxt(campaign/'external.csv',delimiter=',',names=True,dtype=None,encoding='utf8');external_summary={};fig,axes=plt.subplots(1,2,figsize=(10.8,4.5),layout='constrained')
for method,color,marker in [('original_boris',blue,'o'),('mean_velocity',orange,'s')]:
 selected=external[external['method']==method];external_summary[method]={'maximum_sampled_absolute_lab_energy_error_eV':float(np.max(selected['max_lab_energy_change_eV'])),'maximum_terminal_absolute_lab_energy_error_eV':float(np.max(abs(selected['final_lab_energy_change_eV'])))}
 xx=[];yy=[]
 for n in [5,10,20]:xx.append(n);yy.append(np.max(selected['max_lab_energy_change_eV'][selected['substeps']==n]))
 axes[0].semilogy(xx,yy,color=color,marker=marker,label=method,lw=1.7)
 for sample in range(3):
  q=selected[selected['sample']==sample];q=q[np.argsort(q['substeps'])];u=np.array([q['lab_ux'],q['lab_uy'],q['lab_uz']]).T
  error=np.linalg.norm(u[:-1]-u[-1],axis=1);axes[1].loglog(q['substeps'][:-1],np.maximum(error,1e-16),color=color,marker=marker,ls=['-','--',':'][sample],label=method+f' sample {sample}')
for ax in axes:ax.set_xlabel('Particle substeps per unchanged field step');ax.legend(frameon=False,fontsize=8)
axes[0].set_ylabel('Maximum sampled |lab ΔK| (eV / electron)');axes[0].set_title('Static lab undulator: exact ΔK = 0')
axes[1].set_ylabel('|u_end − u_end,20|');axes[1].set_title('Exit momentum convergence; 20 is numerical reference')
for ext in ['png','pdf']:fig.savefig(root/'figures'/('external.'+ext),dpi=180)
plt.close(fig);shutil.copy2(campaign/'external.csv',root/'external.csv');save('EXTERNAL.json',external_summary)

timing={};fig,axes=plt.subplots(2,2,figsize=(10.5,8),layout='constrained')
for row,ledger in enumerate(['','ledger-']):
 for index,size in enumerate(['field-heavy','particle-heavy']):
  prefix='timing-'+ledger+size
  if prefix+'-r0-baseline' not in results:continue
  pairs=[]
  for repeat in range(3):
   b=results[f'{prefix}-r{repeat}-baseline'];q=results[f'{prefix}-r{repeat}-coupled']
   pairs.append({'repeat':repeat,'baseline_seconds_per_step':b['seconds_per_step'],'candidate_seconds_per_step':q['seconds_per_step'],'percent_cost_increase':100*(q['seconds_per_step']/b['seconds_per_step']-1),'baseline_RSS_MiB':b['peak_RSS_KiB']/1024,'candidate_RSS_MiB':q['peak_RSS_KiB']/1024,'candidate_mean_outer_iterations':q['mean_outer_iterations']})
  timing[ledger+size]={'particles':results[prefix+'-r0-baseline']['particles'],'paired_runs':pairs,'median_percent_cost_increase':float(np.median([z['percent_cost_increase'] for z in pairs])),'median_RSS_increase_MiB':float(np.median([z['candidate_RSS_MiB']-z['baseline_RSS_MiB'] for z in pairs])),'confidence_interval_estimated':False}
  for ax,keys in [(axes[row,0],['baseline_seconds_per_step','candidate_seconds_per_step']),(axes[row,1],['baseline_RSS_MiB','candidate_RSS_MiB'])]:
   for j,(key,color) in enumerate(zip(keys,[blue,orange])):
    vals=[z[key] for z in pairs];ax.scatter(np.full(3,index+.13*(j-.5)),vals,color=color,marker='o' if j==0 else 's',s=30,label=('baseline' if j==0 else 'coupled') if index==0 else None)
 for ax in axes[row]:ax.set_xticks([0,1],['4,096 particles','131,072 particles']);ax.legend(frameon=False);ax.set_ylim(bottom=0)
 axes[row,0].set_ylabel('Physical update time / step (s)');axes[row,1].set_ylabel('Peak process RSS (MiB)')
 for ax in axes[row]:ax.set_title('Strict energy audit '+('off' if row==0 else 'on'))
fig.suptitle('Same 1 MPI × 8 cores · three alternating-order pairs · audit off and on reported separately',fontsize=12)
for ext in ['png','pdf']:fig.savefig(root/'figures'/('cost.'+ext),dpi=180)
plt.close(fig);save('TIMING.json',timing)
save('ANALYSIS.json',{'complete':True,'scientific_cases_complete':True,'source_coupling_production_accepted':False,'figures_png_pdf_generated':True,'independent_particle_energy_defect_includes_local_pusher_roundoff':True})
print(json.dumps({'exchange':exchange,'sphere':sphere,'cost':timing,'external':external_summary},indent=2))
