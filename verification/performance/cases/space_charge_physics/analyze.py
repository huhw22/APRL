#!/usr/bin/env python3
"""Compare only closed PIC outputs with independent reference worldlines."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-space-charge-pic-mpl')
import argparse,csv,json,sys
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import Sphere,G_BEAM,B_BEAM,c,m_e,e
from prepare import save
try:
    import report_discrete
except ModuleNotFoundError as error:
    if error.name!='report_discrete':raise
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'discrete_energy_ledger'))
    import report as report_discrete
try:
    from relative_ledger import construct,write_tables,plot as plot_relative
except ModuleNotFoundError as error:
    if error.name!='relative_ledger':raise
    sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'analysis/relative_energy_ledger'))
    from relative_ledger import construct,write_tables,plot as plot_relative

LD=np.longdouble
def read_closed(path,group):
    with h5py.File(path,'r') as h:
        g=h[group];assert int(g['complete'][()])==1,path
        data=g['records'][:];assert int(g['committed_records'][()])==len(data),path
    for n in data.dtype.names:
        if data[n].dtype.kind in 'fc':assert np.isfinite(data[n]).all(),(path,n)
    return data

def rest_u(lab):
    u=np.asarray(lab,dtype=LD);g=np.sqrt(1+np.sum(u*u,axis=1));G=LD(G_BEAM);gb=np.sqrt((G-1)*(G+1))
    plus=(g+u[:,2])/(G+gb);minus=(1+np.sum(u[:,:2]**2,axis=1))/plus
    return np.c_[u[:,:2],(plus-minus)/2]
def kinetic_ev(u):
    s=np.sum(u*u,axis=1);return (LD(m_e)*LD(c)**2/LD(e))*s/(np.sqrt(1+s)+1)
def mean_std(x,w):
    avg=np.sum(w*x);return float(avg),float(np.sqrt(np.sum(w*(x-avg)**2)))

def analyze(directory):
    d=Path(directory);case=json.loads((d/'case.json').read_text());params=case['model']
    model=Sphere(params['kind'],params['electrons'],params['sigma_m'],params['nr'],params['nmu'],params['nphi'])
    expected=np.load(d/'reference_input/expected_events.npz');M=case['macroparticles']
    outputs=[]
    for p in case['planes']:
        arr=read_closed(d/'detectors'/(p['name']+'.h5'),'particle_plane')
        assert len(arr)==M and len(np.unique(arr['source_id']))==M,(p,len(arr))
        arr=arr[np.argsort(arr['source_id'])];assert np.array_equal(arr['source_id'],np.arange(1,M+1))
        assert np.all(arr['crossing_direction']==1)
        outputs.append(arr)
    entry=outputs[0];weights=np.asarray(abs(entry['charge_C']),dtype=LD);weights/=np.sum(weights)
    assert np.max(abs(weights-np.asarray(model.weights,dtype=LD)))<1e-13
    offsets=np.asarray(expected['entry__lab_time_s'],dtype=LD)-np.asarray(entry['time_s'],dtype=LD)
    offset=np.median(offsets);spread=np.max(abs(offsets-offset))
    assert spread<LD('1e-18'),'Entry time mapping is not a constant offset'
    entry_u=rest_u(entry['proper_velocity']);entry_k=kinetic_ev(entry_u)
    ref_entry_u=np.asarray(expected['entry__rest_u'],dtype=LD);ref_entry_k=kinetic_ev(ref_entry_u)
    entry_lab_gamma=np.sqrt(1+np.sum(np.asarray(entry['proper_velocity'],dtype=LD)**2,axis=1))
    ref_entry_lab_gamma=np.sqrt(1+np.sum(np.asarray(expected['entry__lab_u'],dtype=LD)**2,axis=1))
    rest_gamma=LD(G_BEAM);rest_beta=np.sqrt((rest_gamma-1)*(rest_gamma+1))/rest_gamma
    rows=[];per_id={}
    for p,arr in zip(case['planes'],outputs):
        t=np.asarray(arr['time_s'],dtype=LD)+offset;z=np.asarray(arr['position_m'][:,2],dtype=LD)
        T=LD(model.entry_center_time)+rest_gamma*(t-rest_beta*z/LD(c))
        tau=np.asarray(T/LD(model.t0),dtype=float)
        assert tau.min()>=-1e-7 and tau.max()<6
        rad,p_hat=model.shell_state(model.shell,tau)
        ru_ref=np.asarray((np.sqrt(model.eta)*p_hat)[:,None]*model.directions,dtype=LD)
        rp_ref=np.asarray((rad*model.sigma)[:,None]*model.directions,dtype=LD)
        ru=rest_u(arr['proper_velocity']);pos=np.asarray(arr['position_m'],dtype=LD).copy()
        pos[:,2]=rest_gamma*(z-rest_beta*LD(c)*t)
        gain=kinetic_ev(ru)-entry_k;refgain=kinetic_ev(ru_ref)-ref_entry_k
        mean,std=mean_std(gain,weights);rmean,rstd=mean_std(refgain,weights)
        sigma=lambda x:float(np.sqrt(np.sum(weights*(x-np.sum(weights*x))**2)))
        poserr=float(np.sqrt(np.sum(weights*np.sum((pos-rp_ref)**2,axis=1))/np.sum(weights*np.sum(rp_ref**2,axis=1))))
        uerr=float(np.sqrt(np.sum(weights*np.sum((ru-ru_ref)**2,axis=1))/np.sum(weights*np.sum(ru_ref**2,axis=1))))
        actual_lab=np.sqrt(1+np.sum(np.asarray(arr['proper_velocity'],dtype=LD)**2,axis=1))
        expected_lab=rest_gamma*(np.sqrt(1+np.sum(ru_ref**2,axis=1))+rest_beta*ru_ref[:,2])
        rest_mc2_ev=LD(m_e)*LD(c)**2/LD(e)
        lm,ls=mean_std((actual_lab-entry_lab_gamma)*rest_mc2_ev,weights)
        rm,rs=mean_std((expected_lab-ref_entry_lab_gamma)*rest_mc2_ev,weights)
        row={'plane':p['name'],'lab_z_m':p['z'],'mean_rest_tau':float(np.sum(weights*tau)),
          'minimum_rest_tau':float(tau.min()),'maximum_rest_tau':float(tau.max()),
          'PIC_mean_rest_gain_eV':mean,'reference_mean_rest_gain_eV':rmean,
          'PIC_std_rest_gain_eV':std,'reference_std_rest_gain_eV':rstd,
          'PIC_mean_lab_gain_eV':lm,'reference_mean_lab_gain_eV':rm,
          'PIC_std_lab_gain_eV':ls,'reference_std_lab_gain_eV':rs,
          'PIC_sigma_x_m':sigma(pos[:,0]),'reference_sigma_x_m':sigma(rp_ref[:,0]),
          'PIC_sigma_z_rest_m':sigma(pos[:,2]),'reference_sigma_z_rest_m':sigma(rp_ref[:,2]),
          'position_relative_L2':poserr,'proper_velocity_relative_L2':uerr,
          'particle_energy_change_error_rms_eV':float(np.sqrt(np.sum(weights*(gain-refgain)**2)))}
        rows.append(row)
        per_id[p['name']+'__gain_rest_eV']=np.asarray(gain,dtype=float)
        per_id[p['name']+'__reference_gain_rest_eV']=np.asarray(refgain,dtype=float)
    # Each on-disk trajectory rank must be closed, finite and consistent.
    trajectory_rows=0;exit_faces={}
    paths=sorted((d/'trajectory').glob('*.h5'));assert len(paths)==16
    for p in paths:
        tr=read_closed(p,'trajectory');trajectory_rows+=len(tr)
        for face in range(6):exit_faces[str(face)]=exit_faces.get(str(face),0)+int(np.count_nonzero((tr['event_type']!=0)&(tr['boundary_face']==face)))
    assert sum(exit_faces[str(k)] for k in range(4))==0,'Particles entered transverse absorbing layers'
    ledger_rows,sums,ranks,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5')
    index,quantities,budget=report_discrete.calculate(ledger_rows,sums,steps)
    out=d/'analysis';out.mkdir(exist_ok=False)
    with (out/'comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    np.savez(out/'per_particle_gain.npz',**per_id,weights=np.asarray(weights,dtype=float))
    absolute,pct,normalization=construct(ledger_rows,sums,index,quantities)
    write_tables(out,index,absolute,pct,normalization)
    result={'execution_integrity_pass':True,'physics_accepted':False,
      'status':'EXPLORATORY_PHYSICS_COMPARISON_MEASURED_NOT_CONVERGENCE_ACCEPTANCE',
      'entry_time_offset_s':float(offset),'entry_time_offset_spread_s':float(spread),
      'entry_max_xy_error_m':float(np.max(abs(entry['position_m'][:,:2]-expected['entry__xy_m']))),
      'entry_max_rest_u_error':float(np.max(abs(entry_u-ref_entry_u))),
      'entry_rest_mean_kinetic_error_eV':float(np.sum(weights*(entry_k-ref_entry_k))),
      'trajectory_records':trajectory_rows,'trajectory_exit_faces':exit_faces,
      'steps':steps,'particles_per_plane':M,'planes':rows,
      'new_field_operator_maximum_residual_J':budget['max_abs_cumulative_field_operator_residual_J'],
      'new_total_final_residual_J':budget['final']['discrete_total_residual_J'],
      'particle_grid_exchange_final_J':budget['final']['particle_grid_coupling_gap_J'],
      'relative_energy_budget':normalization,
      'time_scope':'Each PIC crossing compared to same source ID reference worldline at its reconstructed rest-frame event time. Energy gain is relative to each ID entry event, not a simultaneous bunch energy.',
      'reference_model':params}
    save(out/'RESULT.json',result)
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained');x=[r['mean_rest_tau'] for r in rows]
    for ax,actual,ref,title,unit in [
      (axes[0,0],'PIC_mean_rest_gain_eV','reference_mean_rest_gain_eV','Mean energy change from entry','eV / electron'),
      (axes[0,1],'PIC_std_rest_gain_eV','reference_std_rest_gain_eV','Std of individual energy changes','eV / electron'),
      (axes[1,0],'PIC_sigma_x_m','reference_sigma_x_m','Transverse rms size','um'),
      (axes[1,1],'PIC_mean_lab_gain_eV','reference_mean_lab_gain_eV','Lab mean energy change from entry','eV / electron')]:
        scale=1e6 if unit=='um' else 1
        ax.plot(x,[r[actual]*scale for r in rows],'o-',color='#4b78a4',label='APRL')
        ax.plot(x,[r[ref]*scale for r in rows],'x--',color='#d38b39',label='Spherical reference')
        ax.set_title(title);ax.set_ylabel(unit);ax.set_xlabel('Mean rest-frame event time / t0');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(params['kind']+': first PIC comparison; fixed lab-plane events, no physical acceptance assumed',fontsize=11)
    fig.savefig(out/'comparison.png',dpi=180);fig.savefig(out/'comparison.pdf');plt.close(fig)
    plot_relative(out,index,pct,normalization,params['kind']+': relative discrete energy budget')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();print(json.dumps(analyze(a.directory),indent=2))
