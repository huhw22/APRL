#!/usr/bin/env python3
"""Common boost-time cohort and measured pusher-field diagnosis, no fitted drift."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-longitudinal-origin-mpl')
import argparse,csv,json,sys
from pathlib import Path
import h5py,numpy as np
from scipy.constants import c,e,m_e,epsilon_0
try:
    from model import Sphere,G_BEAM,G_BOOST,DT
except ImportError:
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'space_charge_physics'))
    from model import Sphere,G_BEAM,G_BOOST,DT

LD=np.longdouble
rap=np.arccosh(LD(G_BEAM))-np.arccosh(LD(G_BOOST))
GR=np.cosh(rap);BR=np.tanh(rap)
REC=np.dtype([('step','<u8'),('source','<u8'),('time','<f8'),('dt','<f8'),('charge','<f8'),
    ('position','<f8',(3,)),('before','<f8',(3,)),('after','<f8',(3,)),('electric','<f8',(3,)),('magnetic','<f8',(3,))])
assert REC.itemsize==160

def write_csv(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')

def boris(u,E,B,dt):
    u=np.asarray(u,dtype=LD);E=np.asarray(E,dtype=LD);B=np.asarray(B,dtype=LD)
    kick=(-LD(e)/LD(m_e)/LD(c))*np.asarray(dt,dtype=LD)/2
    minus=u+kick[:,None]*E
    gamma=np.sqrt(1+np.sum(minus**2,axis=1))
    t=(-LD(e)/LD(m_e))*np.asarray(dt,dtype=LD)[:,None]*B/(2*gamma[:,None])
    prime=minus+np.cross(minus,t);s=2*t/(1+np.sum(t*t,axis=1))[:,None]
    plus=minus+np.cross(prime,s)
    return plus+kick[:,None]*E

def analyze(directory, original, out, kind='uniform'):
    d=Path(directory);original=Path(original);out=Path(out);out.mkdir(exist_ok=True)
    config=json.loads((original/'case.json').read_text())
    result=json.loads((original/'analysis/RESULT.json').read_text())
    model=Sphere(kind)
    with h5py.File(d/'ledger/energy-ledger.h5') as h:
        g=h['energy_ledger'];assert int(g['complete'][()])==1
        rows=g['records'][:];assert len(rows)==int(g['committed_records'][()])
    release=int(np.flatnonzero(abs(rows['released_charge_C'])>=model.electrons*e*(1-1e-12))[0])
    exits=np.flatnonzero(rows['removed_macroparticles']>0);first_exit=int(exits[0]) if len(exits) else len(rows)
    shell=np.repeat(np.arange(model.nr),model.nmu);mu=np.tile(model.mu,model.nr)
    w=model.weights.reshape(model.nr,model.nmu,model.nphi).sum(axis=2).ravel()
    totalrest=LD(model.electrons)*LD(m_e)*LD(c)**2
    origin=LD(config['geometry']['initial_frame_lab_time_s']);offset=LD(result['entry_time_offset_s'])
    comparisons=[]
    for i in sorted(set([release,first_exit-1]+list(range(0,len(rows),20)))):
        if i<release or i>=first_exit:continue
        r=rows[i]
        rhs=(model.entry_center_time+(LD(r['time_box_s'])+G_BOOST*(origin+offset))/GR)/model.t0
        tau=np.full(len(shell),float(rhs))
        for _ in range(8):
            rad,p=model.shell_state(shell,tau);tau=np.asarray(rhs-BR*rad*model.sigma*mu/(c*model.t0),dtype=float)
        rad,p=model.shell_state(shell,tau);uz=np.sqrt(model.eta)*p*mu;gamma=np.sqrt(1+model.eta*p*p)
        expectedpz=totalrest*GR*np.sum(w*(uz+BR*gamma))
        expectedk=totalrest*np.sum(w*(GR*(gamma+BR*uz)-1))
        dp=(LD(r['particle_cPz_box_active_J'])-expectedpz)/totalrest
        dk=(LD(r['particle_kinetic_active_J'])-expectedk)/totalrest
        comparisons.append(dict(step=i,mean_reference_rest_tau=float(np.sum(w*tau)),
          mean_uz_box_error=float(dp),mean_rest_uz_error=float(GR*(dp-BR*dk)),
          mean_box_energy_error_eV=float(dk*LD(m_e)*LD(c)**2/LD(e)),
          lab_energy_error_from_same_box_slice_eV=float(LD(G_BOOST)*(dk+np.sqrt(1-LD(G_BOOST)**-2)*dp)*LD(m_e)*LD(c)**2/LD(e))))
    write_csv(out/'common_boost_time.csv',comparisons)
    force_rows=[];replay_error=0.;coverage=[]
    traces=sorted((d/'momentum').glob('*.bin')) if (d/'momentum').exists() else []
    if traces:
        data=[]
        for f in traces:
            assert f.stat().st_size%REC.itemsize==0
            a=np.fromfile(f,dtype=REC)
            for n in REC.names:
                if a[n].dtype.kind=='f':assert np.isfinite(a[n]).all()
            data.append(a)
        a=np.concatenate(data)
        for step in np.unique(a['step']):
            if step<release or step>=first_exit:continue
            sample=a[a['step']==step];sample=sample[np.argsort(sample['source'])]
            assert np.array_equal(sample['source'],np.arange(1,len(model.weights)+1)),(step,len(sample))
            assert np.max(abs(sample['dt']-DT))<DT*1e-9
            coverage.append(int(step))
            u=sample['before'].astype(LD);after=sample['after'].astype(LD)
            replay=boris(u,sample['electric'],sample['magnetic'],sample['dt'])
            replay_error=max(replay_error,float(np.max(abs(replay-after))))
            assert replay_error<1e-13
            # Position in the numerical box includes its translated spatial origin.
            zphysical=sample['position'][:,2].astype(LD)-LD(config['geometry']['box_translation_m'])
            tau=(model.entry_center_time+LD(G_BEAM)*(origin+offset)
                 +GR*(sample['time'].astype(LD)-BR*zphysical/LD(c)))/model.t0
            gs=LD(G_BOOST);bs=np.sqrt(1-gs**-2);gb=LD(G_BEAM);bb=np.sqrt(1-gb**-2)
            labtime=origin+gs*(sample['time'].astype(LD)+bs*zphysical/LD(c))
            labz=gs*(zphysical+bs*LD(c)*sample['time'])
            tau_direct=(model.entry_center_time+gb*(labtime+offset-bb*labz/LD(c)))/model.t0
            assert np.max(abs(tau-tau_direct))<1e-8
            rad,p=model.shell_state(model.shell,np.asarray(tau,dtype=float))
            estar=(-model.electrons*e/(4*np.pi*epsilon_0)*model.q[model.shell]/(rad*model.sigma)**2)[:,None]*model.directions
            # Evaluate the same continuum density at the actual PIC position too,
            # so centroid/shape displacement is not mislabelled as a field error.
            restpos=sample['position'].astype(LD).copy()
            restpos[:,2]=GR*(zphysical-BR*LD(c)*sample['time'])-gb*bb*LD(c)*(origin+offset)
            rr=np.sqrt(np.sum(restpos**2,axis=1))
            radial=model.solution.sol(np.asarray(tau,dtype=float))[:model.nr].T*model.sigma
            ii=np.sum(radial<=np.asarray(rr)[:,None],axis=1)
            lo=np.clip(ii-1,0,model.nr-2);hi=lo+1;index=np.arange(len(rr))
            qenc=model.q[lo]+(model.q[hi]-model.q[lo])*(rr-radial[index,lo])/(radial[index,hi]-radial[index,lo])
            inner=ii==0;qenc[inner]=model.q[0]*(rr[inner]/radial[inner,0])**3
            qenc=np.clip(qenc,0,1)
            pointE=(-LD(model.electrons)*LD(e)/(4*LD(np.pi)*LD(epsilon_0))*qenc/rr**3)[:,None]*restpos
            refE=pointE.copy();refE[:,:2]*=GR
            refB=np.c_[-BR*refE[:,1]/c,BR*refE[:,0]/c,np.zeros(len(refE))]
            ideal=boris(u,refE,refB,sample['dt'])
            factor=-LD(e)*sample['dt'].astype(LD)/(LD(m_e)*LD(c))
            electric=factor*sample['electric'][:,2];ideal_e=factor*refE[:,2]
            magnetic=after[:,2]-u[:,2]-electric;ideal_b=ideal[:,2]-u[:,2]-ideal_e
            weights=np.abs(sample['charge']).astype(LD);weights/=weights.sum()
            avg=lambda v:float(np.sum(weights*v))
            force_rows.append(dict(step=int(step),mean_duz_electric=avg(electric),mean_duz_magnetic=avg(magnetic),
              reference_field_duz_electric=avg(ideal_e),reference_field_duz_magnetic=avg(ideal_b),
              electric_impulse_error=avg(electric-ideal_e),magnetic_impulse_error=avg(magnetic-ideal_b),
              net_field_impulse_error=avg(after[:,2]-ideal[:,2]),
              mean_Ez_error_V_per_m=avg(sample['electric'][:,2]-refE[:,2]),
              mean_Ez_error_worldline_position_V_per_m=avg(sample['electric'][:,2]-estar[:,2]),
              rms_Ez_error_V_per_m=float(np.sqrt(np.sum(weights*(sample['electric'][:,2]-refE[:,2])**2)))))
        write_csv(out/'measured_force.csv',force_rows)
    summary=dict(all_released_step=release,first_exit_step=first_exit if len(exits) else None,
      steps=int(rows['step'][-1]),frame='Common program boost-time hypersurfaces; full released cohort before any exit',
      final_clean_comparison=comparisons[-1],probe_steps=coverage,probe_records=len(a) if traces else 0,
      independent_boris_replay_max_absolute_u_error=replay_error,force=force_rows,
      field_reference_scope='Spherical reference density at actual PIC positions and actual reconstructed times, with translated box origin included; not a Poisson solve of the perturbed PIC density. Worldline-position comparison also retained.')
    save(out/'ANALYSIS.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='force'},indent=2))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('original',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--kind',default='uniform');a=p.parse_args();analyze(a.directory,a.original,a.output,a.kind)
