#!/usr/bin/env python3
"""Close-read every pusher step and split late mean momentum drift.

The direct term compares actual Boris with Boris fed the independent smooth
sphere field at the *same PIC state*. The remainder includes trajectory and
density feedback, plus reference/discretization differences; it is not assigned
to one numerical component.
"""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
import h5py
from scipy.constants import m_e,e
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text());d=Path(plan['case']);prior=Path(plan['comparison'])
    assert json.loads((d/'COMPLETION.json').read_text())['checks_pass']
    for p,h in json.loads((root/'FROZEN.json').read_text()).items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
    bench=root.parent
    parent=bench/'charge_matched_inlet_20261030'
    case=json.loads((parent/'plan.json').read_text())['cases']['pilot_500'];original=Path(case['original'])
    cfg=json.loads((original/'case.json').read_text());old=json.loads((original/'analysis/RESULT.json').read_text())
    sys.path[:0]=[str(parent/'scripts'),str(root.parents[3]/'source/aprl/verification/performance/cases/surface_residual_origin')]
    from analyze_momentum import REC,LD,GR,BR,G_BEAM,G_BOOST,DT,c,boris,Sphere
    from force_reference import DenseSphereForce
    model=Sphere('converging');dense=DenseSphereForce(model)
    origin=LD(cfg['geometry']['initial_frame_lab_time_s']);offset=LD(old['entry_time_offset_s']);gb=LD(cfg['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    keys=('charge','actual_duz','ideal_duz','direct_duz','actual_dgamma','ideal_dgamma','direct_dgamma',
          'actual_ez','ez_error','b_direct_duz','before_uz','after_uz')
    accum={k:np.zeros(1001,dtype=np.longdouble) for k in keys};counts=np.zeros(1001,dtype=np.int64)
    files=sorted((d/'momentum').glob('*.bin'));assert len(files)==128
    max_replay=0.
    for path in files:
        assert path.stat().st_size%REC.itemsize==0
        with path.open('rb') as handle:
            while True:
                a=np.fromfile(handle,dtype=REC,count=40000)
                if len(a)==0:break
                assert np.all(a['step']%plan['probe_time_divisor']==0)
                step=(a['step']//plan['probe_time_divisor']).astype(np.int64)
                assert np.all((step>=0)&(step<1000))
                time=a['time'].astype(LD);pos=a['position'].astype(LD).copy();pos[:,2]-=LD(cfg['geometry']['box_translation_m'])
                tau=(model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/model.t0
                pos[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
                E=dense.field(pos,tau).astype(LD);E[:,:2]*=GR
                B=np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
                before=a['before'].astype(LD);after=a['after'].astype(LD)
                ideal=boris(before,E,B,a['dt'])
                actual=boris(before,a['electric'],a['magnetic'],a['dt'])
                max_replay=max(max_replay,float(np.max(np.abs(actual-after))))
                w=np.abs(a['charge']).astype(LD)
                gbefore=np.sqrt(1+np.sum(before*before,axis=1))
                gafter=np.sqrt(1+np.sum(after*after,axis=1))
                gideal=np.sqrt(1+np.sum(ideal*ideal,axis=1))
                values={'charge':w,'actual_duz':w*(after[:,2]-before[:,2]),'ideal_duz':w*(ideal[:,2]-before[:,2]),
                    'direct_duz':w*(after[:,2]-ideal[:,2]),
                    'actual_dgamma':w*(gafter-gbefore),'ideal_dgamma':w*(gideal-gbefore),
                    'direct_dgamma':w*(gafter-gideal),'actual_ez':w*a['electric'][:,2],
                    'ez_error':w*(a['electric'][:,2]-E[:,2]),
                    'b_direct_duz':w*((after[:,2]-ideal[:,2])+(LD(e)/LD(m_e)/c)*a['dt']*(a['electric'][:,2]-E[:,2])),
                    'before_uz':w*before[:,2],'after_uz':w*after[:,2]}
                for k,v in values.items():np.add.at(accum[k],step,v)
                np.add.at(counts,step,1)
    assert max_replay<1e-13,max_replay
    assert np.all(counts[123:1000]==32768)
    total=np.longdouble(model.electrons)*np.longdouble(1.602176634e-19)
    rows=[]
    for step in range(1000):
        w=accum['charge'][step]
        if w==0:continue
        row={'step':step,'source_count':int(counts[step]),'charge_C':float(w)}
        for k in keys[1:]:row['mean_'+k]=float(accum[k][step]/w)
        rows.append(row)
        if step>=123:assert abs(w-total)<total*1e-10,(step,w,total)
    out=root/'review';out.mkdir(exist_ok=True)
    with (out/'step_force.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    # The pusher trace is read only; verify both science records and sampled
    # force against the independent sparse run before attributing a drift.
    with h5py.File(d/'ledger/energy-ledger.h5') as h,h5py.File(prior/'ledger/energy-ledger.h5') as g:
        a=h['energy_ledger/records'][:];b=g['energy_ledger/records'][:]
        assert a.shape==b.shape and a.dtype==b.dtype
        science={name:np.ascontiguousarray(a[name]).tobytes()==np.ascontiguousarray(b[name]).tobytes()
                 for name in a.dtype.names if 'seconds' not in name}
        assert all(science.values()),[name for name,v in science.items() if not v]
    sparse=np.genfromtxt(prior/'momentum-analysis/measured_force.csv',names=True,delimiter=',')
    lookup={r['step']:r for r in rows}
    assert max(abs(lookup[int(s['step'])]['mean_ez_error']-s['mean_Ez_error_V_per_m']) for s in sparse)<100.
    # The previous reader used a coarser interpolation of the ideal sphere;
    # compare the raw sampled field itself, which is independent of that choice.
    sparse_actual={}
    for path in sorted((prior/'momentum').glob('*.bin')):
        a=np.fromfile(path,dtype=REC)
        for step in np.unique(a['step']):
            q=a[a['step']==step];v=sparse_actual.setdefault(int(step),[np.longdouble(0),np.longdouble(0)])
            w=np.abs(q['charge']).astype(LD)
            v[0]+=np.sum(w*q['electric'][:,2],dtype=LD);v[1]+=np.sum(w,dtype=LD)
    for step,(ez,w) in sparse_actual.items():
        assert abs(lookup[step]['mean_actual_ez']-float(ez/w))<1e-6,step
    moment=np.genfromtxt(prior/'momentum-analysis/common_boost_time.csv',names=True,delimiter=',')
    bystep={int(v['step']):v for v in moment}
    def interval(lo,hi):
        direct=float(np.sum(accum['direct_duz'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        actual=float(np.sum(accum['actual_duz'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        ideal=float(np.sum(accum['ideal_duz'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        observed=float(bystep[hi]['mean_uz_box_error']-bystep[lo]['mean_uz_box_error'])
        remainder=observed-direct
        direct_g=float(np.sum(accum['direct_dgamma'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        actual_g=float(np.sum(accum['actual_dgamma'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        ideal_g=float(np.sum(accum['ideal_dgamma'][lo:hi]/accum['charge'][lo:hi],dtype=LD))
        mc2=float(m_e*c*c/e)
        observed_g=float((bystep[hi]['mean_box_energy_error_eV']-bystep[lo]['mean_box_energy_error_eV'])/mc2)
        remainder_g=observed_g-direct_g
        return dict(from_step=lo,to_step=hi,observed_mean_uz_error_change=observed,direct_sampled_field_contribution=direct,
                    orbit_density_reference_remainder=remainder,actual_mean_uz_change=actual,ideal_on_PIC_state_mean_uz_change=ideal,
                    observed_mean_gamma_error_change=observed_g,direct_mean_gamma_contribution=direct_g,
                    gamma_orbit_density_reference_remainder=remainder_g,actual_mean_gamma_change=actual_g,
                    ideal_on_PIC_state_mean_gamma_change=ideal_g)
    intervals=[interval(lo,hi) for lo,hi in ((123,180),(180,300),(300,500),(500,700),(700,900),(900,1000),(500,1000),(123,1000))]
    mc2=float(m_e*c*c/e);boostbeta=float(np.sqrt(1-G_BOOST**-2))
    labfactor=float(G_BOOST*boostbeta*mc2)
    for v in intervals:
        for k in ('observed_mean_uz_error_change','direct_sampled_field_contribution','orbit_density_reference_remainder'):
            v[k+'_lab_momentum_eV']=v[k]*labfactor
        v['observed_lab_energy_error_change_eV']=G_BOOST*mc2*(v['observed_mean_gamma_error_change']+boostbeta*v['observed_mean_uz_error_change'])
        v['direct_lab_energy_contribution_eV']=G_BOOST*mc2*(v['direct_mean_gamma_contribution']+boostbeta*v['direct_sampled_field_contribution'])
        v['orbit_density_reference_lab_energy_remainder_eV']=G_BOOST*mc2*(v['gamma_orbit_density_reference_remainder']+boostbeta*v['orbit_density_reference_remainder'])
        expected=float(bystep[v['to_step']]['lab_energy_error_from_same_box_slice_eV']-bystep[v['from_step']]['lab_energy_error_from_same_box_slice_eV'])
        assert abs(v['observed_lab_energy_error_change_eV']-expected)<1e-6,(v['from_step'],v['to_step'])
    assert max(abs(v['observed_mean_uz_error_change']-v['direct_sampled_field_contribution']-v['orbit_density_reference_remainder']) for v in intervals)<1e-18
    result=dict(binary_sha256=plan['binary_sha256'],job_id=json.loads((d/'COMPLETION.json').read_text())['job_id'],
                probe_files=len(files),trace_records=int(counts.sum()),max_boris_replay_error=max_replay,
                science_ledger_prefix_bitwise=science,intervals=intervals,
                interpretation='Direct term uses ideal smooth-sphere E/B at each actual PIC state. Remainder combines changed trajectories/density and any continuum-discrete mismatch; it is not one boundary term. Lab momentum values omit boost kinetic contribution.')
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    x=np.array([r['step'] for r in rows]);direct=np.array([r['mean_direct_duz'] for r in rows]);
    direct_g=np.array([r['mean_direct_dgamma'] for r in rows]);z=np.array([r['mean_ez_error'] for r in rows]);
    seen=moment['step']>=123;anchors=moment[seen]
    per_step_direct=np.zeros(1000);per_step_direct_g=np.zeros(1000)
    per_step_direct[x.astype(int)]=direct;per_step_direct_g[x.astype(int)]=direct_g
    per_step_lab=G_BOOST*mc2*(per_step_direct_g+boostbeta*per_step_direct)
    cum=np.r_[0,np.cumsum(per_step_lab[123:1000])]
    fig,ax=plt.subplots(2,1,figsize=(10,7),layout='constrained')
    anchor_energy=anchors['lab_energy_error_from_same_box_slice_eV']-float(anchors[0]['lab_energy_error_from_same_box_slice_eV'])
    ax[0].plot(anchors['step'],anchor_energy/1000,label='Observed mean energy error change')
    ax[0].plot(np.arange(123,1001),cum/1000,label='Direct field contribution on PIC states')
    ax[0].plot(anchors['step'],(anchor_energy-cum[anchors['step'].astype(int)-123])/1000,label='Orbit/density/reference remainder')
    ax[0].set(ylabel='Lab-equivalent error change (keV/electron)',title='Mean energy error after complete release')
    ax[1].plot(x,z/1000,lw=.7,label='Every step')
    ax[1].plot(sparse['step'],sparse['mean_Ez_error_V_per_m']/1000,'.',ms=3,label='Prior 20-step samples')
    ax[1].set(xlabel='Field step',ylabel='Mean sampled $E_z$ error (kV/m)')
    for a in ax:a.grid(alpha=.2);a.legend()
    fig.savefig(out/'dense_force.png',dpi=180);fig.savefig(out/'dense_force.pdf');plt.close(fig)
    print(json.dumps({'intervals':intervals,'max_boris_replay_error':max_replay,'trace_records':int(counts.sum())},indent=2))
if __name__=='__main__':main(sys.argv[1])
