#!/usr/bin/env python3
"""Compare the first force on every admitted particle with the ideal sphere."""
import csv,hashlib,json,sys
from pathlib import Path
import h5py,numpy as np,yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text());d=Path(plan['case']);prior=Path(plan['comparison'])
    assert json.loads((d/'COMPLETION.json').read_text())['checks_pass']
    for p,h in json.loads((root/'FROZEN.json').read_text()).items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
    bench=root.parent;parent=bench/'charge_matched_inlet_20261030'
    sys.path[:0]=[str(parent/'scripts'),str(root.parents[3]/'source/aprl/verification/performance/cases/surface_residual_origin')]
    from analyze_momentum import REC,LD,GR,BR,G_BEAM,c,boris,Sphere
    from force_reference import DenseSphereForce
    original=Path(json.loads((parent/'plan.json').read_text())['cases']['pilot_500']['original'])
    cfg=json.loads((original/'case.json').read_text());prior_analysis=json.loads((original/'analysis/RESULT.json').read_text())
    origin=LD(cfg['geometry']['initial_frame_lab_time_s']);offset=LD(prior_analysis['entry_time_offset_s']);gb=LD(cfg['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    model=Sphere('converging');dense=DenseSphereForce(model)
    parts=[]
    for p in sorted((d/'momentum').glob('*.bin')):
        assert p.stat().st_size%REC.itemsize==0
        parts.append(np.fromfile(p,dtype=REC))
    a=np.concatenate(parts);raw_records=len(a)
    a=a[np.lexsort((a['time'],a['step'],a['source']))]
    a=a[np.r_[True,a['source'][1:]!=a['source'][:-1]]]
    assert len(a)==32768 and np.array_equal(a['source'],np.arange(1,32769))
    assert np.max(a['step'])<=123
    window=np.genfromtxt(d/'reference/window-rank-0.csv',names=True,delimiter=',')
    cell_z=float(yaml.safe_load((d/'run.yaml').read_text())['mesh']['cell_size'][2])
    assert int(window['step'][0])==1
    initial_left=float(window['origin_box_z_m'][0])+cell_z
    gap_cells=(a['position'][:,2]-(initial_left-a['step']*cell_z))/cell_z
    causal_margin=gap_cells-a['step']
    assert np.min(causal_margin)>0
    time=a['time'].astype(LD);pos=a['position'].astype(LD).copy();pos[:,2]-=LD(cfg['geometry']['box_translation_m'])
    tau=(model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/model.t0
    pos[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
    E=dense.field(pos,tau).astype(LD);E[:,:2]*=GR
    B=np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
    ideal=boris(a['before'],E,B,a['dt']);replay=boris(a['before'],a['electric'],a['magnetic'],a['dt'])
    max_replay=float(np.max(abs(replay-a['after'])));assert max_replay<1e-13
    w=np.abs(a['charge']).astype(LD);avg=lambda v:float(np.sum(w*v,dtype=LD)/np.sum(w,dtype=LD))
    error=a['electric'][:,2]-E[:,2];impulse=a['after'][:,2]-ideal[:,2]
    rows=[]
    for lo,hi in ((0,20),(20,40),(40,60),(60,80),(80,100),(100,124),(0,124)):
        m=(a['step']>=lo)&(a['step']<hi)
        wt=w[m];mean=lambda v:float(np.sum(wt*v[m],dtype=LD)/np.sum(wt,dtype=LD))
        rows.append(dict(first_step_from=lo,first_step_before=hi,count=int(np.sum(m)),charge_C=float(np.sum(wt)),
            mean_Ez_error_V_per_m=mean(error),rms_Ez_error_V_per_m=float(np.sqrt(np.sum(wt*error[m]**2)/np.sum(wt))),
            mean_direct_duz_error=mean(impulse)))
    with h5py.File(d/'ledger/energy-ledger.h5') as h,h5py.File(prior/'ledger/energy-ledger.h5') as g:
        x=h['energy_ledger/records'][:];y=g['energy_ledger/records'][:len(x)]
        assert x.shape==y.shape and x.dtype==y.dtype
        science={name:np.ascontiguousarray(x[name]).tobytes()==np.ascontiguousarray(y[name]).tobytes()
                 for name in x.dtype.names if 'seconds' not in name}
        assert all(science.values()),[name for name,v in science.items() if not v]
    out=root/'review';out.mkdir(exist_ok=True)
    with (out/'first_force.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    fig,ax=plt.subplots(figsize=(8,4.3),layout='constrained')
    x=[r['first_step_from'] for r in rows[:-1]];y=[r['mean_Ez_error_V_per_m']/1000 for r in rows[:-1]]
    ax.bar(x,y,width=19,align='edge',color='#59559c');ax.axhline(0,color='k',lw=.7)
    ax.set(xlabel='First PIC push step',ylabel='Mean sampled Ez error (kV/m)',title='Already present at first force sample (guard 256)')
    ax.grid(axis='y',alpha=.2)
    fig.savefig(out/'first_force.png',dpi=180);fig.savefig(out/'first_force.pdf');plt.close(fig)
    result=dict(job_id=json.loads((d/'COMPLETION.json').read_text())['job_id'],binary_sha256=plan['binary_sha256'],
                raw_records=raw_records,unique_earliest_records=len(a),mean_first_Ez_error_V_per_m=avg(error),mean_first_duz_error=avg(impulse),
                minimum_distance_to_dynamic_left_boundary_cells=float(np.min(gap_cells)),
                minimum_remaining_causal_cells_since_initial_time=float(np.min(causal_margin)),
                max_boris_replay_error=max_replay,science_ledger_prefix_bitwise=science,by_first_step=rows,
                interpretation='Ideal smooth sphere at actual PIC admission positions/times. Dynamic left-layer updates since t=0 cannot reach the first pushes at c*dt=dz. The initial global reference solve, discretization and finite particle representation remain possible; first-push direct impulse has the opposite sign to the later negative mean bias.')
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='science_ledger_prefix_bitwise'},indent=2))
if __name__=='__main__':main(sys.argv[1])
