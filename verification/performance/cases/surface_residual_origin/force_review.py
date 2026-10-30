#!/usr/bin/env python3
"""Closed probe review with dense spherical force and explicit B half-time."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-residual-mpl')
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'longitudinal_momentum_origin'))
from analyze import Sphere,REC,GR,BR,G_BOOST,DT,c,e,m_e,LD,boris
from force_reference import DenseSphereForce

def main(root):
    plan=json.loads((root/'plan.json').read_text());orig=Path(plan['original']);cfg=json.loads((orig/'case.json').read_text())
    prior=json.loads((orig/'analysis/RESULT.json').read_text());model=Sphere('converging');dense=DenseSphereForce(model);coarse=DenseSphereForce(model,256,513)
    origin=LD(cfg['geometry']['initial_frame_lab_time_s']);offset=LD(prior['entry_time_offset_s']);gb=LD(cfg['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    summary={}
    def reference(sample,solver,half=False):
        time=sample['time'].astype(LD)+(LD('.5')*DT if half else 0)
        pos=sample['position'].astype(LD).copy();pos[:,2]-=LD(cfg['geometry']['box_translation_m'])
        tau=(model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/model.t0
        pos[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
        E=solver.field(pos,tau).astype(LD);E[:,:2]*=GR
        return E,np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
    for name in plan['cases']:
        d=root/name
        if not (d/'COMPLETION.json').exists():continue
        assert json.loads((d/'COMPLETION.json').read_text())['checks_pass']
        an=json.loads((d/'analysis/ANALYSIS.json').read_text());a=np.concatenate([np.fromfile(p,dtype=REC) for p in sorted((d/'momentum').glob('*.bin'))])
        rows=[];convergence=0.;replay=0.
        for step in np.unique(a['step']):
            if step<an['all_released_step']:continue
            sample=a[a['step']==step];sample=sample[np.argsort(sample['source'])]
            assert np.array_equal(sample['source'],np.arange(1,len(model.weights)+1))
            w=abs(sample['charge']).astype(LD);w/=w.sum();avg=lambda x:float(np.sum(w*x))
            u=sample['before'].astype(LD);after=sample['after'].astype(LD);E,B=reference(sample,dense);_,Bh=reference(sample,dense,True)
            if int(step) in (140,480):
                Ec,_=reference(sample,coarse);convergence=max(convergence,abs(avg(Ec[:,2]-E[:,2])))
            ideal=boris(u,E,B,sample['dt']);stagger=boris(u,E,Bh,sample['dt'])
            replay=max(replay,float(np.max(abs(boris(u,sample['electric'],sample['magnetic'],sample['dt'])-after))))
            factor=-LD(e)*sample['dt'].astype(LD)/(LD(m_e)*LD(c));de=factor*(sample['electric'][:,2]-E[:,2])
            rows.append(dict(step=int(step),mean_Ez_error_V_per_m=avg(sample['electric'][:,2]-E[:,2]),
                electric_impulse_error=avg(de),magnetic_impulse_error_same_time=avg(after[:,2]-ideal[:,2]-de),
                magnetic_impulse_error_half_time=avg(after[:,2]-stagger[:,2]-de),
                net_impulse_error_same_time=avg(after[:,2]-ideal[:,2]),net_impulse_error_half_time=avg(after[:,2]-stagger[:,2]),
                dense_reference_B_time_difference=avg(stagger[:,2]-ideal[:,2]),
                rms_Ez_error_V_per_m=float(np.sqrt(np.sum(w*(sample['electric'][:,2]-E[:,2])**2)))))
        out=d/'force-review';out.mkdir(exist_ok=True)
        with (out/'force.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
        s=dict(last=rows[-1],coarse_dense_reference_mean_Ez_difference_V_per_m=convergence,replay_max_u_error=replay,
            interpretation='Dense independent spherical continuum at actual positions; not a Poisson solve of actual perturbed PIC density. B half-time and same-time both retained.')
        (out/'SUMMARY.json').write_text(json.dumps(s,indent=2)+'\n');summary[name]=s
    (root/'FORCE_REVIEW.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root.resolve())
