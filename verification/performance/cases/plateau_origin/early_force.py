#!/usr/bin/env python3
"""Actual pusher samples before every source has entered, with analytic force.

This is a subset diagnosis: mean force is over only the already active source
IDs and must not be compared directly to the full-cohort mean-energy curve.
"""
import csv,json,sys
from pathlib import Path
import numpy as np

def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text());parent=Path(plan['parent']);workspace=Path(parent.parents[3])
    b=workspace/'source/aprl/verification/performance/cases'
    sys.path.insert(0,str(root/'scripts'));sys.path.insert(1,str(b/'surface_residual_origin'))
    from analyze_momentum import REC,LD,GR,BR,G_BEAM,G_BOOST,DT,c,e,m_e,boris,Sphere
    from force_reference import DenseSphereForce
    model=Sphere('converging');dense=DenseSphereForce(model)
    original=json.loads((Path(plan['cases']['guard128_500']['original'])/'case.json').read_text())
    prior=json.loads((Path(plan['cases']['guard128_500']['original'])/'analysis/RESULT.json').read_text())
    origin=LD(original['geometry']['initial_frame_lab_time_s']);offset=LD(prior['entry_time_offset_s']);gb=LD(original['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    dirs={32:parent.parent/'charge_matched_inlet_aux_20261030/guard32_500',64:parent/'pilot_500',128:root/'guard128_500',256:root/'guard256_500'}
    rows=[]
    for guard,d in dirs.items():
        if not (d/'momentum').exists() or not list((d/'momentum').glob('*.bin')):continue
        a=np.concatenate([np.fromfile(p,dtype=REC) for p in sorted((d/'momentum').glob('*.bin'))]);assert len(a)>0
        for step in (20,40,60,80,100,120,140,160,180,200,240,280,320,400,480):
            sample=a[a['step']==step]
            if len(sample)==0:continue
            sample=sample[np.argsort(sample['source'])];assert np.array_equal(sample['source'],np.arange(1,len(sample)+1)) if step>=140 else len(np.unique(sample['source']))==len(sample)
            time=sample['time'].astype(LD);pos=sample['position'].astype(LD).copy();pos[:,2]-=LD(original['geometry']['box_translation_m'])
            tau=(model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/model.t0
            pos[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
            E=dense.field(pos,tau).astype(LD);E[:,:2]*=GR
            B=np.c_[-BR*E[:,1]/c,BR*E[:,0]/c,np.zeros(len(E))]
            u=sample['before'].astype(LD);after=sample['after'].astype(LD)
            weights=np.abs(sample['charge']).astype(LD);total=float(np.sum(weights));weights/=weights.sum()
            avg=lambda v:float(np.sum(weights*v))
            ef=-LD(e)*sample['dt'].astype(LD)/(LD(m_e)*LD(c))*(sample['electric'][:,2]-E[:,2]);ideal=boris(u,E,B,sample['dt'])
            replay=boris(u,sample['electric'],sample['magnetic'],sample['dt'])
            assert np.max(abs(replay-after))<1e-13
            rows.append(dict(guard=guard,step=step,source_count=len(sample),source_charge_abs_C=total,
                mean_Ez_error_V_per_m=avg(sample['electric'][:,2]-E[:,2]),
                electric_duz_error=avg(ef),net_duz_error=avg(after[:,2]-ideal[:,2]),
                magnetic_duz_error=avg(after[:,2]-ideal[:,2]-ef),
                rms_Ez_error_V_per_m=float(np.sqrt(np.sum(weights*(sample['electric'][:,2]-E[:,2])**2)))))
    out=root/'EARLY_FORCE.csv'
    with out.open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    print(json.dumps([x for x in rows if x['guard'] in (32,64) and x['step'] in (20,40,60,80,100,120,140,180,280,480)],indent=2))
if __name__=='__main__':main(sys.argv[1])
