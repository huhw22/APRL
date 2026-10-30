#!/usr/bin/env python3
"""Read-only transverse mesh contrast at fixed physical aperture and inlet guard."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np,yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(p):return np.genfromtxt(p,names=True,delimiter=',')
def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text());d=root/'xy240_300';base=Path(plan['cases']['xy240_300']['comparison'])
    assert json.loads((d/'COMPLETION.json').read_text())['checks_pass']
    for p,h in json.loads((root/'FROZEN.json').read_text()).items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
    a=yaml.safe_load((base/'run.yaml').read_text());b=yaml.safe_load((d/'run.yaml').read_text())
    physical=lambda c:[c['mesh']['cells'][i]*c['mesh']['cell_size'][i] for i in (0,1)]
    pml=lambda c:[c['boundary']['cells'][i]*c['mesh']['cell_size'][i] for i in (0,1)]
    assert np.allclose(physical(a),physical(b),rtol=0,atol=1e-15)
    assert np.allclose(pml(a),pml(b),rtol=0,atol=1e-15)
    assert a['mesh']['moving_window']['lower_guard_cells']==b['mesh']['moving_window']['lower_guard_cells']==128
    assert a['beam']['input']['file']==b['beam']['input']['file']
    ma=load(base/'momentum-analysis/common_boost_time.csv');mb=load(d/'momentum-analysis/common_boost_time.csv')
    common=np.intersect1d(ma['step'],mb['step']);assert np.array_equal(common,mb['step'])
    fa=load(base/'momentum-analysis/measured_force.csv');fb=load(d/'momentum-analysis/measured_force.csv')
    commonf=np.intersect1d(fa['step'],fb['step']);assert np.array_equal(commonf,fb['step'])
    rows=[]
    for s in common:
        x=ma[ma['step']==s][0];y=mb[mb['step']==s][0]
        rows.append(dict(step=int(s),coarse_error_eV=float(x['lab_energy_error_from_same_box_slice_eV']),
                         fine_error_eV=float(y['lab_energy_error_from_same_box_slice_eV']),
                         fine_minus_coarse_eV=float(y['lab_energy_error_from_same_box_slice_eV']-x['lab_energy_error_from_same_box_slice_eV'])))
    out=root/'review';out.mkdir(exist_ok=True)
    with (out/'transverse_resolution.csv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    # The normal momentum reader starts at complete release; inspect the
    # already-active subset at the first admission steps independently.
    bench=root.parent;parent=bench/'charge_matched_inlet_20261030'
    sys.path[:0]=[str(parent/'scripts'),str(root.parents[3]/'source/aprl/verification/performance/cases/surface_residual_origin')]
    from analyze_momentum import REC,LD,GR,BR,G_BEAM,c,Sphere
    from force_reference import DenseSphereForce
    original=Path(json.loads((parent/'plan.json').read_text())['cases']['pilot_500']['original'])
    cfg=json.loads((original/'case.json').read_text());prior=json.loads((original/'analysis/RESULT.json').read_text())
    origin=LD(cfg['geometry']['initial_frame_lab_time_s']);offset=LD(prior['entry_time_offset_s']);gb=LD(cfg['model']['beam_gamma_lab']);bb=np.sqrt(1-gb**-2)
    dense=DenseSphereForce(Sphere('converging'))
    early=[]
    for name,path in [('coarse',base),('fine',d)]:
        for step in (20,40,60,80,100,120):
            parts=[]
            for f in sorted((path/'momentum').glob('*.bin')):
                x=np.fromfile(f,dtype=REC);q=x[x['step']==step]
                if len(q):parts.append(q)
            sample=np.concatenate(parts);time=sample['time'].astype(LD);pos=sample['position'].astype(LD).copy()
            pos[:,2]-=LD(cfg['geometry']['box_translation_m'])
            tau=(dense.model.entry_center_time+gb*(origin+offset)+GR*(time-BR*pos[:,2]/c))/dense.model.t0
            pos[:,2]=GR*(pos[:,2]-BR*c*time)-gb*bb*c*(origin+offset)
            Ez=dense.field(pos,tau)[:,2];w=np.abs(sample['charge']).astype(LD)
            early.append(dict(mesh=name,step=step,active_sources=len(sample),active_charge_C=float(np.sum(w)),
                              mean_Ez_error_V_per_m=float(np.sum(w*(sample['electric'][:,2]-Ez))/np.sum(w))))
    with (out/'early_active_field.csv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=early[0]);w.writeheader();w.writerows(early)
    step300=next(x for x in rows if x['step']==300)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    ax[0].plot(ma['step'],ma['lab_energy_error_from_same_box_slice_eV']/1000,label='5 µm / 160²')
    ax[0].plot(mb['step'],mb['lab_energy_error_from_same_box_slice_eV']/1000,label='3.33 µm / 240²')
    ax[0].set(title='Same cohort and physical aperture',xlabel='Field step',ylabel='Lab-equivalent mean error (keV/electron)')
    ax[1].plot(fa['step'],fa['mean_Ez_error_V_per_m']/1000,'.-',label='5 µm / 160²')
    ax[1].plot(fb['step'],fb['mean_Ez_error_V_per_m']/1000,'.-',label='3.33 µm / 240²')
    ax[1].set(title='Actual sampled field - ideal sphere',xlabel='Field step',ylabel='Mean Ez error (kV/m)')
    for x in ax:x.grid(alpha=.25);x.legend(fontsize=8)
    fig.savefig(out/'transverse_resolution.png',dpi=180);fig.savefig(out/'transverse_resolution.pdf');plt.close(fig)
    result=dict(job_id=json.loads((d/'COMPLETION.json').read_text())['job_id'],binary_sha256=plan['binary_sha256'],
                same_input=True,same_physical_aperture_m=physical(a),same_pml_thickness_m=pml(a),
                coarse_dx_m=a['mesh']['cell_size'][0],fine_dx_m=b['mesh']['cell_size'][0],
                step300=step300,relative_improvement_abs_error=(abs(step300['coarse_error_eV'])-abs(step300['fine_error_eV']))/abs(step300['coarse_error_eV']),
                entry_diagnostic=json.loads((d/'ENTRY_DIAGNOSTIC.json').read_text()),early_active_field=early,
                limitation='One refinement level and a 300-step endpoint; does not establish convergence or assign the remaining bias to the initial field alone.')
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(sys.argv[1])
