#!/usr/bin/env python3
"""Separate transverse resolution from outer-interface placement."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np,yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load(p):return np.genfromtxt(p,names=True,delimiter=',')
def main(root):
    root=Path(root).resolve();bench=root.parent;prior=bench/'charge_inlet_plateau_origin_20261030_r3/guard128_500'
    refined=bench/'charge_inlet_plateau_refine_xy_20261030/xy240_300';wide=root/'xywide240_300'
    assert json.loads((wide/'COMPLETION.json').read_text())['checks_pass']
    for p,h in json.loads((root/'FROZEN.json').read_text()).items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
    cases={'base':prior,'fine_dx':refined,'wide_aperture':wide};cards={k:yaml.safe_load((p/'run.yaml').read_text()) for k,p in cases.items()}
    def geometry(c):
        m=c['mesh'];b=c['boundary'];return dict(cells_xy=m['cells'][:2],dxdy_m=m['cell_size'][:2],
            full_width_xy_m=[m['cells'][i]*m['cell_size'][i] for i in (0,1)],
            cpml_xy_m=[b['cells'][i]*m['cell_size'][i] for i in (0,1)],guard=m['moving_window']['lower_guard_cells'])
    geo={k:geometry(v) for k,v in cards.items()}
    assert all(v['guard']==128 for v in geo.values())
    assert all(v['cpml_xy_m']==[5e-5,5e-5] for v in geo.values())
    assert cards['base']['beam']['input']['file']==cards['fine_dx']['beam']['input']['file']==cards['wide_aperture']['beam']['input']['file']
    assert geo['base']['dxdy_m']==geo['wide_aperture']['dxdy_m']==[5e-6,5e-6]
    assert geo['fine_dx']['full_width_xy_m']==geo['base']['full_width_xy_m']
    assert geo['wide_aperture']['cells_xy']==geo['fine_dx']['cells_xy']==[240,240]
    data={k:load(p/'momentum-analysis/common_boost_time.csv') for k,p in cases.items()}
    force={k:load(p/'momentum-analysis/measured_force.csv') for k,p in cases.items()}
    common=np.intersect1d(np.intersect1d(data['base']['step'],data['fine_dx']['step']),data['wide_aperture']['step'])
    assert 300 in common and len(common)>5
    rows=[]
    for s in common:
        row={'step':int(s)}
        for k in cases:row[k+'_lab_error_eV']=float(data[k][data[k]['step']==s]['lab_energy_error_from_same_box_slice_eV'][0])
        rows.append(row)
    out=root/'review';out.mkdir(exist_ok=True)
    with (out/'aperture_vs_resolution.csv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    labels={'base':'5 µm, ±350 µm','fine_dx':'3.33 µm, ±350 µm','wide_aperture':'5 µm, ±550 µm'}
    for k in cases:
        ax[0].plot(data[k]['step'],data[k]['lab_energy_error_from_same_box_slice_eV']/1000,label=labels[k])
        ax[1].plot(force[k]['step'],force[k]['mean_Ez_error_V_per_m']/1000,'.-',label=labels[k])
    ax[0].set(title='Controlled transverse contrasts',xlabel='Field step',ylabel='Mean lab-equivalent error (keV/electron)')
    ax[1].set(title='Particle-sampled field',xlabel='Field step',ylabel='Mean Ez error (kV/m)')
    for a in ax:a.grid(alpha=.25);a.legend(fontsize=8)
    fig.savefig(out/'aperture_vs_resolution.png',dpi=180);fig.savefig(out/'aperture_vs_resolution.pdf');plt.close(fig)
    final=next(r for r in rows if r['step']==300)
    result=dict(job_id=json.loads((wide/'COMPLETION.json').read_text())['job_id'],binary_sha256=json.loads((root/'plan.json').read_text())['binary_sha256'],
                geometry=geo,step300=final,wide_minus_base_eV=final['wide_aperture_lab_error_eV']-final['base_lab_error_eV'],
                fine_minus_base_eV=final['fine_dx_lab_error_eV']-final['base_lab_error_eV'],
                entry_diagnostic=json.loads((wide/'ENTRY_DIAGNOSTIC.json').read_text()),
                limitation='Changing aperture also changes the finite Dirichlet reference domain and surface-current interface. This contrast does not separate those two boundary ingredients or establish convergence.')
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(sys.argv[1])
