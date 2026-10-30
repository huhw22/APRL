#!/usr/bin/env python3
"""Compare post-release momentum bias and sampled force for inlet guards."""
import csv,json,sys,hashlib
from pathlib import Path
import numpy as np
import h5py
from scipy.constants import c,m_e,e
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def read(p):return json.loads(p.read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def load(p):return np.genfromtxt(p,names=True,delimiter=',')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(root,late_root):
    plan=read(root/'plan.json');parent=Path(plan['parent']);dirs={32:parent.parent/'charge_matched_inlet_aux_20261030/guard32_500',64:parent/'pilot_500',128:root/'guard128_500',256:root/'guard256_500'}
    for folder in (root,late_root):
        assert all(sha(p)==h for p,h in read(folder/'FROZEN.json').items())
    late=late_root/'guard64_1000'
    assert read(late/'COMPLETION.json')['checks_pass']
    moment={n:load(p/'momentum-analysis/common_boost_time.csv') for n,p in dirs.items()}
    force={n:load(p/'momentum-analysis/measured_force.csv') for n,p in dirs.items()}
    assert all(read(root/f'guard{n}_500/COMPLETION.json')['checks_pass'] for n in (128,256))
    assert all(np.array_equal(moment[64]['step'],v['step']) for v in moment.values())
    assert all(np.array_equal(force[64]['step'],v['step']) for v in force.values())
    rows=[]
    for step in (123,140,160,180,200,240,280,300,400,500):
        i=int(np.flatnonzero(moment[64]['step']==step)[0]);row=dict(step=step)
        for n in dirs:row[f'guard{n}_lab_error_eV']=float(moment[n]['lab_energy_error_from_same_box_slice_eV'][i])
        rows.append(row)
    with (root/'comparison.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    summary={}
    for n in dirs:
        a=moment[n];f=force[n];summary[n]=dict(error_180_eV=float(a['lab_energy_error_from_same_box_slice_eV'][a['step']==180][0]),error_300_eV=float(a['lab_energy_error_from_same_box_slice_eV'][a['step']==300][0]),error_500_eV=float(a['lab_energy_error_from_same_box_slice_eV'][-1]),max_abs_sampled_mean_Ez_error_V_per_m=float(np.max(abs(f['mean_Ez_error_V_per_m']))),
            source=str(dirs[n]),entry_diagnostic=read(dirs[n]/'ENTRY_DIAGNOSTIC.json') if n in (128,256) else None)
    a=load(late/'momentum-analysis/common_boost_time.csv');original=load(parent/'converging_full/momentum-analysis/common_boost_time.csv')
    common=np.intersect1d(a['step'],original['step']);assert 500 in common and 1000 in common
    for step in common:
        x=a[a['step']==step][0];y=original[original['step']==step][0]
        for name in ('mean_uz_box_error','mean_box_energy_error_eV','lab_energy_error_from_same_box_slice_eV'):
            assert x[name]==y[name],(step,name)
    with h5py.File(late/'ledger/energy-ledger.h5') as h,h5py.File(parent/'converging_full/ledger/energy-ledger.h5') as g:
        x=h['energy_ledger/records'][:];y=g['energy_ledger/records'][:len(x)]
        assert int(h['energy_ledger/complete'][()])==int(g['energy_ledger/complete'][()])==1
        assert x.shape==y.shape and x.dtype==y.dtype
        fields={n:np.ascontiguousarray(x[n]).tobytes()==np.ascontiguousarray(y[n]).tobytes() for n in x.dtype.names if 'seconds' not in n}
        assert all(fields.values()),[n for n,v in fields.items() if not v]
    plt.rcParams.update({'font.size':9})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors={32:'#bd7b2d',64:'#267d70',128:'#59559c',256:'#b34562'}
    for guard in dirs:
        m=moment[guard];f=force[guard]
        ax[0,0].plot(m['step'],m['lab_energy_error_from_same_box_slice_eV']/1000,label=f'{guard} cells',color=colors[guard])
        ax[1,0].plot(f['step'],f['mean_Ez_error_V_per_m']/1000,label=f'{guard} cells',color=colors[guard],marker='.',ms=3)
        if guard!=64:
            ax[0,1].plot(m['step'],(m['lab_energy_error_from_same_box_slice_eV']-moment[64]['lab_energy_error_from_same_box_slice_eV'])/1000,label=f'{guard} - 64 cells',color=colors[guard])
    ax[0,0].set(title='Same released cohort after entry',ylabel='Lab-equivalent mean energy error (keV / electron)')
    ax[0,1].set(title='Guard effect at the same boost time',ylabel='Error difference from 64 cells (keV / electron)')
    ax[1,0].set(title='Actual sampled field minus independent sphere',ylabel='Mean Ez error (kV/m)',xlabel='Field step')
    ax[1,0].axvline(123,color='k',ls=':',lw=1,label='All released')
    ax[1,1].plot(a['step'],a['lab_energy_error_from_same_box_slice_eV']/1000,label='Total lab-equivalent error',color='#267d70')
    ax[1,1].plot(original['step'],original['lab_energy_error_from_same_box_slice_eV']/1000,'k:',lw=1,label='Original full record')
    g=929.1305492858543;beta=np.sqrt(1-g**-2);mc=m_e*c*c/e
    ax[1,1].plot(a['step'],g*a['mean_box_energy_error_eV']/1000,label='Boost kinetic contribution',color='#737d92')
    ax[1,1].plot(a['step'],g*beta*mc*a['mean_uz_box_error']/1000,label='Boost longitudinal momentum contribution',color='#c68135')
    ax[1,1].set(title='64-cell error into the later propagation',ylabel='Lab-equivalent error (keV / electron)',xlabel='Field step')
    for axes in ax.flat:axes.axhline(0,color='#999',lw=.6);axes.grid(alpha=.2);axes.legend(fontsize=7)
    out=root/'figures';out.mkdir(exist_ok=True);fig.savefig(out/'plateau_origin.png',dpi=180);fig.savefig(out/'plateau_origin.pdf');plt.close(fig)
    tail=load(late/'momentum-analysis/measured_force.csv')
    summary_late={int(s):dict(lab_error_eV=float(a['lab_energy_error_from_same_box_slice_eV'][a['step']==s][0]),mean_uz_box_error=float(a['mean_uz_box_error'][a['step']==s][0])) for s in (500,700,900,1000)}
    force_late={int(s):dict(mean_Ez_error_V_per_m=float(tail['mean_Ez_error_V_per_m'][tail['step']==s][0]),net_field_impulse_error=float(tail['net_field_impulse_error'][tail['step']==s][0])) for s in (500,600,700,800,900,980) if s in tail['step']}
    result=dict(summary=summary,late=summary_late,late_sampled_force=force_late,ledger_prefix_science_bitwise=fields,
        interpretation='Short guard sensitivity changes left boundary distance and active Nz. A common error at all guards cannot be assigned uniquely to the left boundary. The 1000-step sparse force trace and original full ledger prefix locate ongoing drift; exact source requires further controlled ablation.')
    save(root/'RESULT.json',result)
    print(json.dumps({k:{n:v[n] for n in ('error_180_eV','error_300_eV','error_500_eV')} for k,v in summary.items()},indent=2))
if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
