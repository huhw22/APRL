"""Small static read: one central closed input cell, never a running output."""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c

from prepare import save


def main(root):
    root=Path(root).resolve();p=json.loads((root/'PLAN.json').read_text())
    cell=p['active_slices']//2
    with h5py.File(p['exact_input'],'r') as f:
        a=f['particles/records'][cell*1024:(cell+1)*1024]
    phase=-2*np.pi*c*a['arrival_time_offset_s'].reshape(64,16)/p['original_plan']['physics']['wavelength_m']
    w=a['macro_weight'].reshape(64,16)
    floor=float(1/w.sum(axis=1)[0]);rows=[]
    for m in [1,2,3,8,16,20,32]:
        b=np.mean(np.exp(1j*m*phase),axis=1)
        rows.append(dict(harmonic=m,mean_conditional_family_b_squared=float(np.mean(abs(b)**2)),
                         independent_electron_Poisson_reference_b_squared=floor))
    result=dict(complete_static_sample=True,input_path=p['exact_input'],
        input_sha256=p['exact_input_sha256'],original_optical_cell_zero_based=cell,
        input_records=1024,represented_electrons=float(w.sum()),families=64,clones_per_family=16,
        unique_transverse_states=int(len(np.unique(a['plane_position_m'],axis=0))),
        rows=rows,scope='One central input cell only. Conditional-family statistics, not the global beam BF or a radiated-energy fraction. No claim of the fraction of enhanced FEL power explained.')
    save(root/'STATIC_SAMPLE.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(8.4,4.6),layout='constrained')
    x=np.array([r['harmonic'] for r in rows]);v=np.array([r['mean_conditional_family_b_squared'] for r in rows])
    ax.plot(x,v,'o',color='#337ba6',markersize=7)
    ax.axhline(floor,color='#555555',ls='--',label='Independent-electron reference, same family charge')
    ax.set_yscale('log');ax.set_xticks(x);ax.set_ylim(1e-6,2.)
    ax.set_xlabel('Harmonic of the optical phase');ax.set_ylabel('Mean conditional family |b|²')
    ax.set_title('Current input: 64 transverse families, 16 phase clones per family')
    ax.grid(axis='y',which='major',alpha=.2);ax.set_axisbelow(True)
    ax.legend(loc='upper left',bbox_to_anchor=(0.,.72),fontsize=9)
    fig.savefig(root/'input_family_harmonics.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main(sys.argv[1])
