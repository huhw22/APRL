#!/usr/bin/env python3
"""Compare the original field code with an independently translated exact pulse."""
import argparse
import csv
import json
import os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-matplotlib-free-work')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LD=np.longdouble


def exact_diagnostic(start,width):
    # Closed translated profiles, not a replay of the Maxwell update loops.
    nz=64; k=np.arange(nz+1,dtype=LD)
    eps=LD('8.8541878128e-12');mu=LD('1.25663706212e-6');c=LD('299792458')
    dz=LD('1.2543262415359034e-6');dt=dz/c;area=(16*LD('5e-6'))**2
    def pulse(x):
        phase=(x-start)/width
        return np.where((phase>0)&(phase<1),LD('1e6')*np.sin(np.arccos(LD(-1))*phase)**2,LD(0))
    def state(n):return pulse(k-2*n),pulse(k[:-1]+1-2*n)/c
    def energy(e,b):return area*dz/2*np.sum(eps*((e[:-1]+e[1:])/2)**2+b*b/mu)
    def power(e,b):return area/mu*((e[-2]+e[-1])/2*b[-1]-(e[0]+e[1])/2*b[0])
    e,b=state(0);initial=energy(e,b);q=LD(0);remap=LD(0)
    values=[(0,initial,0.,0.,0.)]
    for n in range(28):
        e,b=state(n);new_e=pulse(k-1-2*n);new_e[0]=0;new_e[-1]=e[-1]
        new_b=pulse(k[:-1]-2*n)/c
        q+=dt/2*(power(e,b)+power(new_e,new_b))
        shifted_e,shifted_b=state(n+1)
        remap+=energy(shifted_e,shifted_b)-energy(new_e,new_b)
        u=energy(shifted_e,shifted_b);values.append((n+1,u,q,remap,u-initial+q-remap))
    return np.array(values,dtype=LD)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve()
    assert (root/'EXIT_CODE').read_text().strip()=='0'
    rows=[json.loads(x) for x in (root/'probe.jsonl').read_text().splitlines() if x.startswith('{')]
    assert len(rows)==4 and all(x['pass'] and x['particles']==0 and x['current']==0 for x in rows)
    results=[];fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for r in rows:
        path=root/f"pulse-start-{r['start']}-width-{r['width_cells']}.csv"
        with path.open() as f:trace=list(csv.DictReader(f))
        actual=np.array([[LD(x[n]) for n in ('step','field_J','outward_J','remap_J','residual_J')] for x in trace])
        expected=exact_diagnostic(r['start'],r['width_cells'])
        assert len(actual)==29 and np.isfinite(actual).all()
        initial=LD(str(r['initial_energy_J']))
        errors=np.max(abs(actual[:,1:]-expected[:,1:]),axis=0)
        assert np.all(errors<abs(initial)*LD('1e-8')),errors
        results.append({**r,'exact_translation_model_residual_J':float(expected[-1,-1]),
                        'maximum_ledger_difference_from_exact_translation_model_over_initial':float(np.max(errors)/abs(initial)),
                        'scope':'Same mesh and operator, different pulse widths/positions; not a fixed-physics PIC convergence study.'})
        ax.plot(actual[:,0],100*actual[:,-1]/initial,label=f"start {r['start']}, width {r['width_cells']} cells")
        ax.plot(expected[:,0],100*expected[:,-1]/initial,color='black',linestyle=':',linewidth=.8)
    result={'checks_pass':True,'results':results,'frame':'boost',
        'model':'E(k,n)=pulse(k-2n), cB_y(k,n)=pulse(k+1-2n); evaluated without Maxwell iteration',
        'diagnostic':'Same cell-centred energy, time-trapezoid cell-centred Poynting proxy and measured remap energy as the production ledger',
        'conclusion':'A nonzero reported energy residual is reproduced by the exact translated pulse; it is not by itself proof of incorrect wave propagation.',
        'limits':['Synthetic 1-D vacuum test; it does not establish correctness of all 3-D radiation or of the moving inlet.',
                  'The tested central aperture is causally separated from transverse CPML over the finite numerical stencil.',
                  'Pulse-width variation changes the physical pulse; it is not called fixed-physics grid convergence.']}
    (root/'VACUUM_REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    ax.axhline(0,color='#888888',linewidth=.7);ax.set_xlabel('Moving-window step')
    ax.set_ylabel('Reported residual / initial field energy (%)')
    ax.set_title('Original Cowan/crop code: exact vacuum propagation with a nonzero ledger residual')
    ax.grid(alpha=.2);ax.legend(fontsize=9,title='Black dots: exact translation model',title_fontsize=9)
    fig.savefig(root/'vacuum_ledger_residual.png',dpi=180);fig.savefig(root/'vacuum_ledger_residual.pdf');plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
