#!/usr/bin/env python3
"""Compare post-release changes without conflating their initial offsets."""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(bench):
    bench=Path(bench).resolve();root=bench/'charge_inlet_plateau_origin_20261030_r3'
    geometry=np.genfromtxt(bench/'charge_inlet_plateau_wide_xy_20261030/review/aperture_vs_resolution.csv',names=True,delimiter=',')
    guards=np.genfromtxt(root/'comparison.csv',names=True,delimiter=',')
    curves=[('128-cell guard, base',geometry,'base_lab_error_eV'),
            ('256-cell guard',guards,'guard256_lab_error_eV'),
            ('128-cell guard, finer xy',geometry,'fine_dx_lab_error_eV'),
            ('128-cell guard, wider xy',geometry,'wide_aperture_lab_error_eV')]
    rows=[];fig,ax=plt.subplots(figsize=(8.5,4.8),layout='constrained')
    for label,a,key in curves:
        a=a[(a['step']>=123)&(a['step']<=300)]
        value={int(r['step']):float(r[key]) for r in a}
        rows.append(dict(case=label,energy_error_123_eV=value[123],energy_error_180_eV=value[180],
                         change_123_to_180_eV=value[180]-value[123],change_123_to_300_eV=value[300]-value[123]))
        ax.plot(a['step'],(a[key]-value[123])/1000,'.-',label=label)
    ax.set(xlabel='Field step',ylabel='Change from step 123 (keV/electron)',
           title='Similar post-release transient despite different initial offsets')
    ax.grid(alpha=.2);ax.legend()
    fig.savefig(root/'figures/common_transient.png',dpi=180);fig.savefig(root/'figures/common_transient.pdf');plt.close(fig)
    with (root/'common_transient.csv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (root/'COMMON_TRANSIENT.json').write_text(json.dumps({'rows':rows,'scope':'Changes on a common boost-time slice. Similarity does not identify a unique numerical cause.'},indent=2)+'\n')
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main(sys.argv[1])
