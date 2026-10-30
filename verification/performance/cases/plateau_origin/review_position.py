#!/usr/bin/env python3
"""Bounded estimate of smooth-field change caused by orbit displacement."""
import csv,json,sys
from pathlib import Path
import numpy as np
from scipy.constants import e,m_e,c
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(bench):
    bench=Path(bench).resolve();late=bench/'charge_inlet_plateau_late_force_20261030/guard64_1000'
    dense=bench/'charge_inlet_plateau_dense_force_20261030_r2'
    f=np.genfromtxt(late/'momentum-analysis/measured_force.csv',names=True,delimiter=',')
    m=(f['step']>=500)&(f['step']<=980);step=f['step'][m]
    # E_ref(actual PIC event) - E_ref(the same source on its ideal sphere worldline).
    effect=f['mean_Ez_error_worldline_position_V_per_m'][m]-f['mean_Ez_error_V_per_m'][m]
    assert len(effect)==25 and np.all(np.diff(step)==20)
    nominal_dt=4.1839819784122236e-15;g=929.1305492858543;beta=np.sqrt(1-g**-2)
    impulse=-e*nominal_dt/(m_e*c)*np.trapz(effect,step)
    approximate_lab_pz=g*beta*m_e*c*c/e*impulse
    result=json.loads((dense/'review/RESULT.json').read_text())
    late_interval=next(r for r in result['intervals'] if r['from_step']==500 and r['to_step']==1000)
    out=dense/'review'
    rows=[dict(step=int(s),reference_field_at_pic_minus_ideal_worldline_V_per_m=float(v)) for s,v in zip(step,effect)]
    with (out/'position_feedback.csv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(8,4.2),layout='constrained')
    ax.plot(step,effect/1000,'.-',color='#b34562')
    ax.set(xlabel='Field step',ylabel='Mean smooth-reference Ez difference (kV/m)',
           title='Reference force changes as PIC orbits depart from ideal worldlines')
    ax.grid(alpha=.25)
    fig.savefig(out/'position_feedback.png',dpi=180);fig.savefig(out/'position_feedback.pdf');plt.close(fig)
    payload=dict(scope='Sparse 20-step diagnostic; same source IDs and rest time, actual PIC position vs ideal worldline position.',
                 samples=len(rows),step_range=[500,980],monotonic_samples=bool(np.all(np.diff(effect)>0)),
                 first_V_per_m=float(effect[0]),last_V_per_m=float(effect[-1]),
                 trapezoid_lab_pz_eV_500_to_980=float(approximate_lab_pz),
                 exact_remainder_lab_pz_eV_500_to_1000=late_interval['orbit_density_reference_remainder_lab_momentum_eV'],
                 limitation='Not an exact decomposition of the 500-1000 remainder: samples omit the final 20 steps and can alias field oscillations; reference trajectory and numerical time-discretization effects remain.')
    (out/'POSITION_FEEDBACK.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps(payload,indent=2))
if __name__=='__main__':main(sys.argv[1])
