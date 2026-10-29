#!/usr/bin/env python3
"""Read only the closed 258-row free-drift ledger; no particle trajectory scan."""
import argparse
import csv
import json
from pathlib import Path
import h5py
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ledger', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a=p.parse_args(); a.output.mkdir(exist_ok=True,parents=True)
    with h5py.File(a.ledger) as f:
        g=f['energy_ledger']; r=g['records'][:]
        assert int(g['complete'][()])==1 and int(g['committed_records'][()])==len(r)
    assert len(r)==258 and np.array_equal(r['step'],np.arange(258))
    assert not np.any(r['prescribed_boris_work_J'])
    release=int(np.flatnonzero(r['released_charge_C']==r['released_charge_C'][-1])[0])
    first_exit=int(np.flatnonzero(r['removed_macroparticles']>0)[0])
    all_exit=int(np.flatnonzero(r['active_macroparticles']==0)[0])
    names=['field_exchange_residual_J','particle_work_balance_residual_J',
           'delta_particle_kinetic_box_J','delta_field_energy_J',
           'particle_grid_work_J','moving_remap_field_change_J',
           'inlet_overwrite_work_J','frozen_constraint_work_J','carrier_constraint_work_J',
           'outward_field_energy_J','boundary_source_work_J']
    assert all(n in r.dtype.names for n in names)
    stages={}
    for label,i,j in [('injection',0,release),('released_no_exit',release,first_exit-1),
                      ('particle_exit',first_exit-1,all_exit),('after_all_particle_exit',all_exit,257)]:
        stages[label]={'start_step':i,'end_step':j,
             'changes_J':{n:float(np.longdouble(r[n][j])-np.longdouble(r[n][i])) for n in names}}
    result={'ledger':str(a.ledger.resolve()),'rows':len(r),'frame':'boost',
            'all_released_step':release,'first_particle_exit_step':first_exit,'all_particles_exited_step':all_exit,
            'stages':stages,'final':{n:float(r[n][-1]) for n in names},
            'fraction_final_residual_accumulated_from_first_exit_interval':
              float((r['field_exchange_residual_J'][-1]-r['field_exchange_residual_J'][first_exit-1])/r['field_exchange_residual_J'][-1]),
            'fraction_final_residual_after_all_particles_exit':
              float((r['field_exchange_residual_J'][-1]-r['field_exchange_residual_J'][all_exit])/r['field_exchange_residual_J'][-1]),
            'interpretation':'Temporal localization only. It does not identify a causal numerical operator.'}
    (a.output/'DATA_DECOMPOSITION.json').write_text(json.dumps(result,indent=2)+'\n')
    with (a.output/'ledger_values.csv').open('w') as f:
        w=csv.writer(f); w.writerow(['step','time_box_s']+names)
        w.writerows([[int(x['step']),x['time_box_s']]+[x[n] for n in names] for x in r])
    print(json.dumps({k:v for k,v in result.items() if k not in ('stages','final')},indent=2))


if __name__=='__main__': main()
