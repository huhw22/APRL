#!/usr/bin/env python3
"""Read a small deterministic sample; never alter the production particles."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--input',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
if a.out.exists():
    raise RuntimeError('Use a fresh output path')
with h5py.File(a.input,'r') as f:
    d=f['particles/records']
    count=len(d)
    if count<65536:
        raise RuntimeError('This sampling convention expects the original large input')
    starts=np.linspace(0,count-1024,64,dtype=np.int64)
    records=np.concatenate([d[int(s):int(s)+1024] for s in starts])
u=records['proper_velocity'];xy=records['plane_position_m'];w=records['macro_weight']
g=np.sqrt(1+(u*u).sum(axis=1))
mean=np.average(u,axis=0,weights=w)
result={'format':'aprl-held-cut-source-moment-sample-v1','input':str(a.input.resolve()),
        'input_records':count,'sample_records':len(records),
        'sample_method':'64 evenly spaced blocks of1024 records,read-only deterministic sample; sample moments,not full-population guarantee',
        'block_starts':starts.tolist(),'lab_u_mean':mean.tolist(),
        'lab_u_sigma':np.sqrt(np.average((u-mean)**2,axis=0,weights=w)).tolist(),
        'lab_gamma_mean':float(np.average(g,weights=w)),
        'lab_gamma_sigma':float(np.sqrt(np.average((g-np.average(g,weights=w))**2,weights=w))),
        'xy_sigma_m':np.std(xy,axis=0).tolist(),
        'u_x_y_z_correlations':np.corrcoef(u.T).tolist(),
        'x_ux_correlation':float(np.corrcoef(xy[:,0],u[:,0])[0,1]),
        'y_uy_correlation':float(np.corrcoef(xy[:,1],u[:,1])[0,1])}
a.out.write_text(json.dumps(result,indent=2)+'\n')
