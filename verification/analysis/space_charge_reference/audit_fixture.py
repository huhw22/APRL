#!/usr/bin/env python3
"""Read entry-event moments; do not identify them with a common-time cold bunch."""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
import yaml

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--card',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();cfg=yaml.safe_load(a.card.read_text());path=Path(cfg['beam']['input']['file'])
with h5py.File(path,'r') as h:data=h['particles/records'][:]
u=np.asarray(data['proper_velocity'],dtype=np.longdouble)
gamma=np.sqrt(1+np.sum(u*u,axis=1));boost=np.longdouble(cfg['mesh']['boost_gamma'])
gb=np.sqrt((boost-1)*(boost+1));plus=(gamma+u[:,2])/(boost+gb)
minus=(1+np.sum(u[:,:2]**2,axis=1))/plus;gp=(plus+minus)/2
up=np.c_[u[:,:2],(plus-minus)/2];v=up/gp[:,None]*299792458
x=np.asarray(data['plane_position_m'],dtype=np.longdouble);metrics={}
for i,name in enumerate(('x','y')):
    dx=x[:,i]-np.mean(x[:,i]);dv=v[:,i]-np.mean(v[:,i]);cov=np.mean(dx*dv)
    metrics[name]={'sigma_entry_position_m':float(np.sqrt(np.mean(dx*dx))),
      'correlation_entry_position_boost_velocity':float(cov/np.sqrt(np.mean(dx*dx)*np.mean(dv*dv))),
      'conditional_mean_velocity_slope_per_s':float(cov/np.mean(dx*dx))}
result={'scope':'Entry-event moments of existing fixed-plane functional fixture; not a common-time isolated cold Gaussian or an explanation of a measured energy change',
 'input':str(path),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'particles':len(data),
 'represented_electrons':cfg['beam']['input']['electrons'],'boost_gamma':float(boost),
 'mean_boost_longitudinal_beta':float(np.mean(v[:,2])/299792458),'entry_moments':metrics,
 'current_interface':{'moving_window_requires_frozen_plane_injection':True,'strict_discrete_audit_requires_moving_window':True},'new_APRL_jobs':0}
a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
