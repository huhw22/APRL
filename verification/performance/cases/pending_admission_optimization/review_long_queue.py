#!/usr/bin/env python3
"""Read-only independent audit of already completed long-queue HDF5 outputs."""
import argparse, json, pathlib, hashlib, re
import numpy as np
import h5py
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=pathlib.Path,required=True,help='Completed 8x-arrival long-queue campaign; peer mpi-r1/pair-r1 evidence must be siblings')
P=parser.parse_args().root.resolve()
N=31675;EPS=np.finfo(float).eps
report={'scope':'Independent raw HDF5 reread of field plane, all three particle planes, energy ledger and terminal trajectory events; no rerun or threshold change. Interior trajectory comparison is not claimed.', 'pass':True}
def require(v,m):
 if not v: raise RuntimeError(m)
def metrics(a,b):
 x=np.asarray(a,dtype=np.longdouble);y=np.asarray(b,dtype=np.longdouble);d=y-x
 return {'raw_bytes_equal':a.dtype==b.dtype and a.tobytes()==b.tobytes(),'max_abs_difference':float(np.max(abs(d))), 'rms_difference':float(np.sqrt(np.mean(d*d)))}
def ptmetrics(a,b):
 r=metrics(a,b);scale=np.max(abs(a),axis=0);bound=128*EPS*scale
 r['elements_exceeding_128eps_component_gate']=int(np.count_nonzero(abs(np.longdouble(b)-np.longdouble(a))>bound))
 require(r['elements_exceeding_128eps_component_gate']==0,'Particle roundoff gate failed')
 return r
def check_group(g,ds,commit):
 require(int(g['complete'][()])==1,'Incomplete output')
 require(int(g[commit][()])==len(g[ds]),'Uncommitted extent')
report['fields']={}
with h5py.File(P/'baseline/detectors/exit-field.h5') as fa,h5py.File(P/'candidate/detectors/exit-field.h5') as fb:
 a,b=fa['field_plane'],fb['field_plane']
 for g in [a,b]:check_group(g,'fields','committed_samples')
 require(np.array_equal(a['time_s'][:],b['time_s'][:]),'Field sample times differ')
 require(a['fields'].shape==b['fields'].shape==(717,160,160),'Field shape mismatch')
 totals={k:[np.longdouble(0),np.longdouble(0),0.,True] for k in ['electric_V_per_m','magnetic_T']}
 for start in range(0,717,8):
  aa=a['fields'][start:start+8];bb=b['fields'][start:start+8]
  for k,t in totals.items():
   x,y=aa[k],bb[k];require(np.isfinite(x).all() and np.isfinite(y).all(),'Nonfinite field');d=np.longdouble(y)-np.longdouble(x)
   t[0]+=np.sum(d*d,dtype=np.longdouble);t[1]+=np.sum(np.longdouble(x)**2,dtype=np.longdouble);t[2]=max(t[2],float(np.max(abs(d))));t[3]&=x.tobytes()==y.tobytes()
 for k,(ds,rs,mx,eq) in totals.items():
  r=float(np.sqrt(ds/rs));require(r<=2e-12,'Field gate failed');report['fields'][k]={'relative_l2':r,'max_abs_difference':mx,'raw_bytes_equal':eq,'gate':2e-12}
report['planes']={}
for plane in ['release','aux','exit']:
 arr=[]
 for tag in ['baseline','candidate']:
  with h5py.File(P/tag/'detectors'/f'{plane}.h5') as f:
   g=f['particle_plane'];check_group(g,'records','committed_records');x=g['records'][:];x=x[np.argsort(x['source_id'])];require(np.array_equal(x['source_id'],np.arange(1,N+1)),'ID coverage failure');arr.append(x)
 a,b=arr
 for name in ['source_id','particle_id','charge_C','mass_kg','weight','crossing_direction']:require(a[name].tobytes()==b[name].tobytes(),'Particle identity differs')
 report['planes'][plane]={'count':len(a),'identity_exact':True,**{n:ptmetrics(a[n],b[n]) for n in ['position_m','proper_velocity','time_s']}}
report['ledger']={};arr=[]
for tag in ['baseline','candidate']:
 with h5py.File(P/tag/'ledger/energy-ledger.h5') as f:
  g=f['energy_ledger'];check_group(g,'records','committed_records');x=g['records'][:];require(len(x)==717,'ledger length');require(x['active_macroparticles'][0]==N,'initial pending count omitted');require(np.all(x['active_macroparticles']+x['removed_macroparticles']==N),'cohort incomplete');arr.append(x)
a,b=arr
for n in a.dtype.names:
 require(np.isfinite(a[n]).all() and np.isfinite(b[n]).all(),'Nonfinite ledger')
 if a[n].dtype.kind not in 'fc' or n=='time_box_s':require(np.array_equal(a[n],b[n]),'Integer/time ledger differs')
 report['ledger'][n]=metrics(a[n],b[n])
report['terminal_events']={};arr=[]
for tag in ['baseline','candidate']:
 blocks=[];records=0
 for path in sorted((P/tag/'trajectory').glob('*.h5')):
  with h5py.File(path) as f:
   g=f['trajectory'];check_group(g,'records','committed_records');d=g['records'];records+=len(d)
   for start in range(0,len(d),32768):
    x=d[start:start+32768];x=x[x['event_type']!=0]
    if len(x):blocks.append(x)
 x=np.concatenate(blocks);x=x[np.argsort(x['source_id'])]
 require(np.array_equal(x['source_id'],np.arange(1,N+1)),'Terminal IDs incomplete');require(np.all(x['event_type']==2) and np.all(x['boundary_face']==5),'Unexpected terminal event');arr.append(x)
 report['terminal_events'][tag]={'all_trajectory_records':records,'terminal_records':len(x),'one_right_exit_per_id':True}
a,b=arr
report['terminal_events']['comparison']={n:ptmetrics(a[n],b[n]) for n in ['position_m','proper_velocity','time_s']}
report['runtime']={}
for tag in ['baseline','candidate']:
 text=(P/tag/'aprl.log').read_text();m=re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+)',text);require(tuple(map(int,m.groups()))==(4096,289,2274),'Runtime geometry mismatch')
 require('Configured stop reached:' in text and 'phase=mpi-finalize-complete' in text,'Incomplete normal shutdown')
 m=re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)',text);require(tuple(map(int,m.groups()))==(716,N,25688),'Incomplete lifecycle')
 report['runtime'][tag]={'steps':716,'particle_exits':N,'cropped_carriers':25688,'reference_nz':4096,'field_nz':289,'initial_offset':2274,'normal_stop_and_finalize':True}
text=(P/'candidate/aprl.log').read_text();pending=int(re.search(r'Pending admission queue: waiting=(\d+)',text)[1]);require(pending==31669,'pending unexpected');report['initial_pending']=pending
# Cross-host comparison has long-double temporary-padding flags in sigma gate;
# inspect actual serialized float64 fields before stating bitwise agreement.
report['two_node_raw_ledger']={}
mpi=P.parent/'mpi-r1';short=P.parent/'pair-r1'
for case in ['free-drift','undulator']:
 with h5py.File(mpi/case/'ledger/energy-ledger.h5') as fa,h5py.File(short/'science'/f'{case}-16x8-candidate'/'ledger/energy-ledger.h5') as fb:
  a,b=fa['energy_ledger/records'][:],fb['energy_ledger/records'][:]
  report['two_node_raw_ledger'][case]={'all_fields_bytes_equal':all(a[n].tobytes()==b[n].tobytes() for n in a.dtype.names),'nonidentical_fields':[n for n in a.dtype.names if a[n].tobytes()!=b[n].tobytes()]}
report['proof_hashes']={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['LONG_QUEUE_PASS.json','COMPARISON.json','THRESHOLDS.json','GEOMETRY.json','frozen.json']}
(P/'INDEPENDENT_REVIEW.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'pass':report['pass'],'fields':report['fields'],'initial_pending':pending,'two_node_raw_ledger':report['two_node_raw_ledger'],'output':str(P/'INDEPENDENT_REVIEW.json')},indent=2))
