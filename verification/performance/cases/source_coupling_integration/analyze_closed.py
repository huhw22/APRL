import json
from pathlib import Path
import h5py,numpy as np
def verify(directory):
 d=Path(directory);files=sorted(d.glob('detectors/*.h5'))+sorted(d.glob('trajectory/*.h5'))+[d/'ledger/energy-ledger.h5'];result={}
 for p in files:
  with h5py.File(p) as h:
   g=h[next(iter(h))];assert int(g['complete'][()])==1,p
   for key,values in [('committed_records','records'),('committed_samples','time_s')]:
    if key in g:assert int(g[key][()])==len(g[values]),p
   def visit(name,ds):
    if not isinstance(ds,h5py.Dataset):return
    if ds.shape and not ds.size:return
    for start in range(0,len(ds) if ds.shape else 1,8 if name=='fields' else 65536):
     x=ds[start:start+(8 if name=='fields' else 65536)] if ds.shape else ds[()]
     if x.dtype.names:
      for n in x.dtype.names:
       if x[n].dtype.kind in 'fc':assert np.isfinite(x[n]).all(),(p,n)
     elif x.dtype.kind in 'fc':assert np.isfinite(x).all(),(p,name)
   g.visititems(visit);result[str(p.relative_to(d))]=True
 return result
def bitwise(first,second):
 a=Path(first);b=Path(second);result={}
 for p in sorted(a.glob('detectors/*.h5'))+sorted(a.glob('trajectory/*.h5'))+[a/'ledger/energy-ledger.h5']:
  q=b/p.relative_to(a)
  with h5py.File(p) as h,h5py.File(q) as z:
   g=h[next(iter(h))];f=z[next(iter(z))]
   def visit(name,ds):
    if not isinstance(ds,h5py.Dataset) or name not in ['records','fields','time_s','coordinates']:return
    other=f[name];assert ds.shape==other.shape and ds.dtype==other.dtype,(p,name)
    for start in range(0,len(ds),8 if name=='fields' else 65536):
     x=ds[start:start+(8 if name=='fields' else 65536)];y=other[start:start+(8 if name=='fields' else 65536)]
     if x.dtype.names:
      for n in x.dtype.names:
       if 'seconds' not in n:assert np.ascontiguousarray(x[n]).tobytes()==np.ascontiguousarray(y[n]).tobytes(),(p,name,n)
     else:assert x.tobytes()==y.tobytes(),(p,name)
   g.visititems(visit);result[str(p.relative_to(a))]=True
 return result
