from decimal import Decimal as D,getcontext
from pathlib import Path
import json, hashlib
import h5py
getcontext().prec=60
import argparse
p=argparse.ArgumentParser(description='Independent Decimal audit of scalar plane module output')
p.add_argument('--unit',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
args=p.parse_args(); root=args.unit.resolve()
c=D(299792458);mu=D('1.25663706212e-6');eps=D('8.8541878128e-12');g=D(2);beta=D(3).sqrt()/2;area=D('.96');dt=D('.01')/c/g
power=[];density=[];integral=[D(0)]
for n in range(6):
 q=D('16.125')-D('.2')*beta*n
 a=1+D('.1')*q+D('.2')*n;b=2+D('.15')*q-D('.07')*n
 el=g*(a+beta*b);cb=g*(b+beta*a)
 power.append(area*el*cb/(mu*c));density.append(area/2*(eps*el*el+cb*cb/(mu*c*c)))
 if n: integral.append(integral[-1]+((power[n-1]+power[n])*dt/2 if n-1 in (0,3,4) else 0))
checks=[];artifacts={}
for ranks in (2,4):
 for cadence in ('dense','sparse'):
  path=root/f'energy-plane-{ranks}rank/{cadence}-energy.h5'
  with h5py.File(path) as f:
   group=f['energy_plane'];rows=group['records'][:]
   assert int(group['complete'][()])==1 and int(group['committed_records'][()])==len(rows)
  artifacts[str(path.relative_to(root))]={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(rows)}
  for row in rows:
   n=round(float(row['time_box_s'])/(.01/299792458))
   for prefix,factor in (('',D(1)),('physical_',D('.5'))):
    quantities=[('field_energy_lab_J',integral[n])]
    if row['field_available']:
     assert n != 2
     quantities += [('field_power_lab_W',power[n]),('field_energy_per_length_lab_J_per_m',density[n])]
    else: assert n==2
    for field,reference in quantities:
     expected=reference*factor;actual=D.from_float(float(row[prefix+field]))
     difference=abs(actual-expected);relative=difference/abs(expected) if expected else difference
     bound=D('3e-12')*max(abs(expected),D('1e-29'))
     checks.append(dict(ranks=ranks,cadence=cadence,step=n,field=prefix+field,actual=str(actual),reference=str(expected),
       absolute_error=str(difference),relative_error=str(relative),pass_check=difference<=bound))
result={'status':'PASS' if all(x['pass_check'] for x in checks) else 'FAIL','checks':checks,'check_count':len(checks),
 'maximum_relative_error':str(max(D(x['relative_error']) for x in checks)),'artifacts':artifacts,
 'method':'60-digit Decimal, explicit inverse Lorentz formulas, independent area and allowed trapezoid intervals; no APRL transforms or sampler invoked',
 'scope':'All saved available power and J/m samples plus all cumulative field-energy records, full and physical aperture, dense/sparse rhythm, 2/4 MPI ranks',
 'comparison_bound':'3e-12*max(abs(reference),1e-29), fixed module tolerance; observed error reported, not inflated to the bound'}
args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(result['status'],result['check_count'],result['maximum_relative_error'])
assert result['status']=='PASS'
