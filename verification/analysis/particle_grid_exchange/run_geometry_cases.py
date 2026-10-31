#!/usr/bin/env python3
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np,yaml
import report_discrete,transport_review
root=Path(sys.argv[1]).resolve();assert os.environ.get('SLURM_JOB_ID')
plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text())
def check():
 for p,h in frozen.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
def save(p,d):p.write_text(json.dumps(d,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
check();all_results={}
for case in plan['cases']:
 directory=Path(case['directory']);env=os.environ.copy();env.pop('APRL_PARTICLE_WORK_GEOMETRY',None)
 env.update(OMP_NUM_THREADS=str(case['threads']),OMP_DYNAMIC='FALSE',OMP_PLACES='cores',OMP_PROC_BIND='close',APRL_PROFILE_INTERVAL='0')
 if case['mode']=='on':env['APRL_PARTICLE_WORK_GEOMETRY']=str(directory/'geometry')
 cmd=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks='+str(case['ranks']),'--cpus-per-task='+str(case['threads']),'--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30',case['binary'],str(directory/'run.yaml')]
 start=time.monotonic();print('RUN',case['name'],flush=True)
 with (directory/'aprl.log').open('x') as out:subprocess.run(['timeout','--signal=TERM','--kill-after=10',str(case['limit_seconds']),*cmd],cwd=directory,env=env,check=True,stdout=out,stderr=subprocess.STDOUT,timeout=case['limit_seconds']+20)
 seconds=time.monotonic()-start;log=(directory/'aprl.log').read_text()
 assert all(s in log for s in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete'])
 assert not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log)
 scientific=[]
 for path in sorted(directory.glob('detectors/*.h5'))+sorted(directory.glob('trajectory/*.h5'))+[directory/'ledger/energy-ledger.h5']:
  with h5py.File(path) as f:
   g=f[next(iter(f))];assert int(g['complete'][()])==1
   for committed,data in [('committed_records','records'),('committed_samples','time_s')]:
    if committed in g:assert int(g[committed][()])==len(g[data])
   if path.name.endswith('energy.h5'):continue
   if path.parent.name=='detectors' and 'records' in g:
    rows=g['records'][:];assert len(rows)==plan['particles'] and len(np.unique(rows['particle_id']))==plan['particles']
   for name,ds in g.items():
    if not isinstance(ds,h5py.Dataset):continue
    for offset in range(0,len(ds) if ds.shape else 1,8 if name=='fields' else 65536):
     x=ds[()] if not ds.shape else ds[offset:offset+(8 if name=='fields' else 65536)]
     if x.dtype.names:
      for n in x.dtype.names:
       if x[n].dtype.kind in 'fc':assert np.isfinite(x[n]).all()
     elif x.dtype.kind in 'fc':assert np.isfinite(x).all()
   scientific.append(str(path.relative_to(directory)))
 rows,sums,ranks,steps=report_discrete.load(directory/'ledger',directory/'ledger/energy-ledger.h5')
 plane_review=transport_review.plane_pair_review(directory,plan,yaml.safe_load((directory/'run.yaml').read_text()))
 assert plane_review['checks_pass']
 index,q,summary=report_discrete.calculate(rows,sums,steps)
 data={'step':index,'particle_grid_work_J':rows['particle_grid_work_J'],'constraints_J':rows['boundary_source_work_J']+rows['frozen_constraint_work_J']+rows['carrier_constraint_work_J'],
       'mesh_old_J':np.r_[np.longdouble(0),np.cumsum(sums['applied_current_old_J'],dtype=np.longdouble)][index],
       'mesh_mid_J':np.r_[np.longdouble(0),np.cumsum(sums['applied_current_mid_J'],dtype=np.longdouble)][index],
       'particle_change_J':rows['delta_particle_kinetic_box_J'],'mean_gamma_box':rows['mean_gamma_box']}
 np.savez(directory/'mesh_work.npz',**data)
 elapsed={int(s):float(t) for s,t in re.findall(r'\[resource\] phase=(?:progress|final) step=(\d+) elapsed_s=([\d.eE+\-]+)',log)}
 result={'integrity_passed':True,'steps':steps,'rank_count':ranks,'launcher_seconds':seconds,'progress_elapsed_s':elapsed,'strict_summary':summary,'energy_planes':plane_review,'scientific_files':scientific}
 if case['mode']=='on':
  paths=list((directory/'geometry').glob('*.csv'));assert len(paths)==ranks
  for path in paths:
   assert path.with_suffix('.complete').read_text().split()==['1',str(steps)]
   values=np.genfromtxt(path,delimiter=',',names=True);assert np.array_equal(values['step'],np.arange(1,steps+1))
   assert all(np.isfinite(values[n]).all() for n in values.dtype.names)
   assert np.all(values['eligible_pushes']+values['excluded_release_pushes']+values['excluded_boundary_interface_pushes']==values['observed_pushes'])
 save(directory/'RESULT.json',result);all_results[case['name']]=result
 print('COMPLETE',case['name'],steps,seconds,flush=True)
comparisons={}
for layout in ['native','mpi-small']:
 for kind in ['undulator','free-drift']:
  off=root/f'{layout}-{kind}-off';on=root/f'{layout}-{kind}-on'
  pairs=[(off,on)]
  if layout=='native':pairs.insert(0,(root/f'{layout}-{kind}-baseline',off))
  for first,last in pairs:
   item=transport_review.compare(first,last);assert item['checks_pass'],(first,last,item)
   with h5py.File(first/'ledger/energy-ledger.h5') as f,h5py.File(last/'ledger/energy-ledger.h5') as g:
    x=f['energy_ledger/records'][:];y=g['energy_ledger/records'][:]
   ledger_bits={n:x[n].tobytes()==y[n].tobytes() for n in x.dtype.names if 'seconds' not in n}
   item['ledger_fields_bitwise']=ledger_bits
   # Scientific state equality is strict. Pure summed diagnostics retain an
   # independently scaled floating reduction bound rather than feed state.
   for n,exact in ledger_bits.items():
    if not exact:
     bound=512*np.finfo(float).eps*np.maximum(abs(x[n])+abs(y[n]),1e-100)
     assert np.all(abs(x[n]-y[n])<=bound),n
   comparisons[first.name+'__'+last.name]=item
check();save(root/'COMPLETE.json',{'integrity_passed':True,'analysis_job':os.environ['SLURM_JOB_ID'],'physical_state_bitwise_pass':True,'comparisons':comparisons,'cases':all_results,'physics_method_changes':False,'scope':'Diagnostic correctness and unchanged physical outputs; not a coupling correction'})
print('ALL COMPLETE',flush=True)
