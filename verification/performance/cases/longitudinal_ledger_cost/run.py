import csv,json,os,re,sys,time
from pathlib import Path
import numpy as np,h5py
from prepare import sha,save,require
from launcher import launch
from binding import binding

def common_outputs(first,second):
 out={}
 for p in sorted((first/'detectors').glob('*.h5')):
  with h5py.File(p) as h:
   root=h[next(iter(h))]
   if 'energy_plane' in h:continue
   with h5py.File(second/'detectors'/p.name) as z:
    other=z[next(iter(z))];require(int(root['complete'][()])==int(other['complete'][()])==1,'Output not closed')
    for name in ['records','fields','time_s','coordinates']:
     if name not in root:continue
     a,b=root[name],other[name];require(a.shape==b.shape and a.dtype==b.dtype,'Output shape changed '+str(p))
     for k in range(0,len(a),8 if name=='fields' else 65536):
      x=a[k:k+(8 if name=='fields' else 65536)];y=b[k:k+(8 if name=='fields' else 65536)]
      require(x.tobytes()==y.tobytes(),'Scientific output changed '+str(p)+' '+name)
    out[p.name]=True
 return out

def main():
 r=Path(sys.argv[1]).resolve();p=json.loads((r/'PLAN.json').read_text());frozen=json.loads((r/'FROZEN.json').read_text());require(os.environ.get('SLURM_JOB_ID'),'Compute allocation required');require(all(sha(q)==h for q,h in frozen.items()),'Frozen input changed');results={};env=os.environ.copy();env['APRL_PROFILE_INTERVAL']='0'
 for key in ['APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_WORK_GEOMETRY','APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT','APRL_INLET_CONSTRAINT_PROBE','APRL_SOURCE_EXIT_AUDIT','APRL_DIAGNOSTIC_STEPS']:env.pop(key,None)
 sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
 for c in p['cases']:
  d=Path(c['directory']);print('START',c['name'],time.time(),flush=True);launch([*sr,str(Path(p['workspace'])/'build/cluster-tools/binding_probe'),'16','8'],d,env,60,'binding.log');b=binding(d/'binding.log',1,16,8);save(d/'BINDING.json',b)
  launch([*sr,p['binary'],str(d/'run.yaml')],d,env,c['limit_seconds']);log=(d/'aprl.log').read_text();require(all(s in log for s in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),'Normal completion missing');require(not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT',log),'Runtime failure')
  points={int(n):float(t) for n,t in re.findall(r'phase=progress step=(\d+) elapsed_s=([^ ]+)',log)};require(80 in points and 160 in points,'Timing interval missing');seconds=(points[160]-points[80])/80
  for f in (d/'detectors').glob('*.h5'):
   with h5py.File(f) as h:
    g=h[next(iter(h))];require(int(g['complete'][()])==1,'Detector incomplete');
    if 'particle_plane' in h:
     a=g['records'][:];require(len(a)==len(np.unique(a['particle_id']))==31675,'Particle coverage changed')
  iterations=[]
  for f in (d/'reference').glob('source-coupling-rank-*.csv'):
   require(f.with_suffix('.complete').exists(),'Coupling log incomplete');a=np.genfromtxt(f,delimiter=',',names=True);require(np.all(a['self_force_defect']<1e-13) and np.all(a['relative_fixed_point_work_defect']<1e-8),'Coupling tolerance changed');iterations.append(float(np.max(a['iterations'])))
  result=dict(mode=c['mode'],repeat=c['repeat'],seconds_per_step=seconds,interval=[80,160],max_iterations=max(iterations),binding=b,closed=True);save(d/'RESULT.json',result);results[c['name']]=result;save(r/'PARTIAL.json',results);print('COMPLETE',c['name'],seconds,flush=True)
  # Compare every closed configuration to exactly the same first off run.
  if c['name']!='r0-off':save(d/'BITWISE.json',common_outputs(r/'r0-off',d))
 summary=[]
 for repeat in [0,1]:
  t={mode:results[f'r{repeat}-{mode}']['seconds_per_step'] for mode in ['off','basic','strict','all']};den=t['all'];summary.append(dict(repeat=repeat,**t,basic_increment_s=t['basic']-t['off'],strict_increment_s=t['strict']-t['basic'],planes_increment_s=t['all']-t['strict'],all_increment_s=t['all']-t['off'],ledger_fraction_of_all=(t['all']-t['off'])/den,ledger_increase_over_off=(t['all']-t['off'])/t['off']))
 def stats(key):
  values=[v[key] for v in summary];return dict(median=float(np.median(values)),minimum=min(values),maximum=max(values))
 final=dict(passed=True,job_id=os.environ['SLURM_JOB_ID'],same_binary=True,same_physical_outputs_bitwise=True,runs=results,paired=summary,statistics={key:stats(key) for key in ['off','basic','strict','all','basic_increment_s','strict_increment_s','planes_increment_s','all_increment_s','ledger_fraction_of_all','ledger_increase_over_off']},scope='Two reverse-order paired repeats; one-node 16MPIx8 short drift, production ledger cadence. Not a production-run fraction.',production_jobs_changed=False)
 require(all(sha(q)==h for q,h in frozen.items()),'Frozen input changed');save(r/'RESULT.json',final)
 with (r/'COST.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
 print(json.dumps(final['statistics'],indent=2),flush=True)
if __name__=='__main__':main()
