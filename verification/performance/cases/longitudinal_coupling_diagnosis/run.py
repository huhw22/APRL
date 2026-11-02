import os,json,subprocess,concurrent.futures,time
from pathlib import Path
from prepare import sha,save,require
from launcher import launch
R=Path(__file__).resolve().parents[1]
def main():
 for p,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(p)==h,'Frozen mismatch '+p)
 p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();cmd=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes=1','--ntasks=1','--cpus-per-task=1','--mem=2G','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=15'];launch(cmd+[p['programs']['axial_nondispersion']],R,env,60,'axial.log');axial=json.loads((R/'axial.log').read_text());require(axial['complete'],'Axial propagation failure');save(R/'AXIAL_RESULT.json',axial)
 def work(c):
  name=c['name'];print('START',name,time.time(),flush=True);launch(cmd+[p['programs']['periodic_diagnosis'],str(R/(name+'.csv')),*map(str,c['args'])],R,env,c['timeout_s'],name+'.log');lines=(R/(name+'.log')).read_text().splitlines();base=json.loads(next(x for x in lines if x.startswith('{')));diag=json.loads(next(x[len('LONGITUDINAL_DIAGNOSTICS '):] for x in lines if x.startswith('LONGITUDINAL_DIAGNOSTICS ')));result=dict(configuration=c,coupling=base,diagnosis=diag);save(R/(name+'.json'),result);print('COMPLETE',name,time.time(),flush=True);return name,result
 results={}
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for future in concurrent.futures.as_completed([pool.submit(work,c) for c in p['cases']]):
   name,result=future.result();results[name]=result;save(R/'PARTIAL.json',results)
 for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Frozen mismatch '+f)
 save(R/'RESULT.json',dict(cases=results,axial=axial,Cowan_propagator_changed=False,particle_force_changed=p["particle_force_changed"],production_replaced=False,physical_acceptance=False));print('COMPLETE',flush=True)
if __name__=='__main__':main()
