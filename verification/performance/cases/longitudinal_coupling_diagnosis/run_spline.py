import os,json,concurrent.futures,time
from pathlib import Path
from prepare import sha,save,require
from launcher import launch
R=Path(__file__).resolve().parents[1]
def main():
 for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Frozen mismatch '+f)
 p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();cmd=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes=1','--ntasks=1','--cpus-per-task=1','--mem=1G','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=15'];reports={}
 def work(c):
  n=c['name'];print('START',n,time.time(),flush=True);launch(cmd+[p['program'],str(R/(n+'.csv')),*map(str,c['args'])],R,env,c['timeout_s'],n+'.log');d=json.loads((R/(n+'.log')).read_text());save(R/(n+'.json'),dict(configuration=c,result=d));print('COMPLETE',n,time.time(),flush=True);return n,d
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for x in concurrent.futures.as_completed([pool.submit(work,c) for c in p['cases']]):n,d=x.result();reports[n]=d;save(R/'PARTIAL.json',reports)
 for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Frozen mismatch '+f)
 save(R/'RESULT.json',dict(cases=reports,physical_acceptance=False,production_modified=False));print('COMPLETE',flush=True)
if __name__=='__main__':main()
