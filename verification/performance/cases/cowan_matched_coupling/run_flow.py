"""Run bounded short full-program controls, then read closed evidence only."""
import json,os,subprocess,time
from pathlib import Path
from prepare import require,sha,save
from launcher import launch,binding
R=Path(__file__).resolve().parents[1]
def main():
 for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Changed frozen '+f)
 p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();srun=['srun','--mpi=pmix_v3','--exact','--nodes=2','--ntasks=16','--ntasks-per-node=8','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];launch(srun+[p['binding_probe'],'16','8'],R,env,45,'binding.log');save(R/'BINDING.json',binding(R/'binding.log',2))
 for case in p['cases']:
  d=R/case['name'];print('START',case['name'],time.time(),flush=True);ev=env.copy();ev.update(APRL_LAB_FIELD_AUDIT=str(d/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(d/'lab-audit'),APRL_PARTICLE_MOMENTUM_PROBE=str(d/'momentum'),APRL_PARTICLE_MOMENTUM_PROBE_DT='4.1839819784122236e-15',APRL_PARTICLE_MOMENTUM_PROBE_STRIDE='1');launch(srun+[p['programs']['wave'],str(d/'run.yaml')],d,ev,case['timeout_seconds'],'aprl.log');text=(d/'aprl.log').read_text();require(all(x in text for x in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),'Case did not close normally');print('COMPLETE',case['name'],time.time(),flush=True)
 save(R/'PIC_COMPLETE.json',dict(passed=True,candidate_algorithm_changed=True));save(R/'CALCULATION_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID']));small=env.copy();small['OMP_NUM_THREADS']='1';subprocess.run([os.environ['PYTHON_BIN'],str(R/'frozen_source/post.py')],env=small,check=True);subprocess.run([os.environ['PYTHON_BIN'],str(R/'frozen_source/review.py'),str(R)],env=small,check=True)
 for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Changed frozen '+f)
if __name__=='__main__':main()
