import json,os,subprocess,time
from pathlib import Path
from prepare import require,sha,save
from launcher import launch,binding
R=Path(__file__).resolve().parents[1]
def main():
    for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Changed frozen '+f)
    p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();small=env.copy();small['OMP_NUM_THREADS']='1';srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];one=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=1','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(one+[p['programs']['check']],R,small,60,'gather-check.log');check=json.loads((R/'gather-check.log').read_text());require(check['adjoint_relative_max']<2e-12 and check['magnetic_no_work_normalized_max']<1e-14,'Magnetic adjoint check');save(R/'GATHER_CHECK.json',check)
    launch(srun+[p['binding_probe'],'16','8'],R,env,45,'binding.log');save(R/'BINDING.json',binding(R/'binding.log',1))
    for case in p['cases']:
        d=R/case['name'];print('START',case['name'],time.time(),flush=True);e=env.copy();e.update(APRL_LAB_FIELD_AUDIT=str(d/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(d/'lab-audit'),APRL_PARTICLE_MOMENTUM_PROBE=str(d/'momentum'),APRL_PARTICLE_MOMENTUM_PROBE_DT='4.1839819784122236e-15',APRL_PARTICLE_MOMENTUM_PROBE_STRIDE='1');launch(srun+[p['programs']['wave'],str(d/'run.yaml')],d,e,case['timeout_seconds'],'aprl.log');print('COMPLETE',case['name'],time.time(),flush=True)
    save(R/'PIC_COMPLETE.json',dict(passed=True,candidate_algorithm_changed=True));save(R/'CALCULATION_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID']));subprocess.run([os.environ['PYTHON_BIN'],str(R/'frozen_source/post.py')],env=small,check=True);subprocess.run([os.environ['PYTHON_BIN'],str(R/'frozen_source/review.py'),str(R)],env=small,check=True)
if __name__=='__main__':main()
