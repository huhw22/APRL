import json,os,subprocess,time
from pathlib import Path
from prepare import require,save,sha
from launcher import launch,binding
R=Path(__file__).resolve().parents[1]
def main():
    for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(Path(f))==h,'Frozen input changed '+f)
    p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();test=env.copy();test['OMP_NUM_THREADS']='1'
    with open(R/'GEOMETRY_CHECK.json','w') as out:subprocess.run(['timeout','60',p['programs']['check']],stdout=out,stderr=subprocess.STDOUT,env=test,check=True)
    check=json.loads((R/'GEOMETRY_CHECK.json').read_text());require(check['energy_surface_SBP_relative']<1e-12 and check['finite_step_H_relative']<1e-12 and check['finite_step_P_relative']<1e-12 and check['stress_time_pairing_relative']<1e-12,'Independent geometry failed')
    srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(srun+[p['binding_probe'],'16','8'],R,env,45,'binding.log');save(R/'BINDING.json',binding(R/'binding.log',1))
    for case in p['cases']:
        d=R/case['name'];d.mkdir();print('START',case['name'],time.time(),flush=True)
        launch(srun+[p['programs']['wave'],str(d),'yee',case['mode'],str(case['refine'])],d,env,case['timeout_seconds'],'run.log');require((d/'complete').exists(),'Incomplete case');print('COMPLETE',case['name'],time.time(),flush=True)
    save(R/'CALCULATION_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID']))
    subprocess.run([os.environ['PYTHON_BIN'],str(R/'frozen_source/post.py')],check=True,env=test)
if __name__=='__main__':main()
