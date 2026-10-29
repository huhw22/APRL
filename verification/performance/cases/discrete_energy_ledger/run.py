#!/usr/bin/env python3
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent/'transport'))
from prepare import require,save,sha
from run import launch,binding
import review
sys.path.insert(0,str(Path(__file__).parent))
import review_discrete

def main():
    root=Path(sys.argv[1]).resolve();plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'frozen.json').read_text())
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NUM_NODES')=='1','Require one real Slurm node')
    require(os.environ.get('SLURM_NTASKS')=='16' and os.environ.get('SLURM_CPUS_PER_TASK')=='8','Require16x8')
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen campaign changed')
    execution=root/'execution';execution.mkdir();env=os.environ.copy()
    for flag in ('APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT'):env.pop(flag,None)
    env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    def srun(ranks=16):return ['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1','--ntasks='+str(ranks),'--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    try:
        launch(srun()+[plan['binding_probe'],'16','8'],execution,env,45,'binding.log');save(execution/'BINDING.json',binding(execution/'binding.log',1))
        single=env.copy();single['OMP_NUM_THREADS']='1'
        launch(['ctest','--test-dir',plan['build'],'-L','required','--output-on-failure','-j','1'],execution,single,300,'required.log')
        for ranks in (1,4):
            for mode in ('manufactured','injected','pulse'):
                name=f'{mode}-{ranks}rank';out=execution/name
                launch(srun(ranks)+[plan['operator_checks'],str(out),mode],execution,env,90,name+'.log')
                require(json.loads((out/'summary.json').read_text())['pass'],'Operator probe failed')
        for name in ('off','on'):
            case=Path(plan['cases'][name]);launch(srun()+[plan['binary'],str(case/'run.yaml')],case,env,600)
            item={'steps':257,'ledger':True,'energy_planes':False,'binary':'new','magnet':'free-drift'}
            save(case/'CASE_REVIEW.json',review.inspect(case,item,plan,export_csv=False))
        comparison=review.compare(Path(plan['cases']['off']),Path(plan['cases']['on']));require(comparison['checks_pass'],'Flag changed physics')
        save(root/'PHYSICS_BITWISE_REVIEW.json',comparison)
        old=review.compare(Path(plan['baseline']),Path(plan['cases']['off']));require(old['checks_pass'],'Disabled candidate changed accepted physics')
        save(root/'BASELINE_BITWISE_REVIEW.json',old)
        result=review_discrete.inspect(root,plan);require(result['checks_pass'],'Discrete field ledger review failed')
        require(all(sha(p)==v for p,v in frozen.items()),'Frozen inputs changed during computation')
        save(root/'COMPLETION.json',{'checks_pass':True,'job_id':os.environ['SLURM_JOB_ID'],'binary_sha256':plan['binary_sha256'],
             'scope':'Independent discrete-operator identity; numerical exchange channels remain explicit; not a general physical-conservation acceptance'})
    except BaseException as error:save(root/'FAILED_OR_INCOMPLETE.json',{'checks_pass':False,'error':str(error)});raise
if __name__=='__main__':main()
