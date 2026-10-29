#!/usr/bin/env python3
"""Two-node validation gated by the completed single-node science campaign."""
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
    require(os.environ.get('SLURM_JOB_NUM_NODES')=='2' and os.environ.get('SLURM_NTASKS')=='16','Require real 2-node16x8 Slurm')
    gate=json.loads((Path(plan['single_node_root'])/'COMPLETION.json').read_text())
    require(gate['checks_pass'] and gate['binary_sha256']==plan['binary_sha256'],'Same-binary single-node science gate failed')
    require(all(sha(p)==v for p,v in frozen.items()),'MPI inputs changed')
    exe=root/'execution';exe.mkdir();env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    for n in ('APRL_GRID_WORK_AUDIT','APRL_FIELD_STAGE_AUDIT'):env.pop(n,None)
    def srun(ranks=16):return ['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=2','--ntasks='+str(ranks),'--ntasks-per-node='+str(ranks//2),'--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    try:
        launch(srun()+[plan['binding_probe'],'16','8'],exe,env,45,'binding.log');save(exe/'BINDING.json',binding(exe/'binding.log',2))
        for mode in ('manufactured','injected','pulse'):
            directory=exe/mode
            launch(srun(4)+[plan['operator_checks'],str(directory),mode],exe,env,90,mode+'.log')
            require(json.loads((directory/'summary.json').read_text())['pass'],'MPI operator module failed')
        for name in ('off','on'):
            case=Path(plan['cases'][name]);launch(srun()+[plan['binary'],str(case/'run.yaml')],case,env,600)
            item={'steps':257,'ledger':True,'energy_planes':False,'binary':'new','magnet':'free-drift'}
            save(case/'CASE_REVIEW.json',review.inspect(case,item,plan,export_csv=False))
        pair=review.compare(Path(plan['cases']['off']),Path(plan['cases']['on']));require(pair['checks_pass'],'MPI flag changed scientific output')
        save(root/'PHYSICS_BITWISE_REVIEW.json',pair)
        require(review_discrete.inspect(root,plan)['checks_pass'],'MPI discrete review failed')
        require(all(sha(p)==v for p,v in frozen.items()),'MPI frozen inputs changed')
        save(root/'COMPLETION.json',{'checks_pass':True,'job_id':os.environ['SLURM_JOB_ID'],'binary_sha256':plan['binary_sha256'],'nodes':2})
    except BaseException as e:save(root/'FAILED_OR_INCOMPLETE.json',{'checks_pass':False,'error':str(e)});raise
if __name__=='__main__':main()
