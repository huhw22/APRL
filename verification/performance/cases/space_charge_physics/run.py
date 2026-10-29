#!/usr/bin/env python3
"""Run the frozen physical cases within one real bounded Slurm allocation."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,signal
from pathlib import Path
from prepare import sha,save
from launcher import launch,binding
from analyze import analyze

def main():
    root=Path(sys.argv[1]).resolve();plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'frozen.json').read_text())
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NUM_NODES')=='1'
    assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
    assert all(sha(p)==h for p,h in frozen.items())
    exe=root/'execution';exe.mkdir();env=os.environ.copy();env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    srun=['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    started=time.monotonic();results=[]
    try:
        launch(srun+[plan['binding_probe'],'16','8'],exe,env,45,'binding.log');save(exe/'BINDING.json',binding(exe/'binding.log',1))
        for item in plan['cases']:
            assert time.monotonic()-started+item['timeout_s']+40<plan['active_budget_s'],'Insufficient bounded budget'
            d=Path(item['directory']);print('CASE_START '+item['name']+' '+str(d),flush=True)
            launch(srun+[plan['binary'],str(d/'run.yaml')],d,env,item['timeout_s'])
            print('CASE_SIMULATION_COMPLETE '+item['name'],flush=True)
            result=analyze(d);assert result['execution_integrity_pass']
            results.append(item['name']);print('CASE_ANALYZED '+item['name'],flush=True)
        assert all(sha(p)==h for p,h in frozen.items())
        save(root/'COMPLETION.json',{'execution_integrity_pass':True,'physics_accepted':False,'completed_cases':results,
             'job_id':os.environ['SLURM_JOB_ID'],'scope':'First physical comparisons completed; agreement and convergence require review'})
    except BaseException as error:
        save(root/'FAILED_OR_INCOMPLETE.json',{'execution_integrity_pass':False,'error':str(error),'completed_cases':results});raise
if __name__=='__main__':main()
