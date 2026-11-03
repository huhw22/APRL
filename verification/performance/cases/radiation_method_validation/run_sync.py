"""Three methods share one binary, one input and one allocation."""
import json
import os
import traceback
from pathlib import Path
from launcher import launch,binding
from prepare import save,require,sha

def main():
    root=Path(__file__).resolve().parents[1]
    plan=json.loads((root/'campaign.json').read_text())
    frozen=json.loads((root/'FROZEN.json').read_text())
    require(all(sha(p)==h for p,h in frozen.items()),'Frozen input changed')
    require(os.environ.get('SLURM_NTASKS')=='16','Expected 16 MPI ranks')
    env=os.environ.copy()
    srun=['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1',
          '--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8',
          '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(srun+[plan['binding_probe'],'16','8'],root,env,45,'binding.log')
    save(root/'BINDING.json',binding(root/'binding.log',1))
    records=[]
    for case in plan['cases']:
        d=Path(case['directory'])
        print('START',case['name'],flush=True)
        try:
            rec=launch(srun+[plan['binary'],str(d/'run.yaml')],d,env,case['timeout_seconds'])
            save(d/'PIC_COMPLETE.json',rec)
            from analyze_dipole import analyze
            analyze(d)
            records.append({'case':case['name'],'PIC_and_analysis_complete':True,'run':rec})
        except Exception as exc:
            save(d/'FAILED_OR_INCOMPLETE.json',{'error':str(exc),'physical_acceptance':False})
            traceback.print_exc()
            records.append({'case':case['name'],'PIC_and_analysis_complete':False,'error':str(exc)})
        save(root/'STATUS.json',records)
        require(all(sha(p)==h for p,h in frozen.items()),'Frozen files changed during run')
    save(root/'CAMPAIGN_COMPLETE.json',{'cases':records,'physical_acceptance':False})
    require(all(q['PIC_and_analysis_complete'] for q in records),'Some cases incomplete')

if __name__=='__main__':main()
