import os,json
from pathlib import Path
from prepare import require,sha,save
from launcher import launch,binding
R=Path(__file__).resolve().parents[1]
def main():
    for f,h in json.loads((R/'FROZEN.json').read_text()).items():require(sha(f)==h,'Frozen input changed '+f)
    p=json.loads((R/'PLAN.json').read_text());env=os.environ.copy();srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];launch(srun+[p['binding_probe'],'16','8'],R,env,45,'binding.log');save(R/'BINDING.json',binding(R/'binding.log',1));small=env.copy();small['OMP_NUM_THREADS']='1';reports=[]
    for ranks in [1,2,4]:
        print('START checks',ranks,flush=True);cmd=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks='+str(ranks),'--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];name='checks-'+str(ranks)+'.log';launch(cmd+[p['programs']['checks']],R,small,120,name);reports.append(json.loads((R/name).read_text()))
    cmd=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=1','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];periodic=[]
    for wave,quiet in [('axial',2),('oblique',1),('oblique',2),('oblique',4)]:
        label=wave+'-'+str(quiet);print('START coupled periodic particle test',label,flush=True);launch(cmd+[p['programs']['periodic'],str(R/('periodic-'+label+'.csv')),str(quiet),wave],R,small,120,'periodic-'+label+'.log');item=json.loads((R/('periodic-'+label+'.log')).read_text());item.update(wave=wave,quiet_refinement=quiet);periodic.append(item)
    save(R/'RESULT.json',dict(units=reports,periodic=periodic,production_replaced=False,physical_beam_accuracy_accepted=False));print('COMPLETE',flush=True)
if __name__=='__main__':main()
