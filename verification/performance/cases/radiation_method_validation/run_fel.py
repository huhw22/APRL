import json
import os
import traceback
from pathlib import Path
import h5py
import numpy as np
from launcher import launch,binding
from prepare import save,require,sha

def genesis_summary(directory,expected_length):
    with h5py.File(directory/'reference.out.h5','r') as f:
        z=f['Lattice/zplot'][:]
        lattice=f['Lattice/z'][:]
        dz=f['Lattice/dz'][:]
        require(abs(lattice[-1]+dz[-1]-expected_length)<1e-9,'Incomplete GENESIS lattice')
        p=f['Field/power'][:]
        b=f['Beam/bunching'][:]
        require(np.isfinite(p).all() and np.isfinite(b).all(),'Nonfinite GENESIS')
        current=f['Beam/current'][:].reshape(-1,p.shape[-1])[0]
        require(np.any(current>0),'Reference has no physical beam')
        return dict(final_power_sum_W=float(p[-1].sum()),max_bunching=float(b[:,current>0].max()),
                    power_shape=list(p.shape),z_power_sample_m=float(z[-1]),
                    z_lattice_end_m=float(lattice[-1]+dz[-1]))

def main():
    root=Path(__file__).resolve().parents[1]
    plan=json.loads((root/'campaign.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text())
    require(all(sha(p)==h for p,h in frozen.items()),'Frozen files changed')
    require(os.environ.get('SLURM_NTASKS')=='16','Expected 16 MPI ranks')
    env=os.environ.copy()
    srun=['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1','--ntasks=16',
          '--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    launch(srun+[plan['binding_probe'],'16','8'],root,env,45,'binding.log')
    save(root/'BINDING.json',binding(root/'binding.log',1))
    refs={};genenv=env.copy();genenv['OMP_NUM_THREADS']='1'
    for case in plan['genesis_references']:
        d=Path(case['directory']);print('GENESIS START',case['name'],flush=True)
        try:
            if case.get('reuse_closed'):
                rec=json.loads((d/'genesis.log.launch.json').read_text())
                require(rec['exit_code']==0,'Reused reference is incomplete')
            else:
                rec=launch(srun+[plan['genesis_binary'],'reference.in'],d,genenv,case['timeout_seconds'],'genesis.log')
            refs[case['name']]={**genesis_summary(d,plan['parameters']['periods']*plan['parameters']['period_m']+.3),**rec}
        except Exception as exc:
            refs[case['name']]={'complete':False,'error':str(exc)}
            traceback.print_exc()
        save(root/'GENESIS_STATUS.json',refs)
    require(all('final_power_sum_W' in refs.get(n,{}) for n in ['strong','weak']),'No valid independent FEL reference')
    ratio=refs['strong']['final_power_sum_W']/refs['weak']['final_power_sum_W']/plan['parameters']['peak_current_A']**2
    gate=dict(strong_over_weak_current_squared=ratio,
              max_bunching=refs['strong']['max_bunching'],
              pass_feedback=bool(ratio>1.2 and refs['strong']['max_bunching']>.06),
              scope='Ensures this test probes feedback; not PIC accuracy acceptance')
    save(root/'REFERENCE_FEEDBACK_GATE.json',gate)
    require(gate['pass_feedback'],'Chosen finite pulse does not sufficiently probe FEL feedback')
    results=[]
    for case in plan['cases']:
        d=Path(case['directory']);print('PIC START',case['name'],flush=True)
        try:
            rec=launch(srun+[plan['binary'],str(d/'run.yaml')],d,env,case['timeout_seconds'])
            save(d/'PIC_COMPLETE.json',rec)
            from analyze_fel import analyze
            result=analyze(root,case)
            results.append({**result,'run':rec})
        except Exception as exc:
            results.append(dict(case=case['name'],complete=False,error=str(exc),physical_acceptance=False))
            save(d/'FAILED_OR_INCOMPLETE.json',results[-1]);traceback.print_exc()
        save(root/'STATUS.json',results)
        require(all(sha(p)==h for p,h in frozen.items()),'Frozen files changed during simulation')
    save(root/'CAMPAIGN_COMPLETE.json',dict(cases=results,physical_acceptance=False))
    require(all(q['complete'] for q in results),'Some FEL cases incomplete')

if __name__=='__main__':main()
