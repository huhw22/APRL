#!/usr/bin/env python3
import os,json,sys
from pathlib import Path
import h5py,numpy as np
from prepare import sha,save,require
from launcher import launch,binding
from analyze import analyze

def main():
    root=Path(sys.argv[1]);name=sys.argv[2];plan=json.loads((root/'plan.json').read_text())
    frozen=json.loads((root/'FROZEN.json').read_text());assert all(sha(p)==h for p,h in frozen.items())
    assert os.environ.get('SLURM_JOB_ID') and os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
    d=root/name;d.mkdir(exist_ok=True);env=os.environ.copy()
    env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
    env.pop('APRL_PARTICLE_MOMENTUM_PROBE',None);env.pop('APRL_PARTICLE_MOMENTUM_PROBE_DT',None)
    sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    try:
        launch(sr+[plan['binding_probe'],'16','8'],d,env,45,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
        if name=='preflight':
            if plan.get('required_source'):
                required=Path(plan['required_source'])/'preflight/required.log'
                assert '100% tests passed, 0 tests failed out of 77' in required.read_text()
                save(d/'REUSED_REQUIRED.json',dict(path=str(required),sha256=sha(required)))
            else:
                testenv=env.copy();testenv['OMP_NUM_THREADS']='1'
                launch(['ctest','--test-dir',plan['build'],'-L','required','--output-on-failure','-j','1'],d,testenv,300,'required.log')
            for case in ('probe_off','probe_on'):
                item=plan['cases'][case];p=Path(item['directory']);ev=env.copy()
                if item['probe']:
                    ev['APRL_PARTICLE_MOMENTUM_PROBE']=str(p/'momentum');ev['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s'])
                launch(sr+[item['binary'],str(p/'run.yaml')],p,ev,item['timeout_s'])
            comparisons={}
            for sub in ('ledger','trajectory','detectors'):
                a=root/'probe_off'/sub;b=root/'probe_on'/sub
                for p in a.glob('*.h5'):
                    with h5py.File(p) as h,h5py.File(b/p.name) as g:
                        group=next(iter(h));x=h[group+'/records'][:];y=g[group+'/records'][:]
                        assert x.shape==y.shape and x.dtype==y.dtype
                        fields=[n for n in x.dtype.names if 'seconds' not in n]
                        assert all(np.ascontiguousarray(x[n]).tobytes()==np.ascontiguousarray(y[n]).tobytes() for n in fields)
                        comparisons[str(p.relative_to(root/'probe_off'))]=len(x)
            save(d/'BITWISE.json',dict(pass_all=True,files=comparisons))
            assert sum(p.stat().st_size for p in (root/'probe_on/momentum').glob('*.bin'))>0
        else:
            assert json.loads((root/'preflight/COMPLETION.json').read_text())['pass_all']
            item=plan['cases'][name];p=Path(item['directory'])
            ev=env.copy();ev['APRL_PARTICLE_MOMENTUM_PROBE']=str(p/'momentum');ev['APRL_PARTICLE_MOMENTUM_PROBE_DT']=str(plan['nominal_dt_s'])
            launch(sr+[item['binary'],str(p/'run.yaml')],p,ev,item['timeout_s'])
            summary=analyze(p,Path(item['original']),p/'analysis')
            assert summary['steps']==item['expected_steps'] and summary['first_exit_step'] is None
            if name=='base':
                with h5py.File(Path(item['original'])/'ledger/energy-ledger.h5') as h,h5py.File(p/'ledger/energy-ledger.h5') as g:
                    a=h['energy_ledger/records'][:501];b=g['energy_ledger/records'][:]
                    fields=[n for n in a.dtype.names if 'seconds' not in n]
                    check={n:np.ascontiguousarray(a[n]).tobytes()==np.ascontiguousarray(b[n]).tobytes() for n in fields}
                    save(p/'ORIGINAL_LEDGER_PREFIX.json',check)
                    # End-stop metadata may differ; require scientific momentum and energy exactly.
                    for n in ('particle_cPz_box_active_J','particle_kinetic_active_J','particle_grid_work_J'):
                        assert check[n],n
        assert all(sha(p)==h for p,h in frozen.items())
        save(d/'COMPLETION.json',dict(pass_all=True,job_id=os.environ['SLURM_JOB_ID'],case=name))
    except BaseException as error:
        save(d/'FAILED.json',dict(pass_all=False,error=str(error)));raise

if __name__=='__main__':main()
