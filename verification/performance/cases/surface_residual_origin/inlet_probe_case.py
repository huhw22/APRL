#!/usr/bin/env python3
import argparse,copy,json,os,shutil,sys
from pathlib import Path
import h5py,numpy as np,yaml
from prepare import sha,save,require
from launcher import launch,binding
from model import DT,G_BOOST,c

def prepare(root,comparison,build):
    root.mkdir(parents=True,exist_ok=False);(root/'scripts').mkdir();here=Path(__file__).resolve().parent
    for name in ('inlet_probe_case.py','run_probe.sbatch','prepare.py'):shutil.copy2(here/name,root/'scripts'/name)
    for name in ('model.py','launcher.py'):shutil.copy2(here.parent/'space_charge_physics'/name,root/'scripts'/name)
    prior=json.loads((comparison/'plan.json').read_text());original=Path(prior['original']);meta=json.loads((original/'case.json').read_text())
    d=root/'probe_300';d.mkdir();(d/'diagnostic').mkdir();cfg=copy.deepcopy(yaml.safe_load((comparison/'base/run.yaml').read_text()))
    cfg['output']['manifest']=str(d/'run-manifest.yaml');cfg['mesh']['duration']=float(302*DT)
    beta=float(np.sqrt(1-np.longdouble(G_BOOST)**-2));cfg['stop']['z']=float(meta['geometry']['reference_initial_lab_z_m']+G_BOOST*beta*c*299.75*DT)
    cfg['mesh']['moving_window']['reference_directory']=str(d/'reference')
    for key in ('detectors','trajectory','energy_ledger'):cfg[key]['directory']=str(d/('ledger' if key=='energy_ledger' else key))
    (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    plan=dict(root=str(root),comparison=str(comparison),build=str(build),binary=str(build/'aprl'),binding_probe=prior['binding_probe'])
    save(root/'plan.json',plan)
    files=[root/'plan.json',d/'run.yaml',build/'aprl',*list((root/'scripts').iterdir()),Path(cfg['beam']['input']['file'])]
    save(root/'FROZEN.json',{str(p):sha(p) for p in files})

def run(root):
    plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text());d=root/'probe_300'
    assert os.environ['SLURM_NTASKS']=='16' and os.environ['SLURM_CPUS_PER_TASK']=='8'
    try:
        assert all(sha(p)==h for p,h in frozen.items());env=os.environ.copy()
        for k in ('APRL_PARTICLE_MOMENTUM_PROBE','APRL_PARTICLE_MOMENTUM_PROBE_DT','APRL_INLET_CONSTRAINT_PROBE'):env.pop(k,None)
        env['LD_LIBRARY_PATH']=plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH','')
        sr=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
        launch(sr+[plan['binding_probe'],'16','8'],d,env,60,'binding.log');save(d/'BINDING.json',binding(d/'binding.log',1))
        te=env.copy();te['OMP_NUM_THREADS']='1';launch(['ctest','--test-dir',plan['build'],'-L','required','--output-on-failure','-j','1'],d,te,300,'required.log')
        ev=env.copy();ev['APRL_INLET_CONSTRAINT_PROBE']=str(d/'diagnostic');launch(sr+[plan['binary'],str(d/'run.yaml')],d,ev,500)
        checked={}
        for sub in ('ledger','trajectory','detectors'):
            for p in (d/sub).glob('*.h5'):
                with h5py.File(p) as h,h5py.File(Path(plan['comparison'])/'base'/sub/p.name) as old:
                    g=next(iter(h));assert int(h[g+'/complete'][()])==1
                    a=h[g+'/records'][:];b=old[g+'/records'][:len(a)];assert a.dtype==b.dtype and a.shape==b.shape
                    assert all(np.ascontiguousarray(a[n]).tobytes()==np.ascontiguousarray(b[n]).tobytes() for n in a.dtype.names if 'seconds' not in n)
                    checked[sub+'/'+p.name]=len(a)
        a=np.genfromtxt(d/'diagnostic/inlet-constraint.csv',names=True,delimiter=',',dtype=None,encoding=None)
        assert len(a)==4800 and all(np.isfinite(a[n]).all() for n in a.dtype.names if n!='stage')
        assert all(sha(p)==h for p,h in frozen.items())
        save(d/'COMPLETION.json',dict(checks_pass=True,job_id=os.environ['SLURM_JOB_ID'],science_prefix_bitwise=checked,read_only_probe_rows=len(a)))
    except BaseException as error:
        save(d/'FAILED.json',dict(error=str(error)));raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('root',type=Path);p.add_argument('--comparison',type=Path);p.add_argument('--build',type=Path);a=p.parse_args()
    if a.action=='prepare':prepare(a.root.resolve(),a.comparison.resolve(),a.build.resolve())
    else:run(a.root.resolve())
