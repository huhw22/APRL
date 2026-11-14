"""Generate both independent sources and gate all production on closed checks."""
import os
import re
import subprocess
import sys
import time
import traceback
from pathlib import Path
import h5py
import numpy as np
from common import W,read,save,sha,need,frozen,srun,launch,dump_at
from binding import binding


def check_import(root,plan):
    launch(srun(1,128)+[plan['binary'],'main.in'],root/'import-check',root/'preparation/import.log',1200,True)
    need('Puffin simulation completed' in (root/'preparation/import.log').read_text(),'Native check did not close')
    with h5py.File(root/'particles.h5') as f:a=f['electrons'][:]
    with h5py.File(root/'import-check/main_electrons_0.h5') as f:b=f['electrons'][:]
    need(a.shape==b.shape and np.isfinite(b).all(),'Native population/nonfinite changed')
    ia=np.argsort(a[:,0]);ib=np.argsort(b[:,0])
    need(np.all(np.diff(a[ia,0])>0),'Repeated native transverse key')
    need(np.array_equal(a[ia],b[ib]),'Actual native import changed one of seven columns')
    first=dump_at(root/'import-check','electrons',1)
    with h5py.File(first) as f:v=f['electrons'][:]
    weighted=lambda x:np.sum(x[:,3].astype(np.longdouble)*x[:,6].astype(np.longdouble))/np.sum(x[:,6],dtype=np.longdouble)
    delta=float((weighted(v)-weighted(a))*plan['K'])
    expected=-plan['K']*.5*(1-np.cos(np.pi/80))*np.cos(2*np.pi/80)
    need(abs(delta-expected)<1e-5,'Actual first-step magnetic trajectory misses smooth entrance')
    result=dict(passed=True,records=len(a),all_seven_columns_bitwise=True,
        input_native_sha256=sha(root/'particles.h5'),matching='Exact rows matched by unique transverse x',
        observed_mean_delta_ux=delta,on_axis_expected_delta_ux=expected,
        smooth_magnetic_entrance_verified=True,feedback_disabled_for_check=True)
    save(root/'NATIVE_VALIDATED.json',result)
    save(root/'PREPARATION_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],
        actual_native_import_bitwise=True,independent_random_source=True,
        actual_weight=read(root/'INPUT_VALIDATED.json')['actual_weight'],field_route_and_stock_path_reused=True,
        real_first_step_magnetic_route_passed=True,physical_acceptance=False))


def pic_preflight(root,plan):
    env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
    env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores',
        OMPI_MCA_btl='self,vader,tcp',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',HDF5_USE_FILE_LOCKING='FALSE')
    cmd=['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1','--ntasks=16','--ntasks-per-node=16',
        '--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    with (root/'binding.log').open('x') as f:
        subprocess.run(cmd+[str(W/'build/cluster-tools/binding_probe'),'16','8'],env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
    save(root/'BINDING.json',binding(root/'binding.log',1,16,8))
    results=[]
    for case in plan['preflight']:
        d=Path(case['directory']);e=env.copy()
        if case['observer']:
            e.update(APRL_LAB_WORLDTUBE=str(d/'worldtube'),APRL_LAB_FIELD_AUDIT=str(d/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(d/'lab-audit'))
        started=time.monotonic()
        with (d/'aprl.log').open('x') as f:
            proc=subprocess.run(['timeout','--signal=TERM','--kill-after=30','1800']+cmd+[plan['binary'],str(d/'run.yaml')],
                cwd=d,env=e,stdout=f,stderr=subprocess.STDOUT,timeout=1870)
        log=(d/'aprl.log').read_text()
        need(proc.returncode not in [0,124,137],'Preflight did not reach its expected intentional duration cutoff')
        for marker in ['mesh.duration was exhausted before the configured physical stop',
            'phase=duration-exhaustion-readiness-sync-complete','Moving window completed: steps=40']:
            need(marker in log,'Missing preflight marker: '+marker)
        need(not re.search('Segmentation fault|OpenFabrics device',log),'Preflight runtime failure')
        with h5py.File(d/'ledger/energy-ledger.h5') as f:
            need(int(f['energy_ledger/complete'][()])==0,'Short check must remain physically incomplete')
            need(int(f['energy_ledger/records'][-1]['step'])==40,'Incomplete short check')
        if case['observer']:
            need((d/'worldtube/COMPLETE').read_text().split()==['0','40'],'Worldtube did not close at deliberate cutoff')
            for rank in range(16):
                for prefix in ['lab-field','paired-momentum','whole-cowan']:
                    need((d/'lab-audit'/f'{prefix}-rank-{rank:05d}.complete').read_text().split()==['0','40'],'Missing rank ledger')
            data=np.fromfile(d/'worldtube/field_lines.bin',dtype='<f8')
            need(len(data)==40*plan['active_nz']*6 and np.isfinite(data).all(),'Invalid raw worldtube')
            need(len(list((d/'worldtube').glob('particle-sides-rank-*.csv')))==16,'Missing side particle ledgers')
        results.append(dict(case=case['name'],diagnostic_preflight_passed=True,intentional_duration_cutoff=True,
            physical_run_complete=False,seconds=time.monotonic()-started))
        print('PIC_PREFLIGHT',results[-1],flush=True)
    a=root/'preflight/binomial3-off';b=root/'preflight/binomial3-on';same=[]
    for path in sorted((a/'detectors').glob('*.h5')):
        with h5py.File(path) as f,h5py.File(b/'detectors'/path.name) as g:
            fg=f[next(iter(f))];gg=g[next(iter(g))]
            for key in ['records','fields','time_s']:
                if key in fg:need(np.array_equal(fg[key][:],gg[key][:]),'Observer changed '+path.name+' '+key)
        same.append(path.name)
    with h5py.File(a/'ledger/energy-ledger.h5') as f,h5py.File(b/'ledger/energy-ledger.h5') as g:
        x=f['energy_ledger/records'][:];y=g['energy_ledger/records'][:]
        for key in x.dtype.names:
            if not key.endswith('_rank_seconds'):need(np.array_equal(x[key],y[key]),'Observer changed scientific ledger '+key)
    save(root/'PREFLIGHT_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],binary=plan['binary'],
        binary_sha256=sha(plan['binary']),cases=results,off_on_scientific_h5_bitwise=same,
        same_full_input=True,input_sha256=plan['input_sha256'],independent_random_source=True,
        physics_core_unchanged=True,physical_acceptance=False))


def main(root):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm')
    root=Path(root).resolve();frozen(root/'FROZEN_PREPARATION.json')
    campaign=read(root/'PLAN.json')
    launch(srun(1,128)+[str(W/'build/cluster-tools/binding_probe'),'128','1'],root,root/'preparation-binding.log',120)
    save(root/'PREPARATION_BINDING.json',binding(root/'preparation-binding.log',1,128,1))
    from make_inputs import generate
    for name in ['six','eleven']:
        d=root/name;plan=read(d/'PLAN.json')
        print('GENERATING',name,plan['sample_records'],flush=True)
        generate(d)
        check_import(d,plan)
        info=read(d/'INPUT_VALIDATED.json')
        need(info['zero_transverse_duplicates'] and info['independent_fixed_N_source'],'Wrong source loading')
        save(d/'FROZEN_PHYSICS.json',{str(p):sha(p) for p in [d/'PLAN.json',d/'particles.h5',Path(info['input']),
            d/'INPUT_VALIDATED.json',d/'NATIVE_VALIDATED.json',d/'PREPARATION_COMPLETE.json',Path(plan['binary']),
            *[p for p in (root/'scripts').iterdir() if p.is_file()],*[p for p in (d/'main').iterdir() if p.is_file()]]})
    from pic_cards import prepare
    plan=prepare(root)
    need(plan['input_sha256']==read(root/'six/INPUT_VALIDATED.json')['input_sha256'],'Six input identity differs')
    pic_preflight(root/'six-aprl',plan)
    files=[root/'six-aprl/PLAN.json',root/'six-aprl/INPUT_VALIDATED.json',root/'six-aprl/GEOMETRY.json',
        root/'six-aprl/PREFLIGHT_COMPLETE.json',Path(plan['binary']),Path(plan['input_file']),
        *[p for p in (root/'scripts').iterdir() if p.is_file()],
        *[Path(x['directory'])/'run.yaml' for x in plan['cases']+plan['preflight']]]
    save(root/'six-aprl/FROZEN.json',{str(p):sha(p) for p in files})
    frozen(root/'FROZEN_PREPARATION.json')
    save(root/'PREPARATION_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],
        independent_inputs_validated=True,six_same_input=True,six_same_magnetic_lattice=True,
        puffin_actual_native_imports_and_entrances_passed=True,aprl_full_input_and_observer_checks_passed=True,
        aprl_steps=plan['steps'],aprl_grid=[56,56,plan['active_nz']],physical_acceptance=False))
    print('PREPARATION_COMPLETE',flush=True)


if __name__=='__main__':
    try:main(Path(sys.argv[1]))
    except Exception:
        save(Path(sys.argv[1])/'PREPARATION_FAILED.json',dict(job=os.environ.get('SLURM_JOB_ID'),error=traceback.format_exc()))
        raise
