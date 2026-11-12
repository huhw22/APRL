"""Bounded fixture/import checks followed by one same-sample seven-period run."""
import os
import subprocess
import sys
import time
from pathlib import Path

from prepare import read,save,sha,W


def run(root,directory,binary,arguments,ranks,limit,log,env=None):
    command=['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1',f'--ntasks={ranks}',
             f'--ntasks-per-node={ranks}','--cpus-per-task=1','--cpu-bind=cores',
             '--kill-on-bad-exit=1','--wait=30',str(binary),*arguments]
    with Path(log).open('x') as stream:
        completed=subprocess.run(['timeout','--signal=TERM','--kill-after=30',str(limit),*command],
            cwd=directory,stdout=stream,stderr=subprocess.STDOUT,timeout=limit+45,env=env)
    if completed.returncode:raise RuntimeError(f'Run failed: {log}, exit {completed.returncode}')


def validate(stage,root):
    subprocess.run([sys.executable,'-B',str(root/'scripts/validate.py'),stage,str(root)],check=True,timeout=150)


def main(root):
    plan=read(root/'PLAN.json')
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Use a bounded Slurm allocation')
    for path,digest in read(root/'FROZEN.json').items():
        if sha(path)!=digest:raise RuntimeError('Frozen identity changed: '+path)
    if 'validated_native_root' in plan:
        previous=Path(plan['validated_native_root'])
        for name in ['SMOKE_VALIDATED.json','FIELD_VALIDATED.json','IO_VALIDATED.json','NATIVE_VALIDATED.json']:
            if not read(previous/name)['passed']:raise RuntimeError('Required closed check failed: '+name)
        if sha(root/'particles.h5')!=read(previous/'NATIVE_VALIDATED.json')['input_native_sha256']:
            raise RuntimeError('Native particle source differs from actual imported source')
        os.environ['OMPI_MCA_io']='romio321'
    elif 'passed_software_checks_root' in plan:
        previous=Path(plan['passed_software_checks_root'])
        for name in ['SMOKE_VALIDATED.json','FIELD_VALIDATED.json']:
            if not read(previous/name)['passed']:raise RuntimeError('Prior software check did not pass')
        for component,limit in [('ompio',20),('romio321',60)]:
            env=os.environ.copy();env['OMPI_MCA_io']=component
            directory=root/('io-'+component)
            try:
                run(root,directory,plan['io_probe_binary'],[],8,limit,directory/'probe.log',env)
            except RuntimeError:
                if component=='romio321':raise
        validate('io',root)
        os.environ['OMPI_MCA_io']='romio321'
    else:
        for name,binary in [('stock',plan['stock_binary']),('matched',plan['matched_binary'])]:
            directory=root/('smoke-'+name)
            run(root,directory,binary,['inputs/3D/clara_test.in'],2,180,root/(name+'-smoke.log'))
        validate('smoke',root)
        run(root,root,plan['probe_binary'],[],1,60,root/'field_probe.dat')
        validate('probe',root)
    if 'validated_native_root' not in plan:
        run(root,root/'import-check',plan['matched_binary'],['main.in'],plan.get('native_import_ranks',8),240,root/'import-check/puffin.log')
        validate('native',root)
    run(root,root,W/'build/cluster-tools/binding_probe',['128','1'],128,90,root/'binding.log')
    save(root/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],node=os.environ.get('SLURMD_NODENAME'),
        time_unix=time.time(),ranks=128,threads_per_rank=1,case='w64-a',new_seeds=0,
        field_validated=True,actual_native_import_bitwise=True,stock_path_validated=True))
    run(root,root/'short-u07',plan['matched_binary'],['main.in'],128,2700,root/'short-u07/puffin.log')
    save(root/'TRACKING_CLOSED.json',dict(job=os.environ['SLURM_JOB_ID'],complete=True))
    subprocess.run([sys.executable,'-B',str(root/'scripts/analyze.py'),str(root)],check=True,timeout=240)
    save(root/'COMPLETE.json',dict(job=os.environ['SLURM_JOB_ID'],complete=True,new_seeds=0,new_PIC_runs=0,
        full_77_period_simulation=False,physical_acceptance=False))


if __name__=='__main__':main(Path(sys.argv[1]).resolve())
