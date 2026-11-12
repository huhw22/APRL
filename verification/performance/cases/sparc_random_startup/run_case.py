"""Bounded first-segment PIC and matched native GENESIS controls."""
import os
import re
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

import h5py
import numpy as np

from common import W, frozen, load, require, save
from binding import binding


def srun(ranks, threads):
    return ['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1',
        f'--ntasks={ranks}',f'--ntasks-per-node={ranks}',f'--cpus-per-task={threads}',
        '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']


def launch(command,directory,env,limit,filename):
    command=['timeout','--signal=TERM','--kill-after=20',str(limit),*command]
    record=dict(command=command,job=os.environ['SLURM_JOB_ID'],start_unix=time.time())
    start=time.monotonic();process=None
    try:
        with (directory/filename).open('x') as stream:
            process=subprocess.Popen(command,cwd=directory,env=env,stdin=subprocess.DEVNULL,
                stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            record['exit_code']=process.wait(timeout=limit+40)
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=20)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
        record['seconds']=time.monotonic()-start
        save(directory/(filename+'.launch.json'),record)
    require(record.get('exit_code')==0,'Bounded process failed: '+str(directory/filename))
    return record


def verify_import(directory,export):
    mapping=load(directory/'NATIVE_MAPPING.json')
    with h5py.File(directory/'common.par.h5','r') as a,h5py.File(export,'r') as b:
        for i in mapping['active_groups']:
            key=f'slice{i:06d}'
            for field in ('current','x','y','px','py','gamma','theta'):
                require(np.array_equal(a[key][field][:],b[key][field][:]),
                        f'GENESIS changed imported {key}/{field}')
    return True


def main(root,index):
    root=Path(root);plan=load(root/'PLAN.json');entry=plan['cases'][int(index)]
    d=Path(entry['directory'])
    require(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm')
    frozen(root/'FROZEN_PREPARATION.json');frozen(d/'FROZEN.json')
    require(load(root/'INPUTS_COMPLETE.json')['complete'],'Preparation incomplete')
    require(load(d/'IMPORT_VERIFIED.json')['passed'],'Native import not verified')
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('APRL_'):env.pop(key)
    env['APRL_PROFILE_INTERVAL']='0'
    launch(srun(16,8)+[str(W/'build/cluster-tools/binding_probe'),'16','8'],d,env,120,'binding.log')
    save(d/'BINDING.json',binding(d/'binding.log',1,16,8))
    save(d/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],start_unix=time.time(),case=entry['name']))
    runtime=launch(srun(16,8)+[plan['pic_binary'],str(d/'run.yaml')],d,env,entry['timeout_seconds'],'aprl.log')
    text=(d/'aprl.log').read_text()
    require(all(x in text for x in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),
            'PIC did not normally complete')
    require(not re.search('Simulation failure|OpenFabrics device|Segmentation fault',text),'PIC error')
    from analyze_closed import verify
    closed=verify(d)
    save(d/'PIC_COMPLETE.json',dict(complete=True,runtime=runtime,closed=closed,physical_acceptance=False))
    for mode in ('off','longrange'):
        gd=d/('genesis-'+mode);ge=env.copy();ge['OMP_NUM_THREADS']='1'
        launch(srun(32,1)+[str(W/'build/cluster-tools/binding_probe'),'32','1'],gd,ge,120,'binding.log')
        save(gd/'BINDING.json',binding(gd/'binding.log',1,32,1))
        launch(srun(32,1)+[plan['genesis_binary'],'reference.in'],gd,ge,1800,'genesis.log')
        require('End of Track' in (gd/'genesis.log').read_text(),'GENESIS did not complete')
        verify_import(d,gd/'initial.par.h5')
        with h5py.File(gd/'reference.out.h5','r') as f:
            for key in ('Field/power','Beam/bunching','Beam/bunchingphase','Beam/energy'):
                require(np.isfinite(f[key][:]).all(),'Nonfinite GENESIS output: '+key)
            require(f['Lattice/zplot'][-1]>=2.166-1e-8,'GENESIS ended early')
        save(gd/'COMPLETE.json',dict(complete=True,native_import_bitwise=True))
    frozen(d/'FROZEN.json');frozen(root/'FROZEN_PREPARATION.json')
    from analyze import analyze
    analyze(root,entry)
    save(d/'COMPLETE.json',dict(complete=True,analysis_complete=True,physical_acceptance=False,
        job=os.environ['SLURM_JOB_ID']))


if __name__=='__main__':
    try: main(sys.argv[1],sys.argv[2])
    except Exception:
        d=Path(sys.argv[1]);save(d/f'FAILED-{sys.argv[2]}.json',dict(error=traceback.format_exc()))
        raise
