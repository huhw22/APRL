"""Same physical events and binary, explicit validated undulator route."""
import os
import hashlib
import subprocess
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c
import analyze
from analyze import read,save,need,PUF as PREVIOUS,PIC


def launch(root,binary,args,ranks,limit,log,cwd):
    command=['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1',f'--ntasks={ranks}',
        f'--ntasks-per-node={ranks}','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30',str(binary),*args]
    with Path(log).open('x') as stream:
        status=subprocess.run(['timeout','--signal=TERM','--kill-after=30',str(limit),*command],cwd=cwd,
                              stdout=stream,stderr=subprocess.STDOUT,timeout=limit+45)
    need(status.returncode==0,'Program failed: '+str(log))


def main(root):
    need(bool(os.environ.get('SLURM_JOB_ID')),'Bounded allocation required');root=Path(root).resolve()
    plan=read(root/'PLAN.json');d=root/'short-u07'
    for path,expected in read(root/'FROZEN.json').items():
        h=hashlib.sha256()
        with Path(path).open('rb') as stream:
            for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
        need(h.hexdigest()==expected,'Frozen control identity changed: '+path)
    launch(root,plan['route_probe'],[str(d/'matched.latt')],1,60,root/'route.log',root)
    need('ROUTE_PASS' in (root/'route.log').read_text(),'Actual setupMods/initUndulator route not verified')
    save(root/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],route_verified=True,
        actual_undulator_type='aprl-ramp77',automatic_match_in_out=False,new_seeds=0))
    launch(root,plan['binary'],['main.in'],128,900,d/'puffin.log',d)
    need('Puffin simulation completed' in (d/'puffin.log').read_text(),'Physics did not close')
    save(root/'TRACKING_CLOSED.json',dict(complete=True,job=os.environ['SLURM_JOB_ID']))
    analyze.PUF=root
    first=analyze.puffin_dump('electrons',1);initial=analyze.puffin_dump('electrons',0)
    with h5py.File(root/'particles.h5','r') as f,h5py.File(initial,'r') as g:
        a=f['electrons'][:];b=g['electrons'][:]
        need(a.shape==b.shape and np.array_equal(a[np.argsort(a[:,0])],b[np.argsort(b[:,0])]),'Native import changed')
    with h5py.File(first,'r') as f:
        b=f['electrons'][:];distance=float(f['runInfo'].attrs['zbarTotal'])*float(f['runInfo'].attrs['Lg'])
    K=read(root/'INPUT_VALIDATED.json')['scale']['aw_peak']
    observed=float((b[:,3].mean()-a[:,3].mean())*K)
    expected=-K*.5*(1-np.cos(np.pi*distance/.028))*np.cos(2*np.pi*distance/.028)
    need(abs(observed-expected)<1e-4,'First actual physical step does not follow smooth ramp')
    save(root/'PHYSICAL_ROUTE_VALIDATED.json',dict(passed=True,first_actual_snapshot=str(first),
        physical_step_m=distance,mean_delta_ux=observed,on_axis_expected_delta_ux=float(expected),
        native_all_seven_columns_bitwise=True,automatic_synthetic_match_disabled=True,
        actual_segment_type_validated=True))
    with np.load(analyze.PREV/'analysis/initial_sources.npz') as f:frequency=f['frequency_Hz'][::2]
    out=root/'analysis';out.mkdir()
    analyze.measure(out,frequency);analyze.source_maps(out,frequency[::2]);analyze.propagation(out)
    import phase
    phase.main(root)
    save(root/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],new_seeds=0,new_PIC_runs=0,
        corrected_route=True,physical_acceptance=False))


if __name__=='__main__':main(sys.argv[1])
