"""Prepare one exact-current-crop Genesis first-segment reference; never submit."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

import h5py
import numpy as np
from scipy.constants import c, e


def need(value, message):
    if not value:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


RUNNER = r'''"""Execute only the bounded frozen first-segment Genesis reference."""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib
import h5py,numpy as np

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for q in iter(lambda:f.read(4*1024*1024),b''):h.update(q)
    return h.hexdigest()

root=Path(sys.argv[1]).resolve();plan=json.loads((root/'plan.json').read_text())
frozen=json.loads((root/'FROZEN.json').read_text())
def verify():
    for p,digest in frozen.items():
        if sha(p)!=digest:raise RuntimeError('Frozen reference identity changed: '+p)
verify()
if int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:raise RuntimeError('Use one-node Slurm allocation')
cmd=['timeout','--signal=TERM','--kill-after=15',str(plan['case_limit_seconds']),
     'srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks='+str(plan['ranks']),
     '--ntasks-per-node='+str(plan['ranks']),'--cpus-per-task=1','--cpu-bind=cores',
     '--kill-on-bad-exit=1','--wait=30',plan['binary'],'reference.in']
started=time.time()
with (root/'genesis.log').open('x') as log:
    task=subprocess.run(cmd,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=False)
record=dict(run=dict(command=cmd,exit_code=task.returncode,seconds=time.time()-started,
    started_unix=started),job=os.environ['SLURM_JOB_ID'],new_PIC_runs=0,
    purpose='Exact selected current profile; first original undulator segment only')
(root/'closure_record.json').write_text(json.dumps(record,indent=2)+'\n')
if task.returncode!=0:raise RuntimeError('Genesis exited '+str(task.returncode))
text=(root/'genesis.log').read_text()
if 'Core Simulation done.' not in text or 'End of Track' not in text:raise RuntimeError('Incomplete Genesis log')
with h5py.File(root/'reference.out.h5','r') as f:
    z=f['Lattice/z'][:];dz=f['Lattice/dz'][:]
    endpoint=float(z[-1]+dz[-1]);power=f['Field/power'][:]
    if abs(endpoint-plan['endpoint_m'])>1e-9:raise RuntimeError('Wrong reference endpoint')
    if not np.isfinite(power).all():raise RuntimeError('Nonfinite reference power')
    if len(f['Meta/Fielddumps/intstep'])!=1:raise RuntimeError('Expected one first-segment field dump')
    dump_step=int(f['Meta/Fielddumps/intstep'][0])
for p in [root/('reference.'+str(dump_step)+'.fld.h5'),root/'far.fld.h5']:
    if not p.is_file():raise RuntimeError('Missing completed field dump '+str(p))
verify()
record.update(lattice_endpoint_m=endpoint,dump_step=dump_step,integrity_pass=True,
    physical_acceptance=False,scientific_case='231 selected original slices, no new noise')
(root/'closure_record.json').write_text(json.dumps(record,indent=2)+'\n')
(root/'TRACK_COMPLETE').write_text('0\n')
print(json.dumps(record,indent=2),flush=True)
'''


BATCH = r'''#!/usr/bin/env bash
#SBATCH --job-name=short-genesis-estimate
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:08:00
#SBATCH --no-requeue
set -Eeuo pipefail
reference_root=$(realpath "${1:?Frozen reference directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
ulimit -c 0
scontrol show job "$SLURM_JOB_ID" > "$reference_root/allocation.txt"
"$PYTHON_BIN" "$reference_root/scripts/run_reference.py" "$reference_root" > "$reference_root/driver.log" 2>&1
'''


def prepare(workspace, out):
    workspace, out = Path(workspace).resolve(), Path(out).resolve()
    out.mkdir(exist_ok=False)
    source = workspace/'results/sparc_sase_comparison_20261103/explicit_review_20261104/genesis-phase-r1'
    campaign = workspace/'results/energy_critical_control_20261105/campaign-r4'
    pic = campaign/'inputs/particles.h5'
    selection = json.loads((campaign/'INPUT_SELECTION.json').read_text())
    old_plan = json.loads((workspace/'results/sparc_sase_comparison_20261103/compact-production-r1/plan.json').read_text())
    need((source/'TRACK_COMPLETE').read_text().strip()=='0', 'Preferred reference not closed')
    need(json.loads((source/'PHASE_CLOSURE_VERIFIED.json').read_text())['pass_nominal_closure'], 'Preferred reference phase closure not verified')
    original = source/'common.par.h5'
    original_sha, pic_sha = sha(original), sha(pic)
    new = out/'common.par.h5'
    shutil.copy2(original, new)
    need(original.stat().st_ino != new.stat().st_ino, 'New current profile must be an independent file')
    first, last = selection['first_original_slice'], selection['end_original_slice_exclusive']
    pad, active, nper = 512, old_plan['active_slices'], 1024
    lower, upper = pad+active-(last-1), pad+active-first
    need((first,last,lower,upper)==(577,808,1091,1321), 'Unexpected physical crop mapping')
    with h5py.File(new,'r+') as f:
        total = int(f['slicecount'][0])
        for j in range(1,total+1):
            if not lower<=j<=upper:
                f[f'slice{j:06d}/current'][0] = 0.
    phase_error=0.;gamma_error=0.;charge=np.longdouble(0);groups=0
    with h5py.File(original,'r') as old, h5py.File(new,'r') as f, h5py.File(pic,'r') as p:
        a=p['particles/records'][:]
        need(len(a)==selection['particles'], 'PIC crop particle count differs')
        lam=float(f['slicespacing'][0])
        for key in old:
            if isinstance(old[key],h5py.Dataset):
                need(np.array_equal(old[key][()],f[key][()]), 'Root grid/padding metadata changed: '+key)
        for j in range(1,total+1):
            x,y=old[f'slice{j:06d}'],f[f'slice{j:06d}']
            need(set(x)==set(y), 'Particle dataset names changed')
            for name in ['x','y','px','py','gamma','theta']:
                need(np.array_equal(x[name][:],y[name][:]), 'Original phase space changed')
            current=float(y['current'][0])
            if lower<=j<=upper:
                need(np.array_equal(x['current'][:],y['current'][:]), 'Selected current changed')
                i=pad+active-j
                q=a[(i-first)*nper:(i-first+1)*nper]
                need(np.array_equal(q['source_id'],np.arange(i*nper+1,(i+1)*nper+1)), 'PIC source IDs changed')
                for name,col in [('x',0),('y',1)]:
                    need(np.array_equal(y[name][:],q['plane_position_m'][:,col]), 'PIC transverse positions differ')
                for name,col in [('px',0),('py',1)]:
                    need(np.array_equal(y[name][:],q['proper_velocity'][:,col]), 'PIC transverse momentum differs')
                gamma=np.sqrt(1+(q['proper_velocity']**2).sum(axis=1))
                gamma_error=max(gamma_error,float(np.max(abs(gamma-y['gamma'][:]))))
                phase=np.exp(-2j*np.pi*c*q['arrival_time_offset_s']/lam)
                phase_error=max(phase_error,float(np.max(abs(phase-np.exp(1j*y['theta'][:])))))
                expected_weight=current*lam/(e*c*nper)
                need(np.max(abs(q['macro_weight']/expected_weight-1))<1e-14, 'PIC charge weights differ')
                charge+=np.longdouble(current)*np.longdouble(lam)/np.longdouble(c)
                groups+=1
            else:
                need(current==0., 'Unselected slice still carries charge')
        qpic=np.sum(a['macro_weight'].astype(np.longdouble))*np.longdouble(e)
    need(abs(charge/qpic-1)<1e-14, 'Selected physical charge does not match PIC')
    need(gamma_error<2e-12 and phase_error<1e-10, 'Imported energy or phase mapping differs')
    text=(source/'reference.lat').read_text()
    lines=[x for x in text.splitlines() if x.split(':',1)[0] in ['R0','R1','R2','R3','PREP','U','M']]
    prefix='PREP,R0,R1,R2,R3,U,R3,R2,R1,R0,M'
    need(text.split('FEL: LINE={',1)[1].startswith(prefix+','), 'Original magnetic prefix differs')
    lines.append('FEL: LINE={'+prefix+'};')
    (out/'reference.lat').write_text('\n'.join(lines)+'\n')
    shutil.copy2(source/'reference.in',out/'reference.in')
    card=(out/'reference.in').read_text()
    need('shotnoise=false' in card and 'file=common.par.h5' in card and 'rootname=reference' in card, 'Unexpected Genesis card contract')
    scripts=out/'scripts';scripts.mkdir()
    (scripts/'run_reference.py').write_text(RUNNER)
    (scripts/'run_reference.sbatch').write_text(BATCH)
    shutil.copy2(Path(__file__),scripts/'prepare_genesis_estimate.py')
    proof=dict(passed=True,original_reference=str(source),input_source_sha256=original_sha,
        pic_input_sha256=pic_sha,original_slices=[first,last],
        kept_genesis_groups_inclusive=[lower,upper],selected_slices=groups,
        original_total_slices=total,original_slice_grid_and_padding_retained=True,
        all_phase_space_datasets_bitwise_preserved=True,selected_current_bitwise_preserved=True,
        other_current_exactly_zero=True,new_noise_generated=False,
        selected_charge_C=float(charge),PIC_charge_C=float(qpic),
        charge_relative_difference=float(charge/qpic-1),
        max_gamma_reconstruction_absolute_difference=gamma_error,
        max_optical_phase_unit_complex_difference=phase_error,
        physical_phase_mapping='GENESIS slice index j=512+1386-original_zero_based_slice; theta kept; exp(i theta)=exp(-i2pi c arrival/lambda)',
        lattice_prefix=prefix,first_segment_periods=77,endpoint_m=2.166,
        same_entrance_and_exit_quarter_period_ramps=True)
    save(out/'INPUT_VALIDATION.json',proof)
    plan=dict(directory=str(out),binary=old_plan['genesis'],binary_sha256=sha(old_plan['genesis']),
        original_reference=str(source),input=str(new),PIC_input=str(pic),
        PIC_input_sha256=pic_sha,nodes=1,ranks=64,threads=1,case_limit_seconds=300,
        slurm_limit_seconds=480,endpoint_m=2.166,expected_dump_step=309,
        source_slices=[first,last],genesis_kept_groups_inclusive=[lower,upper],
        same_slice_grid_and_padding=True,selected_charge_C=float(charge),
        purpose='First-segment Genesis radiation estimate for same selected231slices as19325/19326',
        new_PIC_runs=0,new_noise=False,physical_acceptance=False)
    save(out/'plan.json',plan)
    paths=[new,out/'reference.in',out/'reference.lat',out/'INPUT_VALIDATION.json',out/'plan.json',
        *scripts.iterdir(),Path(plan['binary']),pic,campaign/'INPUT_SELECTION.json']
    save(out/'FROZEN.json',{str(p):sha(p) for p in paths})
    need(sha(original)==original_sha and sha(pic)==pic_sha, 'Original input changed')
    print(json.dumps(proof,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace',type=Path)
    parser.add_argument('out',type=Path)
    args=parser.parse_args()
    prepare(args.workspace,args.out)
