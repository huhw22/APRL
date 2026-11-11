import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
import h5py
R=Path(sys.argv[1]).resolve();name=sys.argv[2]
P=json.loads((R/'PLAN.json').read_text());gate=json.loads((R/'PREFLIGHT_COMPLETE.json').read_text())
assert gate['complete'] and hashlib.sha256(Path(P['binary']).read_bytes()).hexdigest()==gate['binary_sha256']
for p,digest in json.loads((R/'FROZEN.json').read_text()).items():
    assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,p
case=next(x for x in P['cases'] if x['name']==name);d=Path(case['directory'])
env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores',
    OMPI_MCA_btl='self,vader,tcp',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',
    HDF5_USE_FILE_LOCKING='FALSE',APRL_LAB_WORLDTUBE=str(d/'worldtube'),
    APRL_LAB_FIELD_AUDIT=str(d/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(d/'lab-audit'))
srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
with (d/'binding.log').open('x') as f:subprocess.run(srun+['/public/home/wfollowm/APRL_Workspace/build/cluster-tools/binding_probe','16','8'],env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=90)
record=dict(job_id=os.environ['SLURM_JOB_ID'],binary=P['binary'],binary_sha256=gate['binary_sha256'],
    input_file=P['input_file'],input_sha256=P['input_sha256'],start_unix=time.time(),
    full_ledgers_enabled=['boost-discrete','lab-field-momentum-stress','whole-cowan','lab-worldtube','lab-particle-sides','substep-scalar-lab-planes'])
(d/'STARTED.json').write_text(json.dumps(record,indent=2)+'\n')
t=time.monotonic()
with (d/'aprl.log').open('x') as f:
    outcome=subprocess.run(['timeout','--signal=TERM','--kill-after=30','170400']+srun+[P['binary'],str(d/'run.yaml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
record.update(exit_code=outcome.returncode,seconds=time.monotonic()-t)
(d/'aprl.log.launch.json').write_text(json.dumps(record,indent=2)+'\n')
assert outcome.returncode==0
log=(d/'aprl.log').read_text()
for marker in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete',
    'Moving window completed: steps='+str(P['steps']), '[lab-worldtube]', '[lab-particle-sides]', '[lab-field-audit]']:
    assert marker in log,marker
assert not re.search('Simulation failure|Segmentation fault|OpenFabrics device',log)
for path in (d/'detectors').glob('*.h5'):
    with h5py.File(path) as f:
        g=f[next(iter(f))];assert int(g['complete'][()])==1,path
with h5py.File(d/'ledger/energy-ledger.h5') as f:
    assert int(f['energy_ledger/complete'][()])==1
    assert len(f['energy_ledger/records'])==P['steps']+1
assert (d/'worldtube/COMPLETE').read_text().split()==['1',str(P['steps'])]
for rank in range(16):
    for prefix in ['lab-field','paired-momentum','whole-cowan']:
        assert (d/'lab-audit'/f'{prefix}-rank-{rank:05d}.complete').read_text().split()==['1',str(P['steps'])]
result=dict(complete=True,case=name,run=record,steps=P['steps'],analysis_pending=True,
    raw_worldtube_bytes=(d/'worldtube/field_lines.bin').stat().st_size)
(d/'PIC_COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
