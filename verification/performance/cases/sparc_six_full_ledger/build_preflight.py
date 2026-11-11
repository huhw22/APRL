import hashlib,json,os,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np

R=Path(sys.argv[1]).resolve();P=json.loads((R/'PLAN.json').read_text())
assert os.environ.get('SLURM_JOB_ID')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
for p,d in json.loads((R/'FROZEN.json').read_text()).items():assert sha(p)==d,p
runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');src=Path(P['source'])/'src'
common=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DMPICH_SKIP_MPICXX','-DOMPI_SKIP_MPICXX','-I'+str(src),'-isystem',str(runtime/'include')]
def command(cmd,log,limit=900,env=None):
    with Path(log).open('x') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,check=True,timeout=limit)
command(common+['-DAPRL_BUILD_DESCRIPTION="GNU-9.4.0,Release,HDF5-1.14.6"','-DAPRL_SOURCE_REVISION="diagnostic-side-events"','-c',str(src/'simulation.cpp'),'-o',str(R/'bin/simulation.o')],R/'compile.log')
libs=[P['library'],'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-lfftw3']
command(common+[P['main_object'],str(R/'bin/simulation.o')]+libs+['-o',P['binary']],R/'link.log',120)
command(common+[str(R/'scripts/side_unit.cpp')]+libs+['-o',str(R/'bin/side_unit')],R/'side_unit_compile.log',120)
unit=R/'side-unit';unit.mkdir();command([str(R/'bin/side_unit'),str(unit)],R/'side-unit.log',60)
with (unit/'particle-sides-rank-00000.csv').open() as f:lines=f.read().splitlines()
assert len(lines)==3,lines
env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')}
env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores',
           OMPI_MCA_btl='self,vader,tcp',OPENBLAS_NUM_THREADS='1',HDF5_USE_FILE_LOCKING='FALSE')
srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
probe='/public/home/wfollowm/APRL_Workspace/build/cluster-tools/binding_probe'
command(srun+[probe,'16','8'],R/'binding.log',90,env)
results=[]
for case in P['preflight']:
    d=Path(case['directory']);e=env.copy()
    if case['observer']:e.update(APRL_LAB_WORLDTUBE=str(d/'worldtube'),APRL_LAB_FIELD_AUDIT=str(d/'lab-audit'),APRL_WHOLE_COWAN_AUDIT=str(d/'lab-audit'))
    t=time.monotonic()
    with (d/'aprl.log').open('x') as f:
        process=subprocess.run(['timeout','--signal=TERM','--kill-after=30','1000']+srun+[P['binary'],str(d/'run.yaml')],stdout=f,stderr=subprocess.STDOUT,env=e,timeout=1070)
    log=(d/'aprl.log').read_text()
    assert process.returncode not in [0,124,137]
    assert 'mesh.duration was exhausted before the configured physical stop' in log
    assert 'phase=duration-exhaustion-readiness-sync-complete' in log
    assert 'Moving window completed: steps=40' in log
    with h5py.File(d/'ledger/energy-ledger.h5') as f:
        assert int(f['energy_ledger/complete'][()])==0
        assert int(f['energy_ledger/records'][-1]['step'])==40
    if case['observer']:
        assert (d/'worldtube/COMPLETE').read_text().split()==['0','40']
        for rank in range(16):
            for prefix in ['lab-field','paired-momentum','whole-cowan']:
                assert (d/'lab-audit'/f'{prefix}-rank-{rank:05d}.complete').read_text().split()==['0','40']
        data=np.fromfile(d/'worldtube/field_lines.bin',dtype='<f8')
        assert len(data)==40*P['active_nz']*6 and np.isfinite(data).all()
        assert len(list((d/'worldtube').glob('particle-sides-rank-*.csv')))==16
    results.append(dict(case=case['name'],diagnostic_preflight_passed=True,
        intentional_duration_cutoff=True,physical_run_complete=False,seconds=time.monotonic()-t))
    print(json.dumps(results[-1]),flush=True)
a=R/'preflight/binomial3-off';b=R/'preflight/binomial3-on';same=[]
for p in sorted((a/'detectors').glob('*.h5')):
    with h5py.File(p) as f,h5py.File(b/'detectors'/p.name) as g:
        group=f[next(iter(f))];other=g[next(iter(g))]
        for key in ['records','fields','time_s']:
            if key in group:assert np.array_equal(group[key][:],other[key][:]),(p,key)
    same.append(p.name)
with h5py.File(a/'ledger/energy-ledger.h5') as f,h5py.File(b/'ledger/energy-ledger.h5') as g:
    x=f['energy_ledger/records'][:];y=g['energy_ledger/records'][:]
    for key in x.dtype.names:
        if not key.endswith('_rank_seconds'):assert np.array_equal(x[key],y[key]),key
save(R/'PREFLIGHT_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],
    binary=P['binary'],binary_sha256=sha(P['binary']),cases=results,
    off_on_scientific_h5_bitwise=same,side_event_units_passed=True,
    original_lab_worldtube_MPI_unit_reused='/public/home/wfollowm/APRL_Workspace/results/lab_worldtube_balance_20261112_r1/BUILD_COMPLETE.json',
    physics_core_unchanged=True))
print('PREFLIGHT_COMPLETE',flush=True)
