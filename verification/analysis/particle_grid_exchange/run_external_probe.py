#!/usr/bin/env python3
"""Run the unchanged external pusher on three actual input records."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
import h5py

root=Path(sys.argv[1]).resolve();assert os.environ.get('SLURM_JOB_ID')
for p,h in json.loads((root/'FROZEN.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
ws=Path('/public/home/wfollowm/APRL_Workspace');src=ws/'source/aprl';runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
command=['/usr/bin/mpic++','-std=c++11','-O2','-fopenmp','-DAPRL_HAS_OPENMP=1','-DAPRL_HAS_MOVING_WINDOW=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX',
  '-I'+str(src/'src'),'-I'+str(runtime/'include'),str(root/'scripts/external_magnet_probe.cpp'),
  str(ws/'build/charge-matched-inlet-20261030/libdirect_eb_core.a'),str(runtime/'lib/libhdf5.so'),str(runtime/'lib/libyaml-cpp.so.0.8.0'),
  '-lfftw3','-Wl,-rpath,'+str(runtime/'lib'),'-o',str(root/'external_magnet_probe')]
(root/'COMPILE_COMMAND.json').write_text(json.dumps(command,indent=2)+'\n')
subprocess.run(command,check=True,timeout=120)
input_file=ws/'results/elegant_4m_analysis_20261024/original/particles.h5';results=[];inputs=[]
with h5py.File(input_file) as f:
    d=f['particles/records']
    for k in [0,len(d)//2,len(d)-1]:
        row=d[k];xy=row['plane_position_m'][:2];u=row['proper_velocity'];inputs.append({'row':k,'source_id':int(row['source_id']),'xy_m':xy.tolist(),'proper_velocity':u.tolist()})
        cmd=[str(root/'external_magnet_probe'),str(k),*[repr(float(v)) for v in [*xy,*u]]]
        out=subprocess.check_output(cmd,text=True,timeout=120)
        rows=[json.loads(t) for t in out.splitlines()];assert len(rows)==8
        assert all(v['lab_z_final_m']>=4.5 for v in rows)
        assert all(abs(v['lab_energy_change_eV'])<1e-3 for v in rows if not v['magnet'])
        results.extend(rows)
        print(out,flush=True,end='')
for p,h in json.loads((root/'FROZEN.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
(root/'RESULT.json').write_text(json.dumps({'integrity_passed':True,'job':os.environ['SLURM_JOB_ID'],'input_samples':inputs,'results':results,
  'new_full_PIC_simulations':0,'physics_method_changes':False,'scope':'Three sampled original particles; external-only trajectories. Not a full beam or self-consistent field test.'},indent=2)+'\n')
