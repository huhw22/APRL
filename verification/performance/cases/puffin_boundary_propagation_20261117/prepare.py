"""Freeze a small source-free native Puffin propagation experiment."""
import hashlib
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

W=Path('/public/home/wfollowm/APRL_Workspace')
B=W/'build/puffin-eleven-four-slip-20261113-r2/puffin'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')

def main(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=False)
    for d in ['scripts','bin','cases']: (root/d).mkdir()
    here=Path(__file__).resolve().parent
    for f in here.iterdir():
        if f.is_file():shutil.copyfile(f,root/'scripts'/f.name)
    source=W/'results/puffin_eleven_four_slip_20261113_r2/source'
    source_files=[source/'puffin/lib/diffraction'/f for f in ['diffraction.f90','masks.f90','transforms.f90']]
    before={str(p):sha(p) for p in source_files+[B/'libpuffin_mpi.a',B/'puffin']}
    assert before[str(B/'puffin')]=='c4ecf60611e86628f7ce77b90c752f49dddf6daedb8ca2877ee1783ad5c5c54e'
    s=json.loads((W/'results/fel_independent_w355_20261114_r1/eleven-single-node-r1/INPUT_VALIDATED.json').read_text())['scale']
    original=shlex.split((B/'CMakeFiles/puffin.dir/link.txt').read_text())
    deps=original[next(i for i,x in enumerate(original) if x.startswith('-Wl,-rpath')):]
    deps=[str(B/'libpuffin_mpi.a') if x=='libpuffin_mpi.a' else x for x in deps]
    cmd=[original[0],'-ffree-line-length-none','-fopenmp','-O2','-I'+str(B),str(root/'scripts/native_probe.f90'),'-o',str(root/'bin/native_probe')]+deps
    p=subprocess.run(cmd,capture_output=True,text=True)
    (root/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    cases=[dict(name='plane-preflight',nx=81,ny=81,beta=0,steps=1,stride=1,mode='plane'),
           dict(name='small-native',nx=81,ny=81,beta=1,steps=300,stride=5,mode='gaussian'),
           dict(name='small-no-absorption',nx=81,ny=81,beta=0,steps=300,stride=5,mode='gaussian'),
           dict(name='wide-reference',nx=401,ny=241,beta=1,steps=300,stride=5,mode='gaussian'),
           dict(name='wider-reference',nx=481,ny=321,beta=1,steps=300,stride=5,mode='gaussian')]
    for c in cases:(root/'cases'/c['name']).mkdir()
    plan=dict(date_local='2026-11-17',cases=cases,scale=s,
      same_frozen_library=True,source_hashes=before,probe_sha256=sha(root/'bin/native_probe'),
      scope='Actual native DiffractionStep and AbsorptionStep propagation, no particles; transverse-boundary test only.',
      initial_field=dict(intensity_sigma_um=200,x0_um=-300,kx_per_m=30000,central_ky_per_m=0),
      grid=dict(dx_um=40,nominal_step_m=.007,crop_nodes=81,roi_half_width_um=520),
      longitudinal=dict(mode='One wavelength periodic, 20 intervals; longitudinal end mask disabled by native periodic branch.',
        why='Isolate transverse absorption and finite-array wrap-around; this is not the production temporal FEL configuration.'),
      resources=dict(nodes=1,exclusive=True,mpi=4,threads_per_rank=1,slurm_minutes=20,case_timeout_s=180),
      no_production_source_or_input_change=True,no_FEL_rerun=True)
    save(root/'PLAN.json',plan)
    save(root/'FROZEN.json',{str(p):sha(p) for p in (root/'scripts').iterdir() if p.is_file()})
    print(root)

if __name__=='__main__':main(sys.argv[1])
