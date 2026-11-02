"""Build a declared test-only magnetic gather, retaining the audited field code."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def replace(s,a,b):
    require(s.count(a)==1,'Patch anchor '+a[:60]);return s.replace(a,b)
def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--build',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();build=a.build.resolve();r.mkdir(parents=True,exist_ok=False);build.mkdir(parents=True,exist_ok=False);t=Path(__file__).resolve().parent;base=w/'results/lab_energy_balance_20261102/coupled_sbp_r2';frozen=r/'frozen_source';frozen.mkdir()
    require(json.loads((base/'RESULT.json').read_text())['cases']['undulator-on']['complete'],'Baseline incomplete')
    require(sha(w/'build/source-coupling-exit-r3-20261031/aprl')=='865ebef5b79b6b8670d40c81080f15f916ab1938033afcabd050bcd3c8198752','Production program changed')
    for p in t.iterdir():
        if p.is_file():shutil.copy2(p,frozen/p.name)
    helpers=w/'source/aprl/verification/performance/cases/coupled_lab_sbp_audit'
    for name in ['launcher.py','post.py','paired_post.py','work_replay_post.py','review.py']:shutil.copy2(helpers/name,frozen/name)
    p=frozen/'paired_post.py';s=p.read_text();s=replace(s,"d=R/name;old=R.parent/'r1'/('undulator-off' if name=='undulator-on' else 'drift-on');assert identity(old,d)","d=R/name # Declared algorithm change; no false bitwise identity assertion.");s=replace(s,'original_physics_bitwise=True','original_physics_bitwise=False');p.write_text(s)
    p=frozen/'post.py';s=p.read_text();s=replace(s,'original_science_bitwise_identity=True','original_science_bitwise_identity=False, candidate_magnetic_gather_changed=True');s=replace(s,'production_physics_changed=False,physical_lab_conservation_accepted=False','production_physics_changed=False,candidate_algorithm_changed=True,physical_lab_conservation_accepted=False');p.write_text(s)
    source=r/'source';shutil.copytree(base/'source',source);shutil.copy2(t/'momentum_z_gather.h',source/'src/momentum_z_gather.h');p=source/'src/source_coupling.cpp';s=p.read_text();s=replace(s,'#include "source_coupling_math.h"','#include "source_coupling_math.h"\n#include "momentum_z_gather.h"');require(s.count('openOrbitMagnetic(')==3,'Two physical magnetic gather calls plus original constraint observer');s=s.replace('openOrbitMagnetic(','momentumOrbitMagnetic(',2);p.write_text(s)
    shutil.copy2(t/'native_two_magnetic_halo.h',source/'src/native_two_magnetic_halo.h');math=source/'src/source_coupling_math.h';v=math.read_text();first=v.index('class NativeHalo {');last=v.index('struct MagneticSample',first);math.write_text(v[:first]+'#include "native_two_magnetic_halo.h"\n'+v[last:])
    archive=build/'libdirect_eb_core.a';shutil.copy2(w/'build/coupled-lab-sbp-r2-20261102/libdirect_eb_core.a',archive);runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');obj=build/'source_coupling.cpp.o';common=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(source/'src'),'-isystem',str(runtime/'include')];compilecmd=common+['-c',str(p),'-o',str(obj)]
    with (r/'build.log').open('w') as out:subprocess.run(compilecmd,stdout=out,stderr=subprocess.STDOUT,check=True)
    subprocess.run(['ar','r',str(archive),str(obj)],check=True);binary=build/'aprl';libs=['-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-fopenmp','-lfftw3'];linkcmd=['/usr/bin/mpic++','-O3',str(w/'build/source-coupling-exit-r3-20261031/CMakeFiles/aprl.dir/src/main.cpp.o'),str(archive),*libs,'-o',str(binary)];subprocess.run(linkcmd,check=True)
    check=build/'checks';checkcmd=common+['-I'+str(frozen),str(frozen/'checks.cpp'),str(w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a'),*libs,'-o',str(check)]
    with (r/'build-check.log').open('w') as out:subprocess.run(checkcmd,stdout=out,stderr=subprocess.STDOUT,check=True)
    cases=[]
    for name in ['undulator-on','drift-on']:
        d=r/name;d.mkdir();(d/'run.yaml').write_text((base/name/'run.yaml').read_text().replace(str(base),str(r)))
        for sub in ['detectors','ledger','trajectory','lab-audit','lab-audit/synchronized','momentum']:(d/sub).mkdir()
        cases.append(dict(name=name,timeout_seconds=1200))
    (r/'bin').mkdir();(r/'bin/replay').symlink_to(w/'results/lab_energy_balance_20261102/r5/bin/replay');plan=dict(workspace=str(w),baseline=str(base),output=str(r),programs={'wave':str(binary),'check':str(check)},binding_probe=str(w/'build/cluster-tools/binding_probe'),cases=cases,algorithm='momentum-z-adjoint-magnetic-gather',candidate_algorithm_changed=True,production_program_replaced=False);save(r/'PLAN.json',plan)
    batch='''#!/usr/bin/env bash
#SBATCH --job-name=lab-adjoint-B-test
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 1440 "$PYTHON_BIN" frozen_source/run.py > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(batch)
    deps=list(frozen.iterdir())+list((source/'src').glob('*.h'))+list((source/'src').glob('*.cpp'))+[r/c['name']/'run.yaml' for c in cases]+[binary,check,archive,r/'PLAN.json',r/'run.sbatch'];save(r/'FROZEN.json',{str(p):sha(p) for p in deps if p.is_file()});save(r/'SOURCE_BUILD.json',dict(compile=compilecmd,link=linkcmd,check=checkcmd,binary_sha256=sha(binary),archive_sha256=sha(archive),base_archive_sha256=sha(w/'build/coupled-lab-sbp-r2-20261102/libdirect_eb_core.a'),changed_object='source_coupling.cpp.o',production_replaced=False));print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
