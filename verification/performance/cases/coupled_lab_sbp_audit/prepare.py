import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();r.mkdir(exist_ok=True,parents=True);require(not (r/'PLAN.json').exists(),'Existing frozen plan')
    template=Path(__file__).resolve().parent;frozen=r/'frozen_source';frozen.mkdir()
    for p in template.iterdir():
        if p.is_file() and p.suffix in ['.h','.cpp','.py','.md']:shutil.copy2(p,frozen/p.name)
    require(sha(w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a')=='1e7939f053ab22370be17718aefe50efa3f818989eac3a5d039c0b5cc11db18a','Production field archive changed')
    from overlay import apply
    source=r/'source';base=w/'results/lab_energy_balance_20261102/r5/source';apply(base,source,frozen)
    build=a.build.resolve() if a.build else w/'build/coupled-lab-sbp-20261102';build.mkdir(exist_ok=True);runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');archive=build/'libdirect_eb_core.a';require(not archive.exists(),'Existing candidate archive');shutil.copy2(w/'results/lab_energy_balance_20261102/r5/bin/libdirect_eb_core.a',archive)
    cmd=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_BUILD_DESCRIPTION="GNU-9.4.0,Release,HDF5-1.14.6"','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DAPRL_SOURCE_REVISION="unknown"','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(source/'src'),'-isystem',str(runtime/'include'),'-c',str(source/'src/simulation.cpp'),'-o',str(build/'simulation.cpp.o')]
    with open(r/'build-simulation.log','w') as out:subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,check=True)
    subprocess.run(['ar','r',str(archive),str(build/'simulation.cpp.o')],check=True)
    binary=build/'aprl';link=['/usr/bin/mpic++','-O3','-DNDEBUG',str(w/'build/source-coupling-exit-r3-20261031/CMakeFiles/aprl.dir/src/main.cpp.o'),'-o',str(binary),'-Wl,-rpath,'+str(runtime/'lib'),str(archive),'-L'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-fopenmp','-lfftw3'];subprocess.run(link,check=True)
    check=build/'check_sources';cmdcheck=['/usr/bin/mpic++','-O3','-std=c++11','-fopenmp','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(w/'source/aprl/src'),'-I'+str(frozen),str(frozen/'check_sources.cpp'),str(w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a'),'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-lfftw3','-o',str(check)]
    with open(r/'build-check.log','w') as out:subprocess.run(cmdcheck,stdout=out,stderr=subprocess.STDOUT,check=True)
    binaries={'wave':str(binary),'check':str(check)};commands={'simulation':cmd,'link':link,'check':cmdcheck};probe=w/'build/cluster-tools/binding_probe';baseline=w/'results/lab_energy_balance_20261102/r5';cases=[]
    for name in ['undulator-on','drift-on']:
        d=r/name;d.mkdir();text=(baseline/name/'run.yaml').read_text().replace(str(baseline),str(r));(d/'run.yaml').write_text(text)
        for sub in ['detectors','ledger','trajectory','lab-audit','momentum']:(d/sub).mkdir()
        (d/'lab-audit/synchronized').mkdir();cases.append(dict(name=name,directory=str(d),timeout_seconds=1200))
    (r/'bin').mkdir();(r/'bin/replay').symlink_to(w/'results/lab_energy_balance_20261102/r5/bin/replay')
    plan=dict(workspace=str(w),output=str(r),programs=binaries,binding_probe=str(probe),baseline=str(baseline),cases=cases,production_physics_changed=False,resources=dict(nodes=1,mpi=16,threads=8,exclusive_cpus=128));save(r/'PLAN.json',plan)
    batch='''#!/usr/bin/env bash
#SBATCH --job-name=coupled-lab-sbp-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:45:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=ROOT_PATH
source WORKSPACE_PATH/cluster.env
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 2640 "$PYTHON_BIN" frozen_source/run.py > driver.log 2>&1
'''.replace('ROOT_PATH',str(r)).replace('WORKSPACE_PATH',str(w));(r/'run.sbatch').write_text(batch)
    deps=list(frozen.glob('*'))+list((source/'src').glob('*.h'))+list((source/'src').glob('*.cpp'))+[(r/c['name']/'run.yaml') for c in cases]+[Path(x) for x in binaries.values()]+[archive,probe,r/'PLAN.json',r/'run.sbatch']+list((w/'source/aprl/src').glob('*.h'))
    save(r/'FROZEN.json',{str(p):sha(p) for p in deps if p.is_file()});save(r/'SOURCE_BUILD.json',dict(commands=commands,archive_sha256=sha(archive),program_sha256={k:sha(Path(v)) for k,v in binaries.items()},physics_kernel_rebuilt=False,production_binary_replaced=False))
    print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
