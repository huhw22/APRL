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
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();r.mkdir(exist_ok=True,parents=True);require(not (r/'PLAN.json').exists(),'Existing frozen plan')
    template=Path(__file__).resolve().parent;frozen=r/'frozen_source';frozen.mkdir()
    for p in template.iterdir():
        if p.is_file() and p.suffix in ['.h','.cpp','.py','.md']:shutil.copy2(p,frozen/p.name)
    build=w/'build/yee-lab-sbp-20261102';build.mkdir(exist_ok=True);archive=w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a';require(sha(archive)=='1e7939f053ab22370be17718aefe50efa3f818989eac3a5d039c0b5cc11db18a','Production archive changed')
    runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');commands={};binaries={}
    for name,source in [('wave','wave.cpp'),('check','check_geometry.cpp')]:
        binary=build/('yee_sbp_'+name);require(not binary.exists(),'Refuse to overwrite an existing program '+str(binary));cmd=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(w/'source/aprl/src'),'-I'+str(frozen),'-I'+str(runtime/'include'),str(frozen/source),str(archive),'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-lfftw3','-o',str(binary)]
        with open(r/('build-'+name+'.log'),'w') as out:subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,check=True)
        commands[name]=cmd;binaries[name]=str(binary)
    probe=w/'build/cluster-tools/binding_probe';baseline=w/'results/lab_energy_balance_20261102/vacuum_propagation_r4/yee-oblique'
    plan=dict(workspace=str(w),output=str(r),programs=binaries,binding_probe=str(probe),baseline=str(baseline),cases=[dict(name='yee-oblique',mode='oblique',refine=1,timeout_seconds=1000)],production_physics_changed=False,resources=dict(nodes=1,mpi=16,threads=8,exclusive_cpus=128));save(r/'PLAN.json',plan)
    batch='''#!/usr/bin/env bash
#SBATCH --job-name=yee-lab-sbp-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=ROOT_PATH
source WORKSPACE_PATH/cluster.env
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 1140 "$PYTHON_BIN" frozen_source/run.py > driver.log 2>&1
'''.replace('ROOT_PATH',str(r)).replace('WORKSPACE_PATH',str(w));(r/'run.sbatch').write_text(batch)
    deps=list(frozen.glob('*'))+[Path(x) for x in binaries.values()]+[archive,probe,r/'PLAN.json',r/'run.sbatch']+list((w/'source/aprl/src').glob('*.h'))
    save(r/'FROZEN.json',{str(p):sha(p) for p in deps if p.is_file()});save(r/'SOURCE_BUILD.json',dict(commands=commands,archive_sha256=sha(archive),program_sha256={k:sha(Path(v)) for k,v in binaries.items()},physics_kernel_rebuilt=False,production_binary_replaced=False))
    print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
