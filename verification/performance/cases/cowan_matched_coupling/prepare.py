"""Materialize a test-only matched Cowan gather; production files unchanged."""
import argparse,json,hashlib,shutil,subprocess
from pathlib import Path
def require(v,msg):
    if not v:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def replace(s,a,b):
    require(s.count(a)==1,'Overlay anchor: '+a[:50]);return s.replace(a,b,1)
def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--build',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();b=a.build.resolve();r.mkdir(parents=True,exist_ok=False);b.mkdir(parents=True,exist_ok=False);t=Path(__file__).resolve().parent;frozen=r/'frozen_source';frozen.mkdir();base=w/'results/lab_energy_balance_20261102/coupled_sbp_r2';source=r/'source';shutil.copytree(base/'source',source)
    require(sha(w/'build/source-coupling-exit-r3-20261031/aprl')=='865ebef5b79b6b8670d40c81080f15f916ab1938033afcabd050bcd3c8198752','Production binary changed')
    for p in t.iterdir():
        if p.is_file():shutil.copy2(p,frozen/p.name)
    for name in ['metric_gather.h','native_metric_halo.h']:shutil.copy2(t/name,source/'src'/name)
    p=source/'src/source_coupling_math.h';s=p.read_text();i=s.index('class NativeHalo {');j=s.index('struct MagneticSample',i);p.write_text(s[:i]+'#include "native_metric_halo.h"\n'+s[j:])
    p=source/'src/source_coupling.cpp';s=p.read_text();s=replace(s,'#include "source_coupling_math.h"','#include "source_coupling_math.h"\n#include "metric_gather.h"');require(s.count('openOrbitElectric(')==2 and s.count('openOrbitMagnetic(')==3,'Gather call count');s=s.replace('openOrbitElectric(','matchedOpenOrbitElectric(').replace('openOrbitMagnetic(','matchedOrbitMagnetic(').replace('orbitElectric(','matchedOrbitElectric(');p.write_text(s)
    runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');archive=b/'libdirect_eb_core.a';shutil.copy2(w/'build/coupled-lab-sbp-r2-20261102/libdirect_eb_core.a',archive);common=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(source/'src'),'-isystem',str(runtime/'include')];libs=['-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-fopenmp','-lfftw3'];obj=b/'source_coupling.cpp.o';cmd=common+['-c',str(p),'-o',str(obj)]
    with (r/'build-coupling.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    subprocess.run(['ar','r',str(archive),str(obj)],check=True);programs={};commands={'coupling':cmd}
    cmd=['/usr/bin/mpic++','-O3',str(w/'build/source-coupling-exit-r3-20261031/CMakeFiles/aprl.dir/src/main.cpp.o'),str(archive),*libs,'-o',str(b/'aprl')];subprocess.run(cmd,check=True);programs['aprl']=str(b/'aprl');commands['link']=cmd
    for name in ['checks','periodic']:
        cmd=common+['-I'+str(frozen),str(frozen/(name+'.cpp')),str(w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a'),*libs,'-o',str(b/name)]
        with (r/('build-'+name+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        programs[name]=str(b/name);commands[name]=cmd
    shutil.copy2(t.parent/'coupled_lab_sbp_audit/launcher.py',frozen/'launcher.py');plan=dict(workspace=str(w),root=str(r),programs=programs,binding_probe=str(w/'build/cluster-tools/binding_probe'),production_replaced=False,Cowan_propagator_changed=False,self_field_gather_changed=True,reference_supply_changed=False,prescribed_fields_changed=False);save(r/'PLAN.json',plan)
    batch='''#!/usr/bin/env bash
#SBATCH --job-name=cowan-matched-unit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 540 "$PYTHON_BIN" frozen_source/run_units.py > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(batch);files=list(frozen.iterdir())+list((source/'src').glob('*.cpp'))+list((source/'src').glob('*.h'))+[Path(x) for x in programs.values()]+[archive,r/'PLAN.json',r/'run.sbatch'];save(r/'FROZEN.json',{str(p):sha(p) for p in files if p.is_file()});save(r/'SOURCE_BUILD.json',dict(commands=commands,archive_sha256=sha(archive),program_sha256={k:sha(v) for k,v in programs.items()},changed_archive_object='source_coupling.cpp.o',production_replaced=False));print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
