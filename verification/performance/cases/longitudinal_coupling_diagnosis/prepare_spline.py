"""Freeze an independent longitudinal source/gather model, not production code."""
import argparse,json,shutil,subprocess
from pathlib import Path
from prepare import sha,save
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--build',required=True,type=Path);p.add_argument('--campaign',choices=['all','individual','remainder'],default='all');a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();b=a.build.resolve();r.mkdir(parents=True,exist_ok=False);b.mkdir(parents=True,exist_ok=False);t=Path(__file__).resolve().parent;f=r/'frozen_source';f.mkdir()
 for name in ['spline_1d.cpp','run_spline.py','prepare.py','prepare_spline.py']:shutil.copy2(t/name,f/name)
 old=w/'results/cowan_matched_coupling_20261102/build_r2/source/src';runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');archive=w/'build/source-coupling-exit-r3-20261031/libdirect_eb_core.a';program=b/'spline_1d';cmd=['/usr/bin/mpic++','-O3','-DNDEBUG','-std=c++11','-fopenmp','-DAPRL_HAS_MOVING_WINDOW=1','-DAPRL_HAS_OPENMP=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(old),'-isystem',str(runtime/'include'),str(f/'spline_1d.cpp'),str(archive),'-L'+str(runtime/'lib'),'-Wl,-rpath,'+str(runtime/'lib'),'-lhdf5','-lyaml-cpp','-lpthread','-fopenmp','-lfftw3','-o',str(program)]
 with (r/'build.log').open('w') as out:subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,check=True)
 shutil.copy2(w/'source/aprl/verification/performance/cases/coupled_lab_sbp_audit/launcher.py',f/'launcher.py');cases=[]
 def case(name,degree=1,comp=0,rz=1,pairs=16,phase=.31,u=.23,duration=1,amplitude=.15):cases.append(dict(name=name,args=[degree,comp,rz,pairs,phase,u,duration,amplitude],steps=160*rz*duration,timeout_s=180))
 case('linear_baseline');case('cubic_uncompensated',degree=3);case('cubic_compensated_3pt',degree=3,comp=1);case('cubic_compensated_5pt',degree=3,comp=2)
 case('linear_dense',pairs=96);case('cubic_compensated_dense',degree=3,comp=2,pairs=96);case('linear_z2',rz=2,pairs=32);case('cubic_compensated_z2',degree=3,comp=2,rz=2,pairs=32)
 case('linear_long',duration=10);case('cubic_compensated_long',degree=3,comp=2,duration=10);case('linear_stationary',u=0.);case('cubic_compensated_stationary',degree=3,comp=2,u=0.)
 case('quintic_compensated_5pt',degree=5,comp=2);case('quintic_compensated_dense',degree=5,comp=2,pairs=96);case('quintic_compensated_stationary',degree=5,comp=2,u=0.);case('quintic_compensated_long',degree=5,comp=2,duration=10)
 case('analytic_linear',u=0.,pairs=96,amplitude=.0015);case('analytic_cubic_continuum',degree=3,comp=3,u=0.,pairs=96,amplitude=.0015);case('analytic_quintic_continuum',degree=5,comp=4,u=0.,pairs=96,amplitude=.0015)
 if a.campaign=='individual':
  cases=[c for c in cases if c['name'].startswith('analytic')]
  case('linear_stationary_long',u=0.,duration=20)
  case('quintic_compensated_stationary_long',degree=5,comp=2,u=0.,duration=20)
 if a.campaign=='remainder':
  cases=[]
  case('linear_stationary',u=0.,duration=20)
  case('quintic_stationary',degree=5,comp=2,u=0.,duration=20)
  case('quintic_stationary_z2',degree=5,comp=2,u=0.,duration=20,rz=2)
  case('quintic_stationary_dense',degree=5,comp=2,u=0.,duration=20,pairs=96)
  case('quintic_stationary_dense_z2',degree=5,comp=2,u=0.,duration=20,pairs=96,rz=2)
  case('quintic_stationary_continuum',degree=5,comp=4,u=0.,duration=20)
  for c in cases:c['args'].append(1);c['timeout_s']=300
 save(r/'PLAN.json',dict(workspace=str(w),program=str(program),cases=cases,scope='Independent periodic longitudinal sheet model; source/gather mathematical exploration only',production_deposition_changed=False,production_pusher_changed=False,Cowan_kernel_changed=False,transverse_energy_or_momentum_correction_applied=False,all_cases_c_dt_over_dz=1.0,original_pusher_used=True))
 batch='''#!/usr/bin/env bash
#SBATCH --job-name=longitudinal-shape-model
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
timeout --signal=TERM --kill-after=20 840 "$PYTHON_BIN" frozen_source/run_spline.py > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(batch);save(r/'SOURCE_BUILD.json',dict(command=cmd,program_sha256=sha(program),archive_sha256=sha(archive)));files=[*f.iterdir(),r/'run.sbatch',r/'PLAN.json',program,archive,*[q for q in old.iterdir() if q.suffix in ['.h','.cpp']]];save(r/'FROZEN.json',{str(q):sha(q) for q in files if q.is_file()});print(r/'run.sbatch')
if __name__=='__main__':main()
