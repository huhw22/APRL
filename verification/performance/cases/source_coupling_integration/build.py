#!/usr/bin/env python3
import hashlib,json,os,subprocess,sys
from pathlib import Path
import numpy as np
root=Path(sys.argv[1]).resolve();assert os.environ.get('SLURM_JOB_ID');plan=json.loads((root/'plan.json').read_text());build=Path(plan['build']);source=Path(plan['source']);runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab');cmake='/public/home/wfollowm/software/cmake-3.31.12/bin/cmake'
frozen=json.loads((root/'FROZEN.json').read_text())
def check():
 for p,h in frozen.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
def run(cmd,name,limit):
 print('RUN',name,flush=True)
 with (root/name).open('x') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=limit)
check();assert not build.exists()
run([cmake,'-S',str(source),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_COMPILER=/usr/bin/mpic++','-DCMAKE_PREFIX_PATH='+str(runtime),'-DHDF5_ROOT='+str(runtime),'-DBUILD_TESTING=ON','-DAPRL_ENABLE_OPENMP=ON','-DAPRL_ENABLE_MOVING_WINDOW=ON','-DAPRL_REQUIRE_PARALLEL_HDF5=ON'],'configure.log',120)
run([cmake,'--build',str(build),'--parallel','12'],'compile.log',900)
here=source/'verification/performance/cases/source_coupling_integration';common=['/usr/bin/mpic++','-std=c++11','-O3','-fopenmp','-DAPRL_HAS_OPENMP=1','-DAPRL_HAS_MOVING_WINDOW=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(source/'src'),'-I'+str(runtime/'include')]
for name in ['checks','field_checks']:
 run([*common,str(here/(name+'.cpp')),str(build/'libdirect_eb_core.a'),str(runtime/'lib/libhdf5.so'),str(runtime/'lib/libyaml-cpp.so.0.8.0'),'-lfftw3','-Wl,-rpath,'+str(runtime/'lib'),'-o',str(build/name)],name+'-compile.log',120)
run([str(build/'checks')],'checks.json',60);assert json.loads((root/'checks.json').read_text())['checks_pass']
for ranks in [1,2,4]:
 run(['mpiexec','-n',str(ranks),str(build/'field_checks'),str(root/('fields-'+str(ranks)))],'fields-'+str(ranks)+'.json',90)
 assert json.loads((root/('fields-'+str(ranks)+'.json')).read_text())['passed']
 for file in (root/('fields-'+str(ranks))).glob('discrete-field-rank-*.csv'):
  x=np.genfromtxt(file,delimiter=',',names=True)
  for field,unit in [('electric','V_per_m'),('magnetic','T')]:assert np.all(x['max_'+field+'_prediction_error_'+unit]<=1024*np.finfo(float).eps*np.maximum(x[field+'_prediction_scale_'+unit],1e-300)),str(file)
run(['ctest','--test-dir',str(build),'-L','required','--output-on-failure','-j','1'],'required.log',600)
check();(root/'COMPLETE.json').write_text(json.dumps({'passed':True,'job':os.environ['SLURM_JOB_ID'],'binary':str(build/'aprl'),'binary_sha256':hashlib.sha256((build/'aprl').read_bytes()).hexdigest(),'frozen_snapshot':True,'physical_acceptance':False},indent=2)+'\n');print('COMPLETE',flush=True)
