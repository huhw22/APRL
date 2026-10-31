#!/usr/bin/env python3
"""Build the diagnostic candidate separately; never replace the accepted binary."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve();assert os.environ.get('SLURM_JOB_ID')
ws=Path('/public/home/wfollowm/APRL_Workspace');src=ws/'source/aprl';build=ws/'build/particle-work-geometry-20261031'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
frozen=json.loads((root/'FROZEN.json').read_text())
def check():
    for p,h in frozen.items():assert sha(p)==h,p
def run(cmd,file,timeout):
    print('RUN',cmd,flush=True)
    with (root/file).open('x') as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
check();assert not build.exists()
runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
cmake='/public/home/wfollowm/software/cmake-3.31.12/bin/cmake'
run([cmake,'-S',str(src),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_COMPILER=/usr/bin/mpic++',
     '-DCMAKE_PREFIX_PATH='+str(runtime),'-DHDF5_ROOT='+str(runtime),'-DBUILD_TESTING=ON',
     '-DAPRL_ENABLE_OPENMP=ON','-DAPRL_ENABLE_MOVING_WINDOW=ON','-DAPRL_REQUIRE_PARALLEL_HDF5=ON'], 'configure.log',120)
run([cmake,'--build',str(build),'--parallel','12'],'compile.log',900)
cmd=['/usr/bin/mpic++','-std=c++11','-O2','-fopenmp','-DAPRL_HAS_OPENMP=1','-DAPRL_HAS_MOVING_WINDOW=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX',
     '-I'+str(src/'src'),'-I'+str(runtime/'include'),str(root/'scripts/geometry_checks.cpp'),str(build/'libdirect_eb_core.a'),
     str(runtime/'lib/libhdf5.so'),str(runtime/'lib/libyaml-cpp.so.0.8.0'),'-lfftw3','-Wl,-rpath,'+str(runtime/'lib'),'-o',str(build/'geometry_checks')]
run(cmd,'geometry-compile.log',120)
run([str(build/'geometry_checks')],'geometry-checks.json',60)
assert json.loads((root/'geometry-checks.json').read_text())['passed']
run(['ctest','--test-dir',str(build),'-L','required','--output-on-failure','-j','1'],'required.log',450)
check()
(root/'COMPLETE.json').write_text(json.dumps({'checks_passed':True,'job':os.environ['SLURM_JOB_ID'],'candidate':str(build/'aprl'),
    'candidate_sha256':sha(build/'aprl'),'production_binary_unchanged':sha(ws/'build/charge-matched-inlet-20261030/aprl')=='79d643bc8ade4db4d9e4488a367cb1bf5ef74e7c602dc3a2094462709aed9435',
    'source_frozen_hashes_pass':True,'scope':'Build, original required regressions, independent geometry-work checks; not physical acceptance'},indent=2)+'\n')
print('COMPLETE',flush=True)
