#!/usr/bin/env python3
"""Frozen, serial small-case campaign. No production binary or input changes."""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
root=Path(sys.argv[1]).resolve();assert os.environ.get('SLURM_JOB_ID')
plan=json.loads((root/'plan.json').read_text());frozen=json.loads((root/'FROZEN.json').read_text())
def check():
 for p,h in frozen.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def run(cmd,log,limit=180):
 print('RUN',log,flush=True)
 with (root/log).open('x') as f:subprocess.run(['timeout','--signal=TERM','--kill-after=10',str(limit),*cmd],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=limit+20)
check();build=Path(plan['build']);build.mkdir();ws=Path(plan['workspace']);src=ws/'source/aprl';runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
common=['/usr/bin/mpic++','-std=c++11','-O3','-fopenmp','-DAPRL_HAS_OPENMP=1','-DAPRL_HAS_MOVING_WINDOW=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX','-I'+str(src/'src'),'-I'+str(root/'sources'),'-I'+str(runtime/'include')]
for name in ['checks','prototype','external']:
 run([*common,str(root/'sources'/(name+'.cpp')),plan['library'],str(runtime/'lib/libhdf5.so'),str(runtime/'lib/libyaml-cpp.so.0.8.0'),'-lfftw3','-Wl,-rpath,'+str(runtime/'lib'),'-o',str(build/name)],f'compile-{name}.log',120)
launcher=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=1','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=20']
run([*launcher,str(build/'checks')],'checks.json',60);assert json.loads((root/'checks.json').read_text())['checks_pass']
run([*launcher,sys.executable,str(root/'scripts/binding.py'),str(root/'binding.json')],'binding.log',30)
binding=json.loads((root/'binding.json').read_text());assert binding['distinct_physical_cores']>=8
results={}
for case in plan['cases']:
 name=case['name'];directory=root/name;directory.mkdir();started=time.monotonic()
 run([*launcher,str(build/'prototype'),str(directory),case['kind'],case['mode'],str(case['particles']),str(case['steps']),str(int(case['audit'])),'8',str(case['dt_divisor']),str(case['nx'])],name+'/aprl.log',case['limit_seconds'])
 result=json.loads((directory/'RESULT.json').read_text());assert result['complete']
 if case['audit']:
  assert result['relative_continuity_residual_max']<5e-12
  assert (directory/'discrete-field-rank-00000.complete').read_text().split()==['1',str(case['steps'])]
 assert result['magnetic_CPML_calls']==case['steps']==result['electric_CPML_calls']
 results[name]=result;print('COMPLETE',name,'wall_s',time.monotonic()-started,'step_s',result['seconds_per_step'],'outer',result['maximum_outer_iterations'],flush=True)
 if name=='vacuum-coupled':assert results['vacuum-baseline']['field_hash']==result['field_hash']
 save(root/'PARTIAL_CLOSED_CASES.json',results)
run([*launcher,str(build/'external'),str(root/'external.csv')],'external.log',180)
check();save(root/'COMPLETE.json',{'checks_passed':True,'job':os.environ['SLURM_JOB_ID'],'candidate_is_independent_prototype':True,'production_core_unchanged':True,'cases':results,'binaries_sha256':{n:hashlib.sha256((build/n).read_bytes()).hexdigest() for n in ['checks','prototype','external']}})
print('ALL CASES COMPLETE',flush=True)
