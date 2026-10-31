#!/usr/bin/env python3
import argparse,ast,hashlib,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);p.add_argument('--build',required=True,type=Path);a=p.parse_args()
ws=a.workspace.resolve();root=a.root.resolve();root.mkdir(parents=True);(root/'sources').mkdir();(root/'scripts').mkdir();here=Path(__file__).resolve().parent
for name in ['source_coupling.h','prototype.cpp','checks.cpp','external.cpp']:shutil.copy2(here/name,root/'sources'/name)
for name in ['prepare.py','run.py','run.sbatch','binding.py']:
 if name.endswith('.py'):ast.parse((here/name).read_text())
 shutil.copy2(here/name,root/'scripts'/name)
cases=[]
def case(name,kind,mode,n=4096,steps=64,audit=True,divide=1,nx=96):cases.append({'name':name,'kind':kind,'mode':mode,'particles':n,'steps':steps,'audit':audit,'dt_divisor':divide,'nx':nx,'limit_seconds':180})
for mode in ['baseline','coupled']:case('vacuum-'+mode,'vacuum',mode,n=0,steps=12)
for kind in ['sphere','drift','undulator']:
 for mode in ['baseline','coupled']:case(kind+'-'+mode,kind,mode)
for mode in ['baseline','coupled']:case('sphere-halfdt-'+mode,'sphere',mode,steps=128,divide=2)
for size,n in [('field-heavy',4096),('particle-heavy',131072)]:
 for repeat in range(3):
  for mode in (['baseline','coupled'] if repeat%2==0 else ['coupled','baseline']):case(f'timing-{size}-r{repeat}-{mode}','undulator',mode,n=n,steps=24,audit=False)
library=ws/'build/particle-work-geometry-r2-20261031/libdirect_eb_core.a'
plan={'workspace':str(ws),'build':str(a.build.resolve()),'library':str(library),'cases':cases,'scope':'Independent 3D Cowan/CPML source-coupling prototype; same parameters per baseline/candidate pair; no production switch','represented_electrons':{'gaussian':21300,'sphere':100000000},'field_geometry_cells':[96,96,256],'particle_substeps':5,'seed':20261031,'allocation_physical_CPUs':128,'used_layout':'1 MPI x 8 threads','outer_tolerances':{'state':1e-9,'relative_exchange_work':1e-8,'maximum_iterations':8}}
(root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
frozen={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in list((root/'sources').iterdir())+list((root/'scripts').iterdir())+[root/'plan.json',library]}
# Freeze every unchanged core header/source used to build the accepted library.
for path,digest in json.loads((ws/'results/particle_work_geometry_20261031/build-r2/FROZEN.json').read_text()).items():
 if '/source/aprl/src/' in path or Path(path).name=='CMakeLists.txt':
  assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
  frozen[path]=digest
(root/'FROZEN.json').write_text(json.dumps(frozen,indent=2)+'\n');print(root)
