#!/usr/bin/env python3
import argparse,copy,hashlib,json,shutil
from pathlib import Path
import h5py,yaml
p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--build-review',type=Path);a=p.parse_args()
ws=a.workspace.resolve();root=a.root.resolve();root.mkdir();(root/'scripts').mkdir()
build=json.loads((a.build_review or ws/'results/particle_work_geometry_20261031/build-r1/COMPLETE.json').read_text());assert build['checks_passed']
src=ws/'source/aprl';here=src/'verification/analysis/particle_grid_exchange';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name in ['run_geometry_cases.py','geometry_cases.sbatch','prepare_geometry_cases.py']:shutil.copy2(here/name,root/'scripts'/name)
shutil.copy2(src/'verification/performance/cases/energy_transport_audit/review.py',root/'scripts/transport_review.py')
shutil.copy2(src/'verification/performance/cases/energy_transport_audit/prepare.py',root/'scripts/prepare.py')
for name,rel in [('report_discrete.py','verification/performance/cases/discrete_energy_ledger/report.py'),('review_discrete.py','verification/performance/cases/discrete_energy_ledger/review_discrete.py')]:shutil.copy2(src/rel,root/'scripts'/name)
base=ws/'results/energy_transport_audit_20261029/campaign-r3/science';cases=[]
for layout,ranks,threads in [('native',16,8),('mpi-small',4,8)]:
 for kind in ['undulator','free-drift']:
  original=base/(kind+'-r1-new-audit-on');template=yaml.safe_load((original/'run.yaml').read_text())
  for mode in (['baseline','off','on'] if layout=='native' else ['off','on']):
   name=f'{layout}-{kind}-{mode}';directory=root/name;directory.mkdir();(directory/'geometry').mkdir()
   card=copy.deepcopy(template);card['output']['manifest']=str(directory/'run-manifest.yaml')
   card['mesh']['moving_window'].update(reference_directory=str(directory/'reference'),reference_surface_current=True,inlet_charge_match=True)
   for key,target in [('detectors','detectors'),('trajectory','trajectory'),('energy_ledger','ledger')]:card[key]['directory']=str(directory/target)
   card['energy_ledger']['discrete_operator_audit']=True
   (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
   binary=ws/'build/charge-matched-inlet-20261030/aprl' if mode=='baseline' else Path(build['candidate'])
   cases.append({'name':name,'directory':str(directory),'ranks':ranks,'threads':threads,'mode':mode,'kind':kind,'binary':str(binary),'binary_sha256':sha(binary),'limit_seconds':300})
input_file=Path(template['beam']['input']['file'])
with h5py.File(input_file) as f:n=len(f['particles/records'])
plan={'root':str(root),'cases':cases,'particles':n,'energy_plane_files':{'entry':'audit-entry-energy.h5','exit':'audit-exit-energy.h5'},'purpose':'Read-only work-path diagnostic; physics off/on and accepted/candidate comparison; no method change','build_review':str(a.build_review) if a.build_review else None,'input':str(input_file),'input_sha256':sha(input_file)}
(root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
frozen=[root/'plan.json',input_file,*[Path(c['binary']) for c in cases],*[Path(c['directory'])/'run.yaml' for c in cases],*[q for q in (root/'scripts').iterdir() if q.is_file()]]
(root/'FROZEN.json').write_text(json.dumps({str(q):sha(q) for q in frozen},indent=2)+'\n')
print(root)
