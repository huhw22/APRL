#!/usr/bin/env python3
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('build',type=Path);a=p.parse_args();root=a.root.resolve();root.mkdir(parents=True);source=Path(__file__).resolve().parents[4]
assert (source/'CMakeLists.txt').exists(),source
snapshot=root/'source';snapshot.mkdir()
files=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=source).decode().split('\0')
for name in files:
 if not name:continue
 q=source/name;target=snapshot/name;target.parent.mkdir(parents=True,exist_ok=True)
 if q.is_symlink():target.symlink_to(q.readlink())
 elif q.is_file():shutil.copy2(q,target)
scripts=root/'scripts';scripts.mkdir()
for name in ['build.py','build.sbatch']:shutil.copy2(Path(__file__).parent/name,scripts/name)
plan={'workspace':'/public/home/wfollowm/APRL_Workspace','source':str(snapshot),'build':str(a.build.resolve()),'live_parent_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),'scope':'Independent full-program conservative-source candidate; no default replacement'}
(root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
paths=[q for q in snapshot.rglob('*') if q.is_file() and not q.is_symlink() and (q.suffix in ['.cpp','.h','.py','.cmake','.txt','.yaml'] or q.name=='CMakeLists.txt')]+list(scripts.iterdir())+[root/'plan.json']
(root/'FROZEN.json').write_text(json.dumps({str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},indent=2)+'\n');print(root)
