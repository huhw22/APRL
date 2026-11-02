import argparse,copy,json,shutil,hashlib,subprocess
from pathlib import Path
import yaml

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(value,message):
 if not value:raise RuntimeError(message)
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('build_record',type=Path);p.add_argument('output',type=Path);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();b=a.build_record.resolve();identity=json.loads((b/'COMPLETE.json').read_text());assert identity['passed'];program=Path(identity['binary']);assert sha(program)==identity['binary_sha256'];r.mkdir(parents=True,exist_ok=False);scripts=r/'scripts';scripts.mkdir();t=Path(__file__).parent
 for n in ['prepare_flow.py','run_flow.py','flow.sbatch']:shutil.copy2(t/n,scripts/n)
 shutil.copy2(t/'prepare_flow.py',scripts/'prepare.py')
 old=t.parent/'source_coupling_integration'
 for n in ['binding.py','analyze_closed.py']:shutil.copy2(old/n,scripts/n)
 for base,n,target in [(t.parent/'space_charge_physics','launcher.py','launcher.py'),(t.parent/'discrete_energy_ledger','report.py','report_discrete.py'),(t.parent/'discrete_energy_ledger','review_discrete.py','review_discrete.py')]:shutil.copy2(base/n,scripts/target)
 subprocess.run(['/public/home/wfollowm/micromamba/envs/physics/bin/python3','-B','-c','import runpy; import launcher, binding, analyze_closed, report_discrete'],cwd=scripts,check=True)
 base=w/'results/particle_work_geometry_20261031/campaign-r2';cases=[];oldbin=w/'build/source-coupling-exit-r3-20261031/aprl'
 for kind,original in [('drift','native-free-drift-off'),('undulator','native-undulator-off')]:
  src=base/original
  for method in ['baseline','linear','quintic']:
   name=kind+'-'+method;d=r/name;d.mkdir();config=yaml.safe_load((src/'run.yaml').read_text());beam=copy.deepcopy(config['beam'])
   def rewrite(x):
    if isinstance(x,dict):return {k:rewrite(v) for k,v in x.items()}
    if isinstance(x,list):return [rewrite(v) for v in x]
    return x.replace(str(src),str(d)) if isinstance(x,str) else x
   config=rewrite(config);config['beam']=beam;config['mesh']['particle_coupling']='conservative-source'
   if method=='quintic':config['mesh']['longitudinal_shape']='quintic-compensated'
   for key in ['detectors','ledger','trajectory']:(d/key).mkdir()
   (d/'run.yaml').write_text(yaml.safe_dump(config,sort_keys=False));cases.append(dict(name=name,kind=kind,method=method,directory=str(d),binary=str(oldbin if method=='baseline' else program),nodes=1 if kind=='drift' else 2,ranks=16 if kind=='drift' else 32,limit_seconds=600))
 save(r/'plan.json',dict(workspace=str(w),binary=str(program),binary_sha256=sha(program),build_record=str(b),cases=cases,binding_probe=str(w/'build/cluster-tools/binding_probe'),particles=31675))
 paths=[r/'plan.json',program,oldbin,*scripts.iterdir(),*[r/c['name']/'run.yaml' for c in cases]];save(r/'FROZEN.json',{str(p):sha(p) for p in paths});print(r)
if __name__=='__main__':main()
