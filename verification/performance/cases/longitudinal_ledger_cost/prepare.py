"""Same-binary, same-particle diagnostic cost decomposition; no solver changes."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import yaml

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def require(value,message):
 if not value:raise RuntimeError(message)
def main():
 p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('output',type=Path);a=p.parse_args();w=a.workspace.resolve();r=a.output.resolve();t=Path(__file__).parent;source=w/'results/longitudinal_main_flow_20261103/flow-r2/drift-quintic';assert json.loads((source/'RESULT.json').read_text())['program_integrity_pass'];record=json.loads((w/'results/longitudinal_main_flow_20261103/build-r1/COMPLETE.json').read_text());binary=Path(record['binary']);assert sha(binary)==record['binary_sha256'];r.mkdir(parents=True,exist_ok=False);scripts=r/'scripts';scripts.mkdir()
 for n in ['prepare.py','run.py','run.sbatch']:shutil.copy2(t/n,scripts/n)
 for parent,n in [(t.parent/'space_charge_physics','launcher.py'),(t.parent/'source_coupling_integration','binding.py')]:shutil.copy2(parent/n,scripts/n)
 subprocess.run(['/public/home/wfollowm/micromamba/envs/physics/bin/python3','-B','-c','import launcher,binding'],cwd=scripts,check=True)
 template=yaml.safe_load((source/'run.yaml').read_text());cases=[]
 for repeat,order in enumerate([['off','basic','strict','all'],['all','strict','basic','off']]):
  for mode in order:
   name=f'r{repeat}-{mode}';d=r/name;d.mkdir();beam=copy.deepcopy(template['beam'])
   def relocate(x):
    if isinstance(x,dict):return {k:relocate(v) for k,v in x.items()}
    if isinstance(x,list):return [relocate(v) for v in x]
    return x.replace(str(source),str(d)) if isinstance(x,str) else x
   cfg=relocate(template);cfg['beam']=beam;cfg['energy_ledger']['enabled']=mode!='off';cfg['energy_ledger']['discrete_operator_audit']=mode in ['strict','all'];cfg['energy_ledger']['sample_interval_steps']=10
   cfg['detectors']['energy_planes']=copy.deepcopy(template['detectors']['energy_planes']) if mode=='all' else []
   for plane in cfg['detectors']['energy_planes']:plane['buffer_samples']=64
   # Raw field and particle planes, window, stepping and source solve remain.
   assert cfg['mesh']['longitudinal_shape']=='quintic-compensated'
   for key in ['detectors','ledger','trajectory']:(d/key).mkdir()
   (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));cases.append(dict(name=name,repeat=repeat,mode=mode,directory=str(d),limit_seconds=300))
 plan=dict(workspace=str(w),binary=str(binary),binary_sha256=sha(binary),source_case=str(source),ranks=16,threads=8,nodes=1,cases=cases,normal_ledger_sample_interval_steps=10,strict_audit_every_step=True,energy_planes=2,profile_enabled=False,production_jobs_unchanged=[19261,19262]);save(r/'PLAN.json',plan)
 paths=[r/'PLAN.json',binary,w/'build/cluster-tools/binding_probe',*scripts.iterdir(),*[r/c['name']/'run.yaml' for c in cases],Path(beam['input']['file'])];save(r/'FROZEN.json',{str(q):sha(q) for q in paths});print(r)
if __name__=='__main__':main()
