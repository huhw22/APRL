#!/usr/bin/env python3
import argparse,copy,hashlib,json,shutil,subprocess,sys
from pathlib import Path
import yaml
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def require(v,msg):
 if not v:raise RuntimeError(msg)
def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('build_record',type=Path);p.add_argument('--nodes',type=int,choices=[1,2],default=1);p.add_argument('--comparison',type=Path);p.add_argument('--resume',type=Path);a=p.parse_args();root=a.root.resolve();root.mkdir(parents=True);scripts=root/'scripts';scripts.mkdir();here=Path(__file__).parent
 build_record=a.build_record.resolve();ws=Path(json.loads((build_record/'plan.json').read_text())['workspace']);built=json.loads((build_record/'COMPLETE.json').read_text());assert built['passed'];binary=Path(built['binary']);assert sha(binary)==built['binary_sha256']
 for n in ['prepare.py','run.py','run.sbatch','analyze_closed.py','binding.py','select_rank_mode.py']:shutil.copy2(here/n,scripts/n)
 for n in ['launcher.py','model.py']:shutil.copy2(here.parent/'space_charge_physics'/n,scripts/n)
 for src,name in [(here.parent/'space_charge_physics/analyze.py','analyze_particles.py'),(here.parent/'discrete_energy_ledger/report.py','report_discrete.py'),(here.parent/'discrete_energy_ledger/review_discrete.py','review_discrete.py'),(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py','relative_ledger.py')]:shutil.copy2(src,scripts/name)
 subprocess.run([sys.executable,'-B','-c','import report_discrete, launcher, analyze_closed, analyze_particles'],cwd=scripts,check=True)
 prior=ws/'results/particle_work_geometry_20261031/campaign-r2';cases={}
 for kind in ['free-drift','undulator']:
  for ranks,layout in ([(16,'native'),(4,'mpi-small')] if a.nodes==1 else [(16,'native')]):
   for mode in (['off','on'] if ranks==16 and a.nodes==1 else ['on']):
    name=f'{kind}-{ranks}-{mode}';old=prior/f'{layout}-{kind}-off';d=root/name
    if a.resume and (a.resume/name/'RESULT.json').exists():
     reused=a.resume.resolve()/name;assert json.loads((reused/'RESULT.json').read_text())['program_integrity_pass'];d.symlink_to(reused,target_is_directory=True);cases[name]={'name':name,'kind':kind,'ranks':ranks,'threads':8,'mode':mode,'directory':str(d),'reference':str(old),'reuse_closed':str(reused)};continue
    d.mkdir();cfg=yaml.safe_load((old/'run.yaml').read_text());input_config=copy.deepcopy(cfg['beam'])
    def rewrite(x):
     if isinstance(x,dict):return {k:rewrite(v) for k,v in x.items()}
     if isinstance(x,list):return [rewrite(v) for v in x]
     return x.replace(str(old),str(d)) if isinstance(x,str) else x
    cfg=rewrite(cfg);cfg['beam']=input_config;cfg['mesh']['particle_coupling']='conservative-source' if mode=='on' else 'boris';(d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));cases[name]={'name':name,'kind':kind,'ranks':ranks,'threads':8,'mode':mode,'directory':str(d),'reference':str(old),'limit_seconds':600}
    if a.comparison:cases[name]['same_method_reference']=str(a.comparison.resolve()/name)
 plan={'workspace':str(ws),'build_record':str(build_record),'build':str(binary.parent),'binary':str(binary),'binary_sha256':sha(binary),'binding_probe':str(ws/'build/cluster-tools/binding_probe'),'cases':cases,'nodes':a.nodes,'particles':json.loads((prior/'plan.json').read_text())['particles'],'scope':'Full-program held release/MPI/window/coupling preflight; full original sphere suite follows only after these gates'}
 save(root/'plan.json',plan);paths=[root/'plan.json',binary,*scripts.iterdir()]
 if 'particle_coupling must agree' in (build_record/'source/src/simulation.cpp').read_text():
  negative=root/'inconsistent-mode';negative.mkdir();configs=[]
  for mode in ['boris','conservative-source']:
   cfg=yaml.safe_load((prior/'native-free-drift-off/run.yaml').read_text());old=str(prior/'native-free-drift-off')
   def relocate(x):
    if isinstance(x,dict):return {k:relocate(v) for k,v in x.items()}
    if isinstance(x,list):return [relocate(v) for v in x]
    return x.replace(old,str(negative/'output')) if isinstance(x,str) else x
   beam=copy.deepcopy(cfg['beam']);cfg=relocate(cfg);cfg['beam']=beam;cfg['mesh']['particle_coupling']=mode;p=negative/(mode+'.yaml');p.write_text(yaml.safe_dump(cfg,sort_keys=False));configs.append(str(p));paths.append(p)
  plan['rank_mode_negative_configs']=configs;save(root/'plan.json',plan)
 for q in cases.values():
  p=Path(q['directory'])/'run.yaml';paths.extend([p,Path(yaml.safe_load(p.read_text())['beam']['input']['file'])])
  if q.get('reuse_closed'):paths.append(Path(q['reuse_closed'])/'RESULT.json')
 save(root/'FROZEN.json',{str(p):sha(p) for p in paths});print(root)
if __name__=='__main__':main()
