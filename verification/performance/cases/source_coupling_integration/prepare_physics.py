#!/usr/bin/env python3
"""Re-use the original physical cards and immutable worldline inputs."""
import argparse,copy,json,shutil,subprocess,sys
from pathlib import Path
import yaml
from prepare import sha,save
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('build_record',type=Path);p.add_argument('preflight',type=Path);a=p.parse_args();root=a.root.resolve();root.mkdir(parents=True);here=Path(__file__).resolve().parent;scripts=root/'scripts';scripts.mkdir()
brec=a.build_record.resolve();built=json.loads((brec/'COMPLETE.json').read_text());assert built['passed'];ws=Path(json.loads((brec/'plan.json').read_text())['workspace']);binary=Path(built['binary']);assert sha(binary)==built['binary_sha256'];preflight=a.preflight.resolve();assert json.loads((preflight/'COMPLETE.json').read_text())['program_preflight_pass']
for n in ['prepare.py','run_physics.py','physics.sbatch','analyze_closed.py','binding.py']:shutil.copy2(here/n,scripts/n)
for n in ['launcher.py','model.py']:shutil.copy2(here.parent/'space_charge_physics'/n,scripts/n)
for src,dst in [(here.parent/'space_charge_physics/analyze.py','analyze_particles.py'),(here.parent/'discrete_energy_ledger/report.py','report_discrete.py'),(here.parent/'discrete_energy_ledger/review_discrete.py','review_discrete.py'),(here.parents[2]/'analysis/relative_energy_ledger/relative_ledger.py','relative_ledger.py')]:shutil.copy2(src,scripts/dst)
subprocess.run([sys.executable,'-B','-c','import report_discrete, launcher, analyze_closed, analyze_particles'],cwd=scripts,check=True)
bench=ws/'results/energy_ledger_test_record_20261030/physical_benchmarks';old=bench/'charge_matched_inlet_20261030';scaling=bench/'space_charge_scaling_20261030';cases={}
for name,reference in [(x,old/(x+'_full')) for x in ['uniform','converging','gaussian']]+[(x,scaling/x) for x in ['q0p1','q0p01','q0p001']]:
 d=root/name;d.mkdir();(d/'momentum').mkdir();cfg=yaml.safe_load((reference/'run.yaml').read_text());beam=copy.deepcopy(cfg['beam'])
 def rewrite(x):
  if isinstance(x,dict):return {k:rewrite(v) for k,v in x.items()}
  if isinstance(x,list):return [rewrite(v) for v in x]
  return x.replace(str(reference),str(d)) if isinstance(x,str) else x
 cfg=rewrite(cfg);cfg['beam']=beam;cfg['mesh']['particle_coupling']='conservative-source';(d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));shutil.copy2(reference/'case.json',d/'case.json');(d/'reference_input').symlink_to(reference/'reference_input',target_is_directory=True)
 audit=bench/'charge_matched_source_audit_20261030_r2'/(name+'_on')/'FIRST_SOURCE_STATES.npy'
 cases[name]={'directory':str(d),'original':str(reference),'steps':json.loads((d/'case.json').read_text())['geometry']['expected_stop_step'],'timeout_s':2700,'first_source_reference':str(audit) if audit.exists() else None}
plan={'workspace':str(ws),'build':str(binary.parent),'binary':str(binary),'binary_sha256':sha(binary),'binding_probe':str(ws/'build/cluster-tools/binding_probe'),'preflight':str(preflight),'cases':cases,'jobs':{'uniform_and_scaling':['uniform','q0p1','q0p01','q0p001'],'converging':['converging'],'gaussian':['gaussian']},'nominal_dt_s':1.2543262415359034e-6/299792458.,'scope':'Original full uniform/converging/Gaussian and low-charge physical cards, only coupling method and output locations changed; no theory retuning'}
save(root/'plan.json',plan);files=[root/'plan.json',binary,*scripts.iterdir()]
for q in cases.values():
 d=Path(q['directory']);files.extend([d/'run.yaml',d/'case.json',d/'reference_input/expected_events.npz',Path(yaml.safe_load((d/'run.yaml').read_text())['beam']['input']['file'])])
save(root/'FROZEN.json',{str(p):sha(p) for p in files});print(root)
