"""Freeze a new full-input pair after main-flow preflight, preserving unit charges."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import h5py,yaml

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def differences(a,b,p=''):
 if isinstance(a,dict) and isinstance(b,dict):
  out={}
  for k in set(a)|set(b):out.update(differences(a.get(k),b.get(k),p+'.'+k if p else k))
  return out
 return {} if a==b else {p:[a,b]}
def main():
 p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('build_record',type=Path);p.add_argument('preflight',type=Path);p.add_argument('root',type=Path);a=p.parse_args();w=a.workspace.resolve();build=a.build_record.resolve();flow=a.preflight.resolve();root=a.root.resolve();t=Path(__file__).parent;old=w/'results/source_coupling_full_comparison_20261101';identity=json.loads((build/'COMPLETE.json').read_text());assert identity['passed'];accepted=json.loads((flow/'COMPLETE.json').read_text());assert accepted['program_preflight_pass'];binary=Path(identity['binary']);assert sha(binary)==identity['binary_sha256']
 assert all(sha(q)==h for q,h in json.loads((build/'FROZEN.json').read_text()).items());assert all(sha(q)==h for q,h in json.loads((flow/'FROZEN.json').read_text()).items());root.mkdir(parents=True,exist_ok=False);scripts=root/'scripts';scripts.mkdir()
 for n in ['prepare_full.py']:shutil.copy2(t/n,scripts/n)
 shutil.copy2(t/'prepare_full.py',scripts/'prepare.py')
 originals=t.parent/'source_coupling_full_pair'
 for n in ['run.sbatch','startup_identity.py','check_source_coupling.py']:shutil.copy2(originals/n,scripts/n)
 for src,target in [(t.parent/'corrected_full_pair/check_completion.py','check_completion.py'),(t.parent/'moving_window_production/check_completion.py','completion_helpers.py'),(t.parent/'discrete_energy_ledger/report.py','report_discrete.py'),(t.parent/'discrete_energy_ledger/review_discrete.py','review_discrete.py')]:shutil.copy2(src,scripts/target)
 first=json.loads((old/'undulator_4p5m/run-plan.json').read_text());inp=Path(first['input_file']);st=inp.stat();assert (st.st_size,st.st_mtime_ns)==(first['input_size_bytes'],first['input_mtime_ns']);print('Checking original unit-electron input SHA...',flush=True);input_sha=sha(inp);assert input_sha==first['input_sha256']
 with h5py.File(inp) as h:
  ds=h['particles/records'];assert len(ds)==41732775
  for start in range(0,len(ds),1048576):
   weights=ds.fields('macro_weight')[start:start+1048576];assert (weights==1.).all(),'Input weights differ from unit electrons'
 revision=yaml.safe_load((flow/'drift-linear/run-manifest.yaml').read_text())['source_revision'];cards={};cases={};audits={};paths=[binary,flow/'COMPLETE.json',build/'COMPLETE.json'];compiled={str(q.relative_to(build/'source')):sha(q) for q in (build/'source/src').glob('*') if q.is_file()}
 save(root/'SOURCE_IDENTITY.json',dict(binary=str(binary),binary_sha256=sha(binary),compiled_revision=revision,compiled_source_sha256=compiled,source_head_at_prepare=subprocess.check_output(['git','rev-parse','HEAD'],cwd=w/'source/aprl',text=True).strip(),actual_frozen_build=str(build)))
 for name in ['undulator_4p5m','free_drift_0p5m']:
  src=old/name;prior=json.loads((src/'run-plan.json').read_text());cfg=yaml.safe_load((src/'run.yaml').read_text());base=copy.deepcopy(cfg);d=root/name;d.mkdir()
  assert cfg['beam']['input']['electrons']==41732775 and prior['records']==41732775
  for obj,k,target in [(cfg['output'],'manifest','run-manifest.yaml'),(cfg['mesh']['moving_window'],'reference_directory','reference'),(cfg['detectors'],'directory','detectors'),(cfg['energy_ledger'],'directory','ledger')]:obj[k]=str(d/target)
  cfg['mesh']['longitudinal_shape']='quintic-compensated';assert cfg['mesh']['particle_coupling']=='conservative-source'
  change=differences(base,cfg);assert set(change)=={'output.manifest','mesh.moving_window.reference_directory','detectors.directory','energy_ledger.directory','mesh.longitudinal_shape'},change
  assert cfg['mesh']['moving_window']['lower_guard_cells']==64 and cfg['mesh']['moving_window']['upper_guard_cells']==128
  assert cfg['energy_ledger']['discrete_operator_audit'] and not cfg['particle_retirement']['enabled'] and 'checkpoint' not in cfg
  geo=copy.deepcopy(prior['geometry']);geo['right_collar_cells']+=4;geo['safe_exit_lab_z_m']-=4*geo['dz_m']*geo['gamma_boost'];geo['notes'].append('Quintic compensated support requires four extra trusted-exit collar cells within the existing 128-cell upper guard. Window/grid/time and requested field/particle planes unchanged.')
  (d/'run.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False));save(d/'geometry-plan.json',geo);shutil.copy2(scripts/'run.sbatch',d/'run.sbatch')
  plan=copy.deepcopy(prior);plan.update(campaign=str(root),binary=str(binary),binary_sha256=sha(binary),compiled_revision=revision,source_core='frozen-source:'+str(build/'source'),geometry=geo,preflight_pass=str(root/'PRECHECK.json'),previous_same_geometry_case=str(src),longitudinal_shape='quintic-compensated',particle_coupling='conservative-source',resource_change_only=False)
  save(d/'run-plan.json',plan);cases[name]=str(d);audits[name]=dict(card_changes=change,unit_electron_charge_preserved=True,trusted_exit_collar_change_cells=4);cards[str(d/'run.yaml')]=sha(d/'run.yaml');paths.extend([src/'run.yaml',src/'run-plan.json',d/'run.yaml',d/'run-plan.json',d/'geometry-plan.json',d/'run.sbatch'])
 save(root/'plan.json',dict(workspace=str(w),root=str(root),cases=cases,binary=str(binary),binary_sha256=sha(binary),compiled_revision=revision,preflight=str(flow),physics_accepted=False))
 save(root/'CARD_AUDIT.json',dict(cases=audits,input_sha256=input_sha,input_records=41732775,unit_electron_charge=True,absorption=False,snapshot=False,window_and_sampling_preserved=True))
 save(root/'PRECHECK.json',dict(pass_=True,**{'pass':True},binary_sha256=sha(binary),card_sha256=cards,input_sha256=input_sha,required_tests=78,quintic_current_and_reference_modules_ranks=[1,2,4],main_flow_ranks=[16,32],physics_accepted=False))
 paths.extend([root/'plan.json',root/'CARD_AUDIT.json',root/'PRECHECK.json',root/'SOURCE_IDENTITY.json',*scripts.iterdir()]);frozen={str(q):sha(q) for q in paths if q.is_file()};save(root/'FROZEN.json',frozen)
 for d in cases.values():save(Path(d)/'frozen.json',frozen)
 print(json.dumps(cases,indent=2),flush=True)
if __name__=='__main__':main()
