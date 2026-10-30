#!/usr/bin/env python3
"""Retarget only the stopped undulator case; retain the running free drift."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import yaml

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def main(workspace,root):
    workspace=workspace.resolve();root=root.resolve();source=workspace/'source/aprl';here=Path(__file__).resolve().parent
    previous=workspace/'results/corrected_full_comparison_20261030/r2';old=previous/'undulator_4p5m'
    original=json.loads((old/'run-plan.json').read_text());cfg=yaml.safe_load((old/'run.yaml').read_text())
    # Verify every original frozen input; cancellation changed only outputs.
    frozen=json.loads((old/'frozen.json').read_text());assert all(sha(p)==h for p,h in frozen.items())
    free_frozen=json.loads((previous/'free_drift_0p5m/frozen.json').read_text());assert all(sha(p)==h for p,h in free_frozen.items())
    inp=Path(original['input_file']);stat=inp.stat();assert (stat.st_size,stat.st_mtime_ns)==(original['input_size_bytes'],original['input_mtime_ns'])
    assert sha(inp)==original['input_sha256']
    root.mkdir(parents=True,exist_ok=False);scripts=root/'scripts';scripts.mkdir();d=root/'undulator_4p5m';d.mkdir()
    for name in ['prepare_two_nodes.py','run_two_nodes.sbatch','check_completion.py']:shutil.copy2(here/name,scripts/name)
    for name in ['completion_helpers.py','report_discrete.py','review_discrete.py']:shutil.copy2(previous/'scripts'/name,scripts/name)
    replacements=[(('output','manifest'),'run-manifest.yaml'),(('detectors','directory'),'detectors'),
                  (('energy_ledger','directory'),'ledger'),(('mesh','moving_window','reference_directory'),'reference')]
    card=copy.deepcopy(cfg)
    for keys,target in replacements:
        current=card
        for key in keys[:-1]:current=current[key]
        current[keys[-1]]=str(d/target)
    normalized=copy.deepcopy(card)
    for keys,_ in replacements:
        current=normalized;prior=cfg
        for key in keys[:-1]:current=current[key];prior=prior[key]
        current[keys[-1]]=prior[keys[-1]]
    assert normalized==cfg,'Physics/config changed beyond output paths'
    (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
    geo=copy.deepcopy(original['geometry']);geo['mpi_ranks']=32;geo['minimum_local_field_cells']=geo['active_cells_z']//32
    assert geo['active_cells_z']==4099 and geo['active_cells_z']>=8*32 and geo['minimum_local_field_cells']==128
    save(d/'geometry-plan.json',geo);shutil.copy2(scripts/'run_two_nodes.sbatch',d/'run.sbatch')
    plan=copy.deepcopy(original);plan.update(campaign=str(root),nodes=2,ranks=32,threads_per_rank=8,geometry=geo,
                preflight_pass=str(root/'PRECHECK.json'),previous_single_node_case=str(old),resource_change_only=True)
    save(d/'run-plan.json',plan)
    pre=dict(pass_=True,binary_sha256=original['binary_sha256'],input_sha256=original['input_sha256'],
             card_sha256={str(d/'run.yaml'):sha(d/'run.yaml')},accepted_physics_binary_reused=True,
             global_grid_and_physics_unchanged=True,local_cells_minimum=128,
             new_32_rank_scientific_preflight=False,
             limitation='Existing cross-node acceptance used 16 total ranks; this 32-rank full run must confirm actual binding/layout and compare completed physical outputs.')
    pre['pass']=True;save(root/'PRECHECK.json',pre)
    save(root/'plan.json',dict(root=str(root),workspace=str(workspace),cases={'undulator_4p5m':str(d)},binary=plan['binary'],
        binary_sha256=plan['binary_sha256'],source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
        free_drift_untouched=str(previous/'free_drift_0p5m'),free_drift_job=19098,stopped_undulator_job=19097))
    save(root/'CARD_AUDIT.json',dict(normalized_card_equal=True,only_output_paths_changed=True,
        old_resources={'nodes':1,'ranks':16,'threads':8},new_resources={'nodes':2,'ranks':32,'threads':8},
        active_global_cells=[160,160,4099],minimum_local_z_cells=128,
        free_drift_frozen_manifest_sha256=sha(previous/'free_drift_0p5m/frozen.json')))
    paths=[Path(plan['binary']),workspace/'build/cluster-tools/binding_probe',root/'PRECHECK.json',root/'plan.json',root/'CARD_AUDIT.json',
           old/'run.yaml',old/'run-plan.json',old/'frozen.json',previous/'PRECHECK.json',*list(scripts.iterdir()),
           d/'run.yaml',d/'run-plan.json',d/'geometry-plan.json',d/'run.sbatch']
    identity={str(p):sha(p) for p in paths};save(root/'FROZEN.json',identity);save(d/'frozen.json',identity)
    assert all(sha(p)==h for p,h in free_frozen.items())
    print(json.dumps({'root':str(root),'resources':[2,32,8],'physical_card_unchanged':True,'free_drift_job_untouched':19098},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace,a.root)
