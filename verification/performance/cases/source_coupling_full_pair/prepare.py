#!/usr/bin/env python3
"""Freeze the accepted conservative-source full-input comparison pair."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import h5py,yaml

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def differences(a,b,prefix=''):
    if isinstance(a,dict) and isinstance(b,dict):
        out={}
        for k in set(a)|set(b):out.update(differences(a.get(k),b.get(k),prefix+'.'+k if prefix else k))
        return out
    return {} if a==b else {prefix:[a,b]}

def main(ws,root):
    ws=ws.resolve();root=root.resolve();here=Path(__file__).resolve().parent
    src=ws/'source/aprl';accepted=ws/'results/source_coupling_exit_20261031'
    identity=json.loads((accepted/'SOURCE_FINAL.json').read_text());binary=Path(identity['binary'])
    expected='865ebef5b79b6b8670d40c81080f15f916ab1938033afcabd050bcd3c8198752'
    assert sha(binary)==identity['binary_sha256']==expected
    final=json.loads((accepted/'FINAL.json').read_text());assert final['tail_handling_accepted_for_these_cases'] and final['physics_cases_complete']
    assert final['source_commit']=='88a5a8d9fbd2af41e6bd022a19d36304418ac364'
    assert json.loads((accepted/'two-nodes-r2/COMPLETE.json').read_text())['program_preflight_pass']
    assert '100% tests passed, 0 tests failed out of 78' in (accepted/'build-r3/required.log').read_text()
    frozen_src=accepted/'build-r3/source/src';compiled={}
    for p in sorted(frozen_src.rglob('*')):
        if p.is_file():
            rel=p.relative_to(frozen_src);assert sha(p)==sha(src/'src'/rel),str(rel);compiled[str(rel)]=sha(p)
    oldroot=ws/'results/corrected_full_comparison_20261030'
    originals={'undulator_4p5m':oldroot/'r3_2nodes/undulator_4p5m','free_drift_0p5m':oldroot/'r2/free_drift_0p5m'}
    first=json.loads((originals['undulator_4p5m']/'run-plan.json').read_text());inp=Path(first['input_file']);st=inp.stat()
    assert (st.st_size,st.st_mtime_ns)==(first['input_size_bytes'],first['input_mtime_ns'])
    print('Hashing unchanged original particle input...',flush=True);input_sha=sha(inp);assert input_sha==first['input_sha256']
    with h5py.File(inp) as h:assert len(h['particles/records'])==41732775
    root.mkdir(parents=True,exist_ok=False);scripts=root/'scripts';scripts.mkdir()
    for p in here.iterdir():
        if p.is_file():shutil.copy2(p,scripts/p.name)
    copies={
        here.parent/'corrected_full_pair/check_completion.py':'check_completion.py',
        here.parent/'moving_window_production/check_completion.py':'completion_helpers.py',
        here.parent/'discrete_energy_ledger/report.py':'report_discrete.py',
        here.parent/'discrete_energy_ledger/review_discrete.py':'review_discrete.py',
        here.parent/'corrected_full_pair/check_id_repair.py':'check_id_repair.py'}
    for p,name in copies.items():shutil.copy2(p,scripts/name)
    save(root/'SOURCE_IDENTITY.json',dict(core_commit=final['source_commit'],compiled_revision=identity['compiled_revision'],compiled_src_sha256=compiled,binary_sha256=expected,diagnostic_force_split_binary_used=False))
    cases={};audits={};cards={};paths=[binary,ws/'build/cluster-tools/binding_probe',root/'SOURCE_IDENTITY.json']
    for name,old in originals.items():
        prior=json.loads((old/'run-plan.json').read_text());base=yaml.safe_load((old/'run.yaml').read_text());d=root/name;d.mkdir()
        assert prior['input_sha256']==input_sha and prior['records']==41732775
        assert prior['nodes']==(2 if name=='undulator_4p5m' else 1)
        for p,h in json.loads((old/'frozen.json').read_text()).items():
            if p==str(old/'run.yaml'):assert sha(p)==h
        card=copy.deepcopy(base)
        for obj,k,target in [(card['output'],'manifest','run-manifest.yaml'),(card['mesh']['moving_window'],'reference_directory','reference'),(card['detectors'],'directory','detectors'),(card['energy_ledger'],'directory','ledger')]:obj[k]=str(d/target)
        card['mesh']['particle_coupling']='conservative-source'
        change=differences(base,card);allowed={'output.manifest','mesh.moving_window.reference_directory','mesh.particle_coupling','detectors.directory','energy_ledger.directory'}
        assert set(change)==allowed,change
        assert card['mesh']['moving_window']['lower_guard_cells']==64 and card['mesh']['moving_window']['upper_guard_cells']==128
        assert card['mesh']['moving_window']['reference_surface_current'] and card['mesh']['moving_window']['inlet_charge_match']
        assert card['energy_ledger']['enabled'] and card['energy_ledger']['discrete_operator_audit']
        assert not card['particle_retirement']['enabled'] and 'checkpoint' not in card
        assert len(card['detectors']['particle_planes'])==2
        geo=copy.deepcopy(prior['geometry']);assert geo['active_cells_z'] in [666,4099]
        assert geo['active_cells_z']//prior['ranks']>=8
        (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False));save(d/'geometry-plan.json',geo);shutil.copy2(scripts/'run.sbatch',d/'run.sbatch')
        p=copy.deepcopy(prior);p.update(campaign=str(root),binary=str(binary),binary_sha256=expected,compiled_revision=identity['compiled_revision'],source_core=final['source_commit'],preflight_pass=str(root/'PRECHECK.json'),previous_same_geometry_case=str(old),particle_coupling='conservative-source',resource_change_only=False)
        save(d/'run-plan.json',p);cases[name]=str(d);audits[name]=change;cards[str(d/'run.yaml')]=sha(d/'run.yaml')
        paths +=[old/'run.yaml',old/'run-plan.json',d/'run.yaml',d/'geometry-plan.json',d/'run-plan.json',d/'run.sbatch']
    evidence=[accepted/'SOURCE_FINAL.json',accepted/'FINAL.json',accepted/'two-nodes-r2/COMPLETE.json',accepted/'build-r3/required.log',ws/'results/cowan_transverse_attribution_20261101/force-split-late-r1/FINAL.json']
    save(root/'plan.json',dict(root=str(root),workspace=str(ws),cases=cases,binary=str(binary),binary_sha256=expected,core_commit=final['source_commit'],source_head_at_prepare=subprocess.check_output(['git','rev-parse','HEAD'],cwd=src,text=True).strip(),compiled_revision=identity['compiled_revision'],physics_accepted=False))
    save(root/'CARD_AUDIT.json',dict(changes=audits,geometry_and_input_and_sampling_and_stop_unchanged=True,only_physics_selector_changed='mesh.particle_coupling',full_input_sha256=input_sha))
    pre=dict(binary_sha256=expected,input_sha256=input_sha,card_sha256=cards,accepted_binary_reused=True,required_tests_reused=78,existing_cross_node_total_ranks=16,new_32_rank_scientific_preflight=False,actual_binding_required=True,compiled_source_files_matched=len(compiled))
    pre['pass']=True;save(root/'PRECHECK.json',pre)
    paths += evidence+list(scripts.glob('*'))+[root/'plan.json',root/'CARD_AUDIT.json',root/'PRECHECK.json']
    frozen={str(p):sha(p) for p in paths if p.is_file()};save(root/'FROZEN.json',frozen)
    for d in cases.values():save(Path(d)/'frozen.json',frozen)
    print(json.dumps(dict(root=str(root),cases=cases,precheck=True)),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace,a.root)
