#!/usr/bin/env python3
"""Freeze the requested full-input pair with accepted inlet/boundary fixes."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import h5py,yaml

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(2**20),b''):h.update(chunk)
    return h.hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def diff(a,b,prefix=''):
    out={}
    if isinstance(a,dict) and isinstance(b,dict):
        for key in set(a)|set(b):out.update(diff(a.get(key),b.get(key),prefix+'.'+key if prefix else key))
    elif a!=b:out[prefix]=[a,b]
    return out

def main(workspace,root):
    workspace=workspace.resolve();root=root.resolve();source=workspace/'source/aprl';here=Path(__file__).resolve().parent
    oldroot=workspace/'results/moving_window_comparison_20261029'
    bench=workspace/'results/energy_ledger_test_record_20261030/physical_benchmarks'
    accepted=bench/'charge_matched_inlet_20261030';build=workspace/'build/charge-matched-inlet-20261030';binary=build/'aprl'
    binary_sha=sha(binary);assert binary_sha=='79d643bc8ade4db4d9e4488a367cb1bf5ef74e7c602dc3a2094462709aed9435'
    assert json.loads((accepted/'preflight/COMPLETION.json').read_text())['checks_pass']
    assert '100% tests passed, 0 tests failed out of 78' in (accepted/'preflight/required.log').read_text()
    known=json.loads((bench/'charge_matched_source_audit_20261030_r2/SOURCE_FINAL.json').read_text())['audit_cpp_and_headers']
    changed=[p for p,h in known.items() if sha(source/p)!=h]
    assert changed==['src/particle_momentum_probe.h'],changed
    compiled=yaml.safe_load((accepted/'pilot_500/run-manifest.yaml').read_text())['source_revision']
    root.mkdir(parents=True,exist_ok=False);scripts=root/'scripts';scripts.mkdir()
    for p in here.iterdir():
        if p.is_file():shutil.copy2(p,scripts/p.name)
    shutil.copy2(here.parent/'moving_window_production/check_completion.py',scripts/'completion_helpers.py')
    for a,b in [('report.py','report_discrete.py'),('review_discrete.py','review_discrete.py')]:
        shutil.copy2(here.parent/'discrete_energy_ledger'/a,scripts/b)
    first=json.loads((oldroot/'undulator_4p5m/run-plan.json').read_text());inp=Path(first['input_file']);st=inp.stat()
    assert (st.st_size,st.st_mtime_ns)==(first['input_size_bytes'],first['input_mtime_ns'])
    input_sha=sha(inp);assert input_sha==first['input_sha256']
    oldpass=json.loads((oldroot/'preflight/PASS.json').read_text());assert oldpass['pass'] and oldpass['input_sha256']==input_sha
    evidence=[accepted/'preflight/COMPLETION.json',accepted/'preflight/required.log',accepted/'SOURCE_BUILD.json',
              bench/'charge_matched_source_audit_20261030_r2/SOURCE_FINAL.json',oldroot/'preflight/PASS.json']
    plan=dict(format='aprl-corrected-full-pair-v1',workspace=str(workspace),root=str(root),binary=str(binary),binary_sha256=binary_sha,
              compiled_revision=compiled,source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
              source_core='81f5f0b8e7fb4b174017c617fd3c6e26a59059ee',cases={},acceptance_reuse=True,
              note='Accepted optimized moving-window physics binary; latest commits only read-only probes and Python reports. No new build.')
    audits={};casefiles=[]
    allowed={'output.manifest','mesh.moving_window.reference_directory','mesh.moving_window.reference_surface_current',
             'mesh.moving_window.inlet_charge_match','detectors.directory','detectors.particle_planes','detectors.energy_planes',
             'energy_ledger.directory','energy_ledger.discrete_operator_audit'}
    for name,z in [('undulator_4p5m',4.5),('free_drift_0p5m',.5)]:
        old=oldroot/name;prior=json.loads((old/'run-plan.json').read_text());base=yaml.safe_load((old/'run.yaml').read_text());geo=prior['geometry']
        assert prior['input_sha256']==input_sha and prior['records']==41732775 and base['detectors']['field_planes'][0]['z']==z
        assert oldpass['card_sha256'][str(old/'run.yaml')]==sha(old/'run.yaml')
        d=root/name;d.mkdir();card=copy.deepcopy(base)
        card['output']['manifest']=str(d/'run-manifest.yaml')
        card['mesh']['moving_window'].update(reference_directory=str(d/'reference'),reference_surface_current=True,inlet_charge_match=True)
        card['detectors']['directory']=str(d/'detectors')
        card['detectors']['particle_planes']=[dict(name='entry-particles',z=-.05,buffer_records=32768,compression=0),
                                                dict(name='exit-particles',z=z,buffer_records=32768,compression=0)]
        card['detectors']['energy_planes']=[dict(name='entry-transport',z=-.05,rhythm=geo['dt_box_s']/geo['gamma_boost'],buffer_samples=64,compression=0),
                                              dict(name='exit-transport',z=z,rhythm=geo['dt_box_s']/geo['gamma_boost'],buffer_samples=64,compression=0)]
        card['energy_ledger'].update(directory=str(d/'ledger'),discrete_operator_audit=True)
        changes=diff(base,card);assert set(changes)<=allowed,changes
        assert card['beam']==base['beam'] and card['sources']==base['sources'] and card['stop']==base['stop']
        assert card['boundary']==base['boundary'] and card['mesh']['cell_size']==base['mesh']['cell_size']
        assert card['detectors']['field_planes']==base['detectors']['field_planes']
        assert not card.get('checkpoint') and not card['particle_retirement']['enabled']
        # All diagnostic planes lie inside the compact window at t_box=0, so
        # they do not force the fallback that scans every pending particle.
        distances={}
        for plane in card['detectors']['particle_planes']+card['detectors']['energy_planes']:
            box=(plane['z']-geo['release_lab_z_m'])/geo['gamma_boost']+geo['translation_z_m']
            distance=(box-geo['active_origin_z_m'])/geo['dz_m'];assert 0<=distance<=geo['active_cells_z']
            distances[plane['name']]=distance
        (d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False));save(d/'geometry-plan.json',geo)
        shutil.copy2(scripts/'run.sbatch',d/'run.sbatch')
        run=dict(format=plan['format'],case=name,campaign=str(root),binary=str(binary),binary_sha256=binary_sha,
                 compiled_revision=compiled,source_core=plan['source_core'],nodes=1,ranks=16,threads_per_rank=8,
                 input_file=str(inp),input_sha256=input_sha,input_size_bytes=st.st_size,input_mtime_ns=st.st_mtime_ns,
                 records=prior['records'],field_plane_z_m=z,geometry=geo,original_case=str(old),
                 preflight_pass=str(root/'PRECHECK.json'),expected_pending_queue=True)
        save(d/'run-plan.json',run)
        audits[name]=dict(changes=changes,diagnostic_distance_from_inlet_cells=distances,
                         original_card_sha256=sha(old/'run.yaml'),new_card_sha256=sha(d/'run.yaml'))
        plan['cases'][name]=str(d);casefiles.extend([d/'run.yaml',d/'run-plan.json',d/'geometry-plan.json',d/'run.sbatch'])
        evidence.extend([old/'run.yaml',old/'run-plan.json',old/'geometry-plan.json'])
    save(root/'plan.json',plan);save(root/'CARD_AUDIT.json',audits)
    pre=dict(pass_=True,binary_sha256=binary_sha,input_sha256=input_sha,records=41732775,
             new_simulation_preflight_run=False,acceptance_reuse=True,required_tests=78,
             prior_cpp_differences=changed,cpp_difference_scope='Only disabled read-only momentum probe modes; no physics kernel change.',
             card_sha256={str(Path(d)/'run.yaml'):sha(Path(d)/'run.yaml') for d in plan['cases'].values()},
             scope='Existing same-binary acceptance plus unchanged completed full-input geometry and fresh whole-input hash; full corrected physics remains to compare.')
    pre['pass']=True;save(root/'PRECHECK.json',pre)
    frozen=[binary,Path(workspace/'build/cluster-tools/binding_probe'),root/'plan.json',root/'CARD_AUDIT.json',root/'PRECHECK.json',*evidence,*casefiles,*scripts.iterdir()]
    # The input was freshly hashed once above; startup/end verify its frozen
    # size/mtime, avoiding two extra concurrent multi-gigabyte reads.
    identity={str(p):sha(p) for p in frozen if p.is_file()}
    for d in plan['cases'].values():save(Path(d)/'frozen.json',identity)
    save(root/'FROZEN.json',identity)
    print(json.dumps({'root':str(root),'input_sha256':input_sha,'compiled_revision':compiled,'cases':audits},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);p.add_argument('root',type=Path);a=p.parse_args();main(a.workspace,a.root)
