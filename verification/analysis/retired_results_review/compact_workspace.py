"""Compact historical raw outputs under the user's <=500GB workspace goal."""
import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
import sys

from execute_approved import current_frozen_checks,fingerprint,match,unlink_verified

W=Path('/public/home/wfollowm/APRL_Workspace')
REPORT=W/'results/workspace_compaction_20261115_r1'
CURRENT=W/'results/fel_independent_w355_20261114_r1'
TARGET=500_000_000_000
RAW={'.h5','.hdf5','.bin','.raw','.npy','.npz'}

# Preserve current science and the minimum complete canonical earlier controls.
SCIENCE=[
 'physics_extension_20261112_r1',
 'accepted_physics_transfer_20261112',
 'synchrotron_independent_electron_100k_seed20261111_20261110_r1',
 'synchrotron_independent_electron_100k_seed20261112_20261110_r1',
 'synchrotron_sample_matched_lw_20261110_r2',
 'puffin_aprl_source_propagation_20261113_r1',
 'sparc_single_seed_gap_20261113_r1',
 'synchrotron_full_ledger_20261112_r1',
 'energy_ledger_test_record_20261030/evidence',
 'energy_ledger_test_record_20261030/physical_benchmarks/space_charge_pic_20261030_r2/uniform',
 'energy_ledger_test_record_20261030/physical_benchmarks/space_charge_pic_20261030_r2/converging',
 'energy_ledger_test_record_20261030/physical_benchmarks/space_charge_pic_gaussian_20261030_r4/gaussian',
 'energy_ledger_test_record_20261030/physical_benchmarks/space_charge_scaling_20261030',
 'retired_results_review_20261114_r1',
]
for case in ['w64-a','w512-a','w512-b']:
    for sub in ['detectors','ledger']:SCIENCE.append('sparc_random_startup_20261112_r2/'+case+'/'+sub)
SIX_LEDGER='sparc_six_full_ledger_20261112_r3'
for case in ['energy-binomial3','energy-momentum']:
    SCIENCE += [SIX_LEDGER+'/'+case+'/'+x for x in ['ledger','lab-audit','worldtube']]


def save(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def under(path,root):return path==root or path.startswith(root+'/')


def input_paths(rows):
    exact=set()
    def add(value):
        if isinstance(value,str) and value.startswith('/'):
            p=Path(value)
            if p.is_relative_to(W) and p.exists():exact.add(str(p));exact.add(str(p.resolve()))
    # Beam.input.file is the simulator's actual particle-input contract.
    import yaml
    for path,*unused in rows:
        p=Path(path)
        if p.name=='run.yaml' and '/source/' not in path and '/scripts/' not in path:
            try:
                value=yaml.safe_load(p.read_text())
                source=value.get('beam',{}).get('input',{})
                if isinstance(source,dict):add(source.get('file'))
            except (OSError,ValueError,AttributeError,yaml.YAMLError):pass
        if p.name in ['INPUT_READY.json','INPUT_VALIDATED.json','MASTER_READY.json'] and '/source/' not in path and '/scripts/' not in path:
            try:
                value=read(p)
                for key in ['input','input_file','master','master_file','source_input','source_file']:add(value.get(key))
            except (OSError,ValueError,AttributeError):pass
    index=W/'inputs/particle_inputs/INDEX.json'
    if index.exists():
        def visit(value):
            if isinstance(value,dict):
                for v in value.values():visit(v)
            elif isinstance(value,list):
                for v in value:visit(v)
            else:add(value)
        visit(read(index))
    return exact


def plan():
    assert read(REPORT/'AUTHORIZATION.json')['execute_old_output_cleanup']
    with sqlite3.connect(REPORT/'inventory.sqlite') as db:
        rows=list(db.execute("SELECT path,dev,ino,nlink,size,allocated,mtime_ns FROM files WHERE kind='file'"))
    frozen=current_frozen_checks(W,REPORT)
    save(REPORT/'CURRENT_FROZEN_BEFORE.json',frozen)
    exact=input_paths(rows)|set(frozen['files'])
    # Keep final scalar substep particles of the independent 60M benchmark,
    # its final field/energy ports and dense ledgers; discard intermediate planes.
    full=W/'results/binomial3_fel_20261106/full-r1/one4one'
    exact.update(str(full/'detectors'/n) for n in ['p37-particles-substeps.h5','p37.h5','p37-transport-energy.h5'])
    prefixes=[str(CURRENT),str(REPORT),*[str(W/n) for n in ['source','build','inputs','tests','logs','reference']],
        *[str(W/'results'/n) for n in SCIENCE],str(full/'ledger')]
    keep_inodes=set();decisions={};all_inodes={};all_paths=collections.defaultdict(list)
    for path,dev,ino,nlink,size,allocated,mtime in rows:
        key=(dev,ino);all_paths[key].append(path)
        all_inodes[key]=max(allocated,all_inodes.get(key,0));p=Path(path);parts=p.relative_to(W).parts
        reason=None
        if not parts or parts[0]!='results':reason='non_results_code_inputs_and_metadata'
        elif any(under(path,q) for q in prefixes):reason='current_or_canonical_science'
        elif path in exact:reason='declared_original_input_or_frozen_dependency'
        elif any(x in ['source','src','scripts','bin','inputs','frozen_source'] or x.startswith('build') for x in parts):reason='source_program_or_input_tree'
        elif p.name in ['particles.h5','master.h5','beam.h5','input.h5','common.par.h5'] or p.name.startswith('input') or p.name.endswith('_electrons_0.h5') or (p.name.endswith('.par.h5') and 'dump' not in p.name):reason='original_or_native_particle_input'
        elif p.suffix.lower() not in RAW:reason='plots_tables_cards_logs_reports_or_other_nonraw'
        elif size<8*1024**2:reason='small_raw_fixture_or_compact_evidence'
        elif p.suffix.lower()=='.npz' and size<32*1024**2:reason='compact_analysis_cache'
        elif p.suffix.lower()=='.npy' and size<16*1024**2:reason='compact_analysis_cache'
        elif os.lstat(path).st_mode & 0o111:reason='executable_program'
        if reason:keep_inodes.add(key);decisions[path]=reason
    candidates=[];retained=[]
    for path,dev,ino,nlink,size,allocated,mtime in rows:
        key=(dev,ino);reason=decisions.get(path)
        if reason is None and key in keep_inodes:reason='hardlink_to_preserved_data'
        if reason is None and len(all_paths[key])!=nlink:reason='unobserved_external_hardlinks'
        item=dict(path=path,dev=dev,inode=ino,nlink=nlink,size_bytes=size,allocated_bytes=allocated,mtime_ns=mtime)
        if reason:
            if size>=8*1024**2:retained.append(dict(reason=reason,**item))
        else:candidates.append(item)
    by_inode=collections.defaultdict(list)
    for row in candidates:by_inode[row['dev'],row['inode']].append(row)
    reclaimed=0
    for values in by_inode.values():
        assert len(values)==values[0]['nlink'];reclaimed+=values[0]['allocated_bytes']
    total=sum(all_inodes.values());expected=total-reclaimed
    groups=collections.defaultdict(lambda:dict(files=0,bytes=0,inodes=set()))
    for row in candidates:
        name=Path(row['path']).relative_to(W/'results').parts[0];g=groups[name];g['files']+=1
        key=row['dev'],row['inode']
        if key not in g['inodes']:g['bytes']+=row['allocated_bytes'];g['inodes'].add(key)
    for g in groups.values():g.pop('inodes')
    result=dict(prepared=True,authorized=True,target_bytes=TARGET,delete_files=len(candidates),
        planned_reclaim_bytes=reclaimed,snapshot_workspace_file_bytes=total,
        expected_remaining_file_bytes=expected,expected_remaining_decimal_GB=expected/1e9,
        margin_for_directories_and_growth_bytes=TARGET-expected,groups=dict(groups),
        canonical_prefixes=prefixes,declared_input_paths=sorted(exact),
        notes=['Old completed raw outputs only; not a withdrawal of valid scientific conclusions.',
        'Current campaign outputs are growing and are all protected.',
        'All plots/CSV/JSON/cards/logs/programs and small regression data remain.',
        'Old six-segment worldtube ledger and final60M substep particles remain.'])
    save(REPORT/'PLAN.json',result)
    for name,items in [('DELETE_FILES.jsonl',candidates),('PRESERVED_LARGE_FILES.jsonl',retained)]:
        with (REPORT/name).open('w') as f:
            for row in items:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['groups','canonical_prefixes','declared_input_paths']},ensure_ascii=False),flush=True)
    for name,g in sorted(groups.items(),key=lambda x:x[1]['bytes'],reverse=True)[:20]:print('RETIRE',round(g['bytes']/1e9,3),'GB',g['files'],name,flush=True)


def preflight():
    p=read(REPORT/'PLAN.json');assert p['prepared'] and p['authorized']
    assert p['expected_remaining_file_bytes']<TARGET-5_000_000_000,'Adjust preservation before executing: insufficient growth margin'
    rows=[json.loads(x) for x in (REPORT/'DELETE_FILES.jsonl').read_text().splitlines()]
    assert len(rows)==p['delete_files']==len({x['path'] for x in rows})
    protected=set(p['declared_input_paths']);prefixes=p['canonical_prefixes'];kept=[]
    for row in rows:
        path=Path(row['path']);match(row,path.lstat())
        assert path.is_relative_to(W/'results') and path.resolve()==path
        assert str(path) not in protected and not any(under(str(path),q) for q in prefixes)
        assert path.suffix.lower() in RAW
        assert not any(x in ['inputs','source','src','scripts','bin','frozen_source'] or x.startswith('build') for x in path.relative_to(W).parts)
    for line in (REPORT/'PRESERVED_LARGE_FILES.jsonl').read_text().splitlines():
        row=json.loads(line);path=Path(row['path'])
        if not path.is_relative_to(CURRENT):match(row,path.lstat());kept.append(row)
    save(REPORT/'PRESERVED_STATIC_BEFORE.json',kept)
    before=read(REPORT/'CURRENT_FROZEN_BEFORE.json');now=current_frozen_checks(W,REPORT);assert now==before
    save(REPORT/'PREFLIGHT.json',dict(passed=True,delete_files=len(rows),preserved_large_static_files=len(kept),
        current_frozen_dependencies_match=True,plan_sha256=sha(REPORT/'PLAN.json'),filelist_sha256=sha(REPORT/'DELETE_FILES.jsonl'),
        identity_and_scope_verified=True,target_bytes=TARGET))
    print('PREFLIGHT_PASS',len(rows),len(kept),flush=True)


def execute():
    check=read(REPORT/'PREFLIGHT.json');p=read(REPORT/'PLAN.json');assert check['passed']
    assert sha(REPORT/'PLAN.json')==check['plan_sha256'] and sha(REPORT/'DELETE_FILES.jsonl')==check['filelist_sha256']
    assert read(REPORT/'AUTHORIZATION.json')['execute_old_output_cleanup']
    rows=[json.loads(x) for x in (REPORT/'DELETE_FILES.jsonl').read_text().splitlines()]
    removed=collections.Counter();freed=0;count=0
    with (REPORT/'DELETE_JOURNAL.jsonl').open('x') as f:
        for row in rows:
            f.write(json.dumps(dict(event='intent',**row),ensure_ascii=False)+'\n');f.flush();os.fsync(f.fileno())
            amount=unlink_verified(W,row,removed);count+=1;freed+=amount
            f.write(json.dumps(dict(event='deleted',reclaimed_bytes=amount,**row),ensure_ascii=False)+'\n');f.flush();os.fsync(f.fileno())
            if count%250==0:print('DELETED',count,'of',len(rows),'freed_GB',round(freed/1e9,3),flush=True)
    assert count==len(rows) and freed==p['planned_reclaim_bytes']
    assert all(not Path(x['path']).exists() for x in rows)
    kept=read(REPORT/'PRESERVED_STATIC_BEFORE.json')
    for row in kept:match(row,Path(row['path']).lstat())
    after=current_frozen_checks(W,REPORT);save(REPORT/'CURRENT_FROZEN_AFTER.json',after)
    assert after==read(REPORT/'CURRENT_FROZEN_BEFORE.json')
    result=dict(complete=True,deleted_files=count,reclaimed_bytes=freed,reclaimed_decimal_GB=freed/1e9,
        preserved_static_large_files=len(kept),current_frozen_dependencies_unchanged=True,
        original_inputs_and_key_science_preserved=True,plots_tables_cards_logs_and_code_preserved=True,
        no_running_output_read=True,no_directories_deleted=True,no_links_followed=True,
        actual_workspace_size_check_pending=True)
    save(REPORT/'EXECUTION_COMPLETE.json',result);print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['plan','preflight','execute']);a=p.parse_args()
    {'plan':plan,'preflight':preflight,'execute':execute}[a.mode]()
