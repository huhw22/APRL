"""Execute only the user's explicit approved groups from the reviewed CSV."""
import argparse
import collections
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
import sys

from build_plan import GROUPS,collect_protection,under

APPROVED=[f'A{i:02d}' for i in range(1,9)]+['B01','B03']


def save(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def read(path):return json.loads(Path(path).read_text())
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def fingerprint(s):
    return dict(dev=s.st_dev,inode=s.st_ino,nlink=s.st_nlink,size_bytes=s.st_size,
        allocated_bytes=s.st_blocks*512,mtime_ns=s.st_mtime_ns)


def match(row,s,removed=0):
    assert stat.S_ISREG(s.st_mode),'Not a regular file'
    actual=fingerprint(s);expected={k:int(row[k]) for k in actual}
    expected['nlink']-=removed
    assert actual==expected,('File identity changed',row['path'],expected,actual)


def current_frozen_checks(workspace,review):
    current=workspace/'results/fel_independent_w355_20261114_r1'
    manifests=[];paths={};cache={}
    for p in current.rglob('*FROZEN*.json'):
        if 'scripts' in p.relative_to(current).parts:continue
        value=read(p)
        if not isinstance(value,dict):continue
        manifests.append(dict(path=str(p),sha256=sha(p)))
        for path,expected in value.items():
            if not Path(path).is_absolute():continue
            s=Path(path).stat();key=(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns)
            if key not in cache:cache[key]=sha(path)
            assert cache[key]==expected,('Current frozen dependency changed',path)
            paths[path]=dict(sha256=expected,**fingerprint(s))
    assert manifests and paths
    return dict(manifests=manifests,files=paths,unique_inodes_hashed=len(cache),
        scope='Immutable cards, scripts, inputs and programs only; no running scientific output')


def verify_archive(review):
    manifest=read(review/'ARCHIVE_MANIFEST.json')
    for item in manifest['files']:assert sha(review/item['archive'])==item['sha256']
    return len(manifest['files'])


def retained_snapshot(review,rows):
    paths={r['path'] for r in rows if r['group'] not in APPROVED}
    for row in csv.DictReader((review/'EXCLUDED_PROTECTED_FILES.csv').open()):paths.add(row['path'])
    return {p:fingerprint(Path(p).lstat()) for p in sorted(paths)}


def prepare(workspace,review,execution):
    approval=read(execution/'APPROVAL.json')
    assert approval['approved_groups']==APPROVED and approval['explicit_human_authorization']
    expected_hashes={}
    for line in (review/'SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);expected_hashes[name]=digest
    for name in ['CANDIDATES.json','CANDIDATE_FILES.csv','ARCHIVE_MANIFEST.json','DELIVERY.json']:
        assert sha(review/name)==expected_hashes[name],('Reviewed manifest changed',name)
    rows=list(csv.DictReader((review/'CANDIDATE_FILES.csv').open()))
    chosen=[r for r in rows if r['group'] in APPROVED]
    assert len(chosen)==len({r['path'] for r in chosen})
    plan=read(review/'CANDIDATES.json');roots={g[0]:[str(workspace/'results'/p) for p in g[3]] for g in GROUPS}
    with sqlite3.connect(review/'inventory.sqlite') as db:prefixes,exact,manifests=collect_protection(workspace,db)
    protected_inodes=set()
    for p in exact:
        try:s=Path(p).stat();protected_inodes.add((s.st_dev,s.st_ino))
        except FileNotFoundError:continue
    by_inode=collections.defaultdict(list)
    for row in chosen:
        p=Path(row['path']);s=p.lstat();match(row,s)
        assert p.is_relative_to(workspace/'results') and p.resolve()==p
        assert any(under(str(p),root) for root in roots[row['group']])
        assert str(p) not in exact and not any(under(str(p),q) for q in prefixes)
        assert (s.st_dev,s.st_ino) not in protected_inodes
        assert p.suffix.lower() in ['.h5','.hdf5','.bin','.raw']
        assert not any(q in ['inputs','source','src','scripts','build','ledger','lab-audit'] or q.startswith('build-') for q in p.parts)
        assert p.name!='particles.h5' and not p.name.endswith('.par.h5') and not p.name.endswith('_electrons_0.h5')
        assert not p.name.startswith('input')
        if row['group']=='A04':assert 'particles' not in p.name and not p.name.endswith('-energy.h5')
        by_inode[s.st_dev,s.st_ino].append(row)
    allocated=sum(int(v[0]['allocated_bytes']) for v in by_inode.values() if len(v)==int(v[0]['nlink']))
    expected=sum(g['exclusive_reclaim_bytes'] for g in plan['groups'] if g['id'] in APPROVED)
    assert allocated==expected==int(approval['approved_reclaim_bytes'])
    archive=verify_archive(review);frozen=current_frozen_checks(workspace,review)
    keep=retained_snapshot(review,rows)
    save(execution/'PRESERVED_BEFORE.json',keep);save(execution/'CURRENT_FROZEN_BEFORE.json',frozen)
    save(execution/'FRESH_PROTECTION_RULES.json',dict(prefixes=prefixes,exact=exact,current_manifests=manifests))
    with (execution/'APPROVED_FILES.jsonl').open('x') as f:
        for row in chosen:f.write(json.dumps(row,ensure_ascii=False)+'\n')
    result=dict(passed=True,approved_groups=APPROVED,files=len(chosen),allocated_reclaim_bytes=allocated,
        protected_intersection=0,all_selected_identities_match=True,archive_sha_verified=archive,
        current_frozen_unique_inodes_hashed=frozen['unique_inodes_hashed'],
        preserved_files_snapshot=len(keep),files_deleted=0)
    save(execution/'PREFLIGHT.json',result);return chosen,result


def unlink_verified(workspace,row,removed):
    path=Path(row['path']);parts=path.relative_to(workspace).parts
    flags=os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW
    fd=os.open(workspace,flags)
    try:
        for part in parts[:-1]:
            new=os.open(part,flags,dir_fd=fd);os.close(fd);fd=new
        s=os.stat(parts[-1],dir_fd=fd,follow_symlinks=False)
        key=s.st_dev,s.st_ino;match(row,s,removed[key])
        os.unlink(parts[-1],dir_fd=fd);os.fsync(fd)
        removed[key]+=1
        return s.st_blocks*512 if s.st_nlink==1 else 0
    finally:os.close(fd)


def execute(workspace,review,execution):
    assert read(execution/'PREFLIGHT.json')['passed']
    assert not (execution/'DELETE_JOURNAL.jsonl').exists(),'Fresh execution only; inspect any existing partial journal before resuming'
    rows=[json.loads(x) for x in (execution/'APPROVED_FILES.jsonl').read_text().splitlines()]
    before=os.statvfs(workspace);removed=collections.Counter();done=[];freed=0
    with (execution/'DELETE_JOURNAL.jsonl').open('x') as journal:
        for row in rows:
            journal.write(json.dumps(dict(event='intent',**row),ensure_ascii=False)+'\n');journal.flush();os.fsync(journal.fileno())
            amount=unlink_verified(workspace,row,removed)
            done.append(row);freed+=amount
            journal.write(json.dumps(dict(event='deleted',reclaimed_allocated_bytes=amount,**row),ensure_ascii=False)+'\n');journal.flush();os.fsync(journal.fileno())
            if len(done)%100==0:print('DELETED',len(done),'of',len(rows),'allocated_GiB',round(freed/2**30,3),flush=True)
    assert all(not Path(r['path']).exists() for r in rows)
    keep=read(execution/'PRESERVED_BEFORE.json')
    for p,value in keep.items():assert fingerprint(Path(p).lstat())==value,('Unselected/protected file changed',p)
    after_frozen=current_frozen_checks(workspace,review);save(execution/'CURRENT_FROZEN_AFTER.json',after_frozen)
    assert after_frozen==read(execution/'CURRENT_FROZEN_BEFORE.json')
    assert verify_archive(review)==read(execution/'PREFLIGHT.json')['archive_sha_verified']
    after=os.statvfs(workspace);by_group=collections.Counter(r['group'] for r in done)
    result=dict(complete=True,approved_groups=APPROVED,deleted_files=len(done),deleted_by_group=dict(by_group),
        reclaimed_allocated_bytes=freed,reclaimed_allocated_GiB=freed/2**30,
        filesystem_available_before_bytes=before.f_bavail*before.f_frsize,
        filesystem_available_after_bytes=after.f_bavail*after.f_frsize,
        filesystem_net_available_change_bytes=after.f_bavail*after.f_frsize-before.f_bavail*before.f_frsize,
        net_change_scope='Concurrent jobs continue writing; allocated-file reclaim is the controlled count',
        preserved_files_verified=len(keep),current_frozen_dependencies_unchanged=True,
        unapproved_groups_retained=True,B04_retained=True,B05_retained=True,
        archived_evidence_unchanged=True,original_inputs_retained=True,source_programs_retained=True,
        no_whole_directory_removal=True,no_symlinks_followed=True,no_running_output_read=True)
    assert freed==read(execution/'PREFLIGHT.json')['allocated_reclaim_bytes']
    save(execution/'RESULT.json',result);print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['preflight','execute'])
    p.add_argument('workspace',type=Path);p.add_argument('review',type=Path);p.add_argument('execution',type=Path)
    a=p.parse_args();workspace=a.workspace.resolve();review=a.review.resolve();execution=a.execution.resolve()
    assert workspace==Path('/public/home/wfollowm/APRL_Workspace')
    assert review==workspace/'results/retired_results_review_20261114_r1'
    assert execution.parent==review
    if a.mode=='preflight':prepare(workspace,review,execution)
    else:execute(workspace,review,execution)
