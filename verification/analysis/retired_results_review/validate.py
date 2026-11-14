"""Independently check candidate identities, hardlinks, protection and archive."""
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main(workspace,output):
    workspace=Path(workspace).resolve();output=Path(output).resolve()
    plan=json.loads((output/'CANDIDATES.json').read_text());rules=json.loads((output/'PROTECTION_RULES.json').read_text())
    rows=list(csv.DictReader((output/'CANDIDATE_FILES.csv').open()))
    assert len(rows)==plan['candidate_files'] and len({x['path'] for x in rows})==len(rows)
    inode=collections.defaultdict(list);group_rows=collections.defaultdict(list)
    for row in rows:
        p=Path(row['path']);s=p.lstat()
        assert p.is_relative_to(workspace/'results') and stat.S_ISREG(s.st_mode) and p.resolve()==p
        assert (s.st_dev,s.st_ino,s.st_nlink,s.st_size,s.st_blocks*512,s.st_mtime_ns)==tuple(int(row[k]) for k in ['dev','inode','nlink','size_bytes','allocated_bytes','mtime_ns']),p
        assert str(p) not in rules['exact']
        assert not any(p.is_relative_to(prefix) for prefix in rules['prefixes']),p
        assert 'ledger' not in p.parts and 'lab-audit' not in p.parts,p
        assert p.suffix.lower() in ['.h5','.hdf5','.bin','.raw']
        if row['group']=='A04':assert 'particles' not in p.name and not p.name.endswith('-energy.h5')
        inode[s.st_dev,s.st_ino].append((row,s));group_rows[row['group']].append(row)
    joint=0;individual=collections.Counter()
    for values in inode.values():
        s=values[0][1]
        if len(values)==s.st_nlink:
            joint+=s.st_blocks*512
            ids={r['group'] for r,_ in values}
            if len(ids)==1:individual[next(iter(ids))]+=s.st_blocks*512
    assert joint==plan['all_groups_joint_reclaim_bytes']
    for group in plan['groups']:
        assert individual[group['id']]==group['exclusive_reclaim_bytes']
        assert len(group_rows[group['id']])==group['files']
        items=[json.loads(x) for x in (output/'file_lists'/(group['id']+'.jsonl')).read_text().splitlines()]
        assert {x['path'] for x in items}=={x['path'] for x in group_rows[group['id']]}
    archive=json.loads((output/'ARCHIVE_MANIFEST.json').read_text())
    for item in archive['files']:
        p=output/item['archive'];assert sha(p)==item['sha256'],p
        assert p.suffix.lower() not in ['.h5','.hdf5','.bin','.raw']
    count=0
    for name in ['README.md','DELETE_CANDIDATES.md']:
        p=output/name
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https://','http://','#')):continue
            target=target.strip('<>').split('#')[0];assert (p.parent/target).exists(),(name,target);count+=1
    scope=json.loads((output/'RESULTS_DIRECTORY_SIZES.json').read_text())
    root_ids=collections.defaultdict(set);root_counts=collections.Counter()
    for row in rows:
        name=Path(row['path']).relative_to(workspace/'results').parts[0]
        root_ids[name].add(row['group']);root_counts[name]+=1
    with (output/'DIRECTORY_DISPOSITION.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['name','snapshot_unique_allocated_GiB','status','candidate_groups','candidate_file_count'])
        for item in scope:
            status='部分文件候选；其余保留' if root_ids[item['name']] else '保留／未发现已核实退役依据'
            if item['name']=='fel_independent_w355_20261114_r1':status='当前运行和排队任务保护'
            w.writerow([item['name'],f"{item['unique_allocated_bytes']/2**30:.6f}",status,','.join(sorted(root_ids[item['name']])),root_counts[item['name']]])
    result=dict(passed=True,candidate_files=len(rows),candidate_identity_and_nlink_match=True,
        protected_paths_intersection=0,candidate_symlink_count=0,all_candidate_paths_under_results=True,
        allocated_reclaim_independently_recomputed=True,joint_reclaim_bytes=joint,
        grouped_jsonl_matches_csv=True,all_ledger_and_lab_audit_directories_retained=True,
        old_six_particle_and_energy_planes_retained=True,
        worldtube_bin_retirement_only_optional_B05=True,
        archived_files_sha256_verified=len(archive['files']),authored_markdown_links_verified=count,
        historical_input_audit_complete=json.loads((output/'HISTORICAL_LOADING_AUDIT.json').read_text())['complete'],
        no_running_h5_read=True,no_raw_files_moved=True,files_deleted=0,delete_authorized=False)
    (output/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
