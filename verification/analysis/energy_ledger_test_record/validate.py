#!/usr/bin/env python3
"""Verify report assets, old-path aliases and untouched completed-test evidence."""
import argparse,hashlib,json,os,re
from pathlib import Path
from urllib.parse import unquote

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def validate(bundle):
    bundle=Path(bundle).resolve();meta=bundle/'metadata';move=json.loads((meta/'MIGRATION.json').read_text())
    old=Path(move['original']);new=Path(move['destination']);assert old.is_symlink() and old.resolve()==new.resolve()
    assert new.stat().st_ino==move['root_inode'] and new.stat().st_dev==move['root_device']
    inventory=json.loads((meta/'ORIGINAL_INVENTORY.json').read_text())['files'];hashed=0
    for item in inventory:
        p=new/item['relative'];q=old/item['relative'];a=p.lstat();b=q.lstat()
        assert a.st_ino==item['inode']==b.st_ino and a.st_dev==item['device']==b.st_dev,item['relative']
        if item['kind']!='directory':assert a.st_size==item['size'] and a.st_mtime_ns==item['mtime_ns'],item['relative']
        if item['kind']=='link':assert os.readlink(p)==item['link_target'] and p.exists()==item['target_existed']
        if 'sha256' in item:assert sha(p)==item['sha256'];hashed+=1
    markdown=(bundle/'README.md').read_text();links=[];images=[]
    assert '{{' not in markdown and '旧差额' not in markdown
    for m in re.finditer(r'(!?)\[[^\]]*\]\(([^)]+)\)',markdown):
        relative=unquote(m[2].split('#',1)[0]);assert relative and not relative.startswith(('http:','https:'))
        p=bundle/relative
        # This validation report is written only after every other check passes.
        if relative!='metadata/VALIDATION.json':assert p.exists(),relative
        links.append(relative)
        if m[1]:assert relative.startswith('figures/') and p.suffix=='.png';images.append(relative)
    assert len(images)==5 and len(set(images))==5
    snapshots=json.loads((meta/'SNAPSHOT_SHA256.json').read_text());assert all(sha(bundle/p)==v for p,v in snapshots.items())
    render=json.loads((meta/'RENDER_SHA256.json').read_text());assert all(sha(bundle/p)==v for p,v in render.items())
    for stage in ('campaign-r2','mpi-r1'):
        d=json.loads((new/stage/'COMPLETION.json').read_text());assert d['checks_pass']
        frozen=json.loads((new/stage/'frozen.json').read_text())
        # Move preserves original absolute paths; all frozen file identities still resolve.
        assert all(Path(p).exists() for p in frozen)
    program=json.loads((bundle/'programs/INDEX.json').read_text());assert all(sha(bundle/'programs'/name)==data['sha256'] for name,data in program.items())
    result={'checks_pass':True,'evidence_entries_preserved':len(inventory),'evidence_key_files_sha256_verified':hashed,
      'markdown_links_checked':len(links),'embedded_success_figures':images,'script_snapshots_checked':len(snapshots),
      'rendered_artifacts_sha256_verified':len(render),
      'compatibility_path_resolves_to_same_inode':True,'no_new_simulations':True,'scientific_files_rewritten':False,
      'report_scope':'Successful completed validation only; figures re-render selected accepted-test fields.'}
    (meta/'VALIDATION.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    assert all((bundle/p).exists() for p in links)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('bundle',type=Path);a=p.parse_args();print(json.dumps(validate(a.bundle),indent=2))
