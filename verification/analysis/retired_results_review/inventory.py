"""Read-only filesystem metadata inventory; never follows links or opens H5."""
import argparse
import json
import os
from pathlib import Path
import sqlite3
import stat
import time


def main(workspace,output):
    workspace=Path(workspace).resolve();output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    db=sqlite3.connect(output/'inventory.sqlite')
    db.execute('CREATE TABLE files(path TEXT PRIMARY KEY,kind TEXT,dev INTEGER,ino INTEGER,nlink INTEGER,size INTEGER,allocated INTEGER,mtime_ns INTEGER,target TEXT)')
    db.execute('CREATE TABLE errors(path TEXT,message TEXT)')
    device=workspace.stat().st_dev;queue=[workspace];batch=[];n=0;nd=0;start=time.monotonic()
    while queue:
        directory=queue.pop();nd+=1
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    path=Path(entry.path)
                    if path==output or entry.name in ['.git','.codex','.agents','__pycache__']:
                        continue
                    try:
                        s=entry.stat(follow_symlinks=False)
                        if s.st_dev!=device:
                            db.execute('INSERT INTO errors VALUES (?,?)',(str(path),'Different filesystem; not traversed'));continue
                        if stat.S_ISDIR(s.st_mode):queue.append(path);continue
                        kind='file' if stat.S_ISREG(s.st_mode) else ('symlink' if stat.S_ISLNK(s.st_mode) else 'other')
                        target=os.readlink(path) if kind=='symlink' else ''
                        batch.append((str(path),kind,s.st_dev,s.st_ino,s.st_nlink,s.st_size,s.st_blocks*512,s.st_mtime_ns,target));n+=1
                        if len(batch)>=5000:
                            db.executemany('INSERT INTO files VALUES (?,?,?,?,?,?,?,?,?)',batch);db.commit();batch=[]
                            if n%50000==0:print('INVENTORY',n,'files',nd,'directories',round(time.monotonic()-start,1),'seconds',flush=True)
                    except OSError as exc:db.execute('INSERT INTO errors VALUES (?,?)',(str(path),str(exc)))
        except OSError as exc:db.execute('INSERT INTO errors VALUES (?,?)',(str(directory),str(exc)))
    if batch:db.executemany('INSERT INTO files VALUES (?,?,?,?,?,?,?,?,?)',batch)
    db.execute('CREATE INDEX inode_index ON files(dev,ino)');db.commit()
    result=dict(complete=True,workspace=str(workspace),files=n,directories=nd,seconds=time.monotonic()-start,
        metadata_only=True,h5_contents_read=False,symlinks_followed=False,
        errors=db.execute('SELECT count(*) FROM errors').fetchone()[0],
        unique_file_allocated_bytes=db.execute("SELECT sum(b) FROM (SELECT max(allocated) b FROM files WHERE kind='file' GROUP BY dev,ino)").fetchone()[0],
        note='Metadata snapshot; running outputs may grow. Reclaim calculations must account for hardlinks and protected dependencies. Directories and excluded git/cache metadata are not counted.')
    (output/'INVENTORY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    roots={}
    for path,dev,ino,size,allocated,nlink in db.execute("SELECT path,dev,ino,size,allocated,nlink FROM files WHERE kind='file'"):
        rel=Path(path).relative_to(workspace)
        if len(rel.parts)<2 or rel.parts[0]!='results':continue
        key=rel.parts[1];item=roots.setdefault(key,dict(files=0,size=0,inodes={}))
        item['files']+=1;item['size']+=size;item['inodes'][dev,ino]=allocated
    summary=[dict(name=name,files=x['files'],apparent_bytes=x['size'],unique_allocated_bytes=sum(x['inodes'].values())) for name,x in roots.items()]
    summary.sort(key=lambda x:x['unique_allocated_bytes'],reverse=True)
    (output/'RESULTS_DIRECTORY_SIZES.json').write_text(json.dumps(summary,indent=2)+'\n')
    for row in summary[:30]:print('LARGEST',row['name'],round(row['unique_allocated_bytes']/2**30,3),'GiB',row['files'],'files',flush=True)
    db.close()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('workspace');p.add_argument('output');a=p.parse_args();main(a.workspace,a.output)
