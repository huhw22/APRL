#!/usr/bin/env python3
"""Copy an immutable historical directory and verify every regular file by SHA256."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
from datetime import datetime, timezone


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('destination', type=Path)
    p.add_argument('--manifest', type=Path, required=True)
    a = p.parse_args()
    src = a.source.resolve(strict=True)
    dst = a.destination.resolve()
    if dst.exists() or dst == src or src in dst.parents:
        raise ValueError('Destination must be a new directory outside the source')
    inventory = []
    for path in sorted(src.rglob('*')):
        stat = path.lstat()
        inventory.append((path, stat.st_size, stat.st_mtime_ns, stat.st_mode))
    dst.mkdir(parents=True)
    report = {'source': str(src), 'destination': str(dst),
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'complete': False, 'files': []}
    def save():
        tmp = a.manifest.with_suffix('.tmp')
        tmp.write_text(json.dumps(report, indent=2) + '\n')
        os.replace(tmp, a.manifest)
    save()
    for path, size, mtime, mode in inventory:
        relative = path.relative_to(src)
        out = dst / relative
        if path.is_symlink():
            out.symlink_to(os.readlink(path))
            report['files'].append({'path': str(relative), 'symlink': os.readlink(path)})
            save()
            continue
        if path.is_dir():
            out.mkdir(exist_ok=True)
            continue
        before = path.stat()
        if (before.st_size, before.st_mtime_ns, before.st_mode) != (size, mtime, mode):
            raise RuntimeError(f'Source changed before copy: {relative}')
        h = hashlib.sha256()
        with path.open('rb') as inp, out.open('xb') as copied:
            for block in iter(lambda: inp.read(8 * 1024 * 1024), b''):
                h.update(block)
                copied.write(block)
        shutil.copystat(path, out)
        after = path.stat()
        if (after.st_size, after.st_mtime_ns) != (size, mtime):
            raise RuntimeError(f'Source changed during copy: {relative}')
        actual = digest(out)
        if out.stat().st_size != size or actual != h.hexdigest():
            raise RuntimeError(f'Copy verification failed: {relative}')
        report['files'].append({'path': str(relative), 'size_bytes': size,
                                'mtime_ns': mtime, 'sha256': actual, 'verified': True})
        save()
        print(f'verified {relative} ({size} bytes)', flush=True)
    current_paths = {str(path.relative_to(src)) for path in src.rglob('*')}
    if current_paths != {str(path.relative_to(src)) for path, *_ in inventory}:
        raise RuntimeError('Source directory inventory changed during copying')
    for path, size, mtime, mode in inventory:
        stat = path.lstat()
        if (stat.st_size, stat.st_mtime_ns, stat.st_mode) != (size, mtime, mode):
            raise RuntimeError(f'Source changed during snapshot: {path}')
    report['complete'] = True
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    report['total_regular_bytes'] = sum(x.get('size_bytes', 0) for x in report['files'])
    save()
    print('Snapshot complete and verified', flush=True)


if __name__ == '__main__':
    main()
