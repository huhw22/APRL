import hashlib
import json
from pathlib import Path

import numpy as np

W = Path('/public/home/wfollowm/APRL_Workspace')


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x))+'\n')


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for a in iter(lambda: f.read(4*1024*1024), b''):
            h.update(a)
    return h.hexdigest()


def require(value, message):
    if not value:
        raise RuntimeError(message)


def frozen(path):
    for name, value in load(path).items():
        require(sha(name) == value, 'Frozen identity changed: '+name)


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k:native(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [native(v) for v in value]
    return value
