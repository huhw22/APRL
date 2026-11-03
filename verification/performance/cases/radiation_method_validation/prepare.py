"""Shared small utilities; importing this module never prepares or runs PIC."""
import hashlib
import json
from pathlib import Path
import numpy as np

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value

def save(path, value):
    Path(path).write_text(json.dumps(native(value), indent=2, allow_nan=False)+'\n')
