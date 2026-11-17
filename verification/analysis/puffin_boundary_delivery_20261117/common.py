"""Shared I/O for the closed Puffin boundary delivery."""
import csv
import hashlib
import json
from pathlib import Path

W=Path('/public/home/wfollowm/APRL_Workspace')
SOURCE=W/'results/puffin_boundary_propagation_20261117_r1'
OLD=W/'results/puffin_boundary_explanation_20261117_r1'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for x in iter(lambda:f.read(1048576),b''):h.update(x)
    return h.hexdigest()

def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def table(path,rows):
    with Path(path).open('w',newline='') as f:
        z=csv.DictWriter(f,fieldnames=list(rows[0]));z.writeheader();z.writerows(rows)
