#!/usr/bin/env python3
"""Read-only dynamic-loader check on each allocated host before simulation."""
import json
from pathlib import Path
import socket
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
plan = json.loads((root / 'plan.json').read_text())
host = socket.gethostname().split('.')[0]
allocated = (root / 'allocation-nodes.txt').read_text().splitlines()
if host not in allocated or len(set(allocated)) != 2:
    raise SystemExit('Runtime check did not run on the recorded two-node allocation')
targets = [Path(item['binary']) for item in plan['variants'].values()]
targets += [Path(plan['workspace']) / 'build/cluster-tools/binding_probe']
failed = []
with (root / ('runtime-links-' + host + '.txt')).open('x') as output:
    for target in targets:
        result = subprocess.run(['ldd', '-r', str(target)], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=20)
        output.write(str(target) + '\n' + result.stdout + '\n')
        if result.returncode or 'not found' in result.stdout or 'undefined symbol' in result.stdout:
            failed.append(str(target))
if failed:
    raise SystemExit('Runtime dependency preflight failed: ' + ', '.join(failed))
print('Runtime dependencies passed on ' + host)
