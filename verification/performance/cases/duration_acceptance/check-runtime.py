#!/usr/bin/env python3
"""Reject loader errors before negative tests can mistake them for success."""
from pathlib import Path
import os
import socket
import subprocess
import sys

build, result = map(Path, sys.argv[1:])
targets = [build / name for name in (
    'aprl', 'particle_text_to_hdf5', 'uniform_lab_plane_bunch_generator',
    'sase_bunch_generator', 'energy_ledger_report')]
targets += list((build / 'verification').glob('aprl_*'))
targets += list((build / 'verification').glob('libaprl_*.so'))
targets += [build / 'postprocess' / name / name for name in (
    'field_reconstruction', 'field_plane_analysis', 'field_power_compare',
    'trajectory_radiation', 'energy_closure')]
failed = []
with (result / ('runtime-links-' + socket.gethostname() + '.txt')).open('w') as log:
    for target in targets:
        if not target.is_file() or not os.access(target, os.R_OK):
            raise SystemExit('Missing runtime target: ' + str(target))
        checked = subprocess.run(['ldd', '-r', str(target)], text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.write(str(target) + '\n' + checked.stdout + '\n')
        if checked.returncode or 'not found' in checked.stdout or 'undefined symbol' in checked.stdout:
            failed.append(str(target))
if failed:
    raise SystemExit('Runtime dependency preflight failed: ' + ', '.join(failed))
print('Runtime dependencies passed on ' + socket.gethostname())
