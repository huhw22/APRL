"""Bounded GENESIS controls and read-only audits, all on the allocation."""
import concurrent.futures
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

from prepare import save, sha

R = Path(sys.argv[1]).resolve()
P = json.loads((R/'PLAN.json').read_text())
assert os.environ.get('SLURM_JOB_ID')
for path, digest in json.loads((R/'FROZEN.json').read_text()).items():
    assert sha(path) == digest, path


def launch(case):
    d = Path(case['directory'])
    srun = ['srun','--mpi=pmix_v3','--exact','--exclusive','--nodes=1',
           '--ntasks=32','--ntasks-per-node=32','--cpus-per-task=1',
           '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
    with (d/'binding.log').open('x') as f:
        subprocess.run(srun+['/public/home/wfollowm/APRL_Workspace/build/cluster-tools/binding_probe','32','1'],
                       cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
    cmd = ['timeout','--signal=TERM','--kill-after=30',str(case['timeout_seconds']),
           *srun,case['binary'],'reference.in']
    record = dict(command=cmd, start_unix=time.time(), job_id=os.environ['SLURM_JOB_ID'])
    save(d/'STARTED.json', record)
    start = time.monotonic()
    with (d/'genesis.log').open('x') as f:
        process = subprocess.run(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT)
    record.update(exit_code=process.returncode, seconds=time.monotonic()-start)
    save(d/'genesis.log.launch.json', record)
    assert process.returncode == 0, case['name']
    text = (d/'genesis.log').read_text()
    assert 'End of Track' in text and 'nan' not in text.lower(), case['name']
    with h5py.File(P['genesis_input'],'r') as a, h5py.File(d/'imported.par.h5','r') as b:
        # The original active input groups are imported with their exact values.
        for index in range(P['active_slices']):
            key = f'slice{512+index+1:06d}'
            for field in ['current','x','y','px','py','gamma','theta']:
                assert np.array_equal(a[key][field][:], b[key][field][:]), (case['name'],key,field)
    with h5py.File(d/'reference.out.h5','r') as f:
        for key in ['Field/power','Beam/bunching','Beam/bunchingphase','Beam/energy','Beam/energyspread']:
            assert np.isfinite(f[key][:]).all(), (case['name'],key)
        z = f['Lattice/zplot'][:]
        assert z[-1] >= P['original_plan']['geometry']['field_planes'][5][1]-.01
    save(d/'COMPLETE.json', dict(complete=True, input_import_bitwise=True,
                               launch=record, physical_acceptance=False))
    print(json.dumps(dict(case=case['name'],complete=True,seconds=record['seconds'])),flush=True)
    return case['name']


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    futures = [pool.submit(launch, case) for case in P['cases']]
    # CPU-light audit runs while the four separate Slurm steps propagate.
    from input_audit import audit
    audit(R, P)
    outcomes = [f.result() for f in futures]

save(R/'GENESIS_COMPLETE.json', dict(complete=True,cases=outcomes,
                                   physical_acceptance=False))
from compare import analyze
analyze(R, P)
print('ANALYSIS_COMPLETE',flush=True)
