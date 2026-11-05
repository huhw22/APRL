"""Execute default/substep invariance and notch integration on the original small card."""
import json
import os
from pathlib import Path
import re
import sys

import h5py
import numpy as np

from prepare import require, save, sha
from launcher import launch
from binding import binding
from analyze_closed import verify, bitwise
import report_discrete


def main(root):
    root = Path(root).resolve()
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Use a bounded Slurm allocation')
    plan = json.loads((root/'plan.json').read_text())
    frozen = json.loads((root/'FROZEN_INPUTS.json').read_text())
    require(all(sha(path)==digest for path,digest in frozen.items()), 'Smoke inputs changed')
    build_file = Path(plan['build_record'])/'COMPLETE.json'
    identity = json.loads(build_file.read_text())
    require(identity.get('passed') is True, 'Build/modules not accepted')
    require(str(Path(identity['binary']).resolve()) == plan['binary'], 'Wrong candidate path')
    require(sha(plan['binary']) == identity['binary_sha256'], 'Wrong candidate identity')
    frozen[plan['binary']] = identity['binary_sha256']
    frozen[str(build_file)] = sha(build_file)
    save(root/'FROZEN.json', frozen)
    results = {}
    for case in plan['cases']:
        directory = Path(case['directory'])
        cmd = ['srun', '--mpi=pmix_v3', '--exclusive', '--exact', '--nodes=1',
               '--ntasks=16', '--ntasks-per-node=16', '--cpus-per-task=8',
               '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']
        env = os.environ.copy()
        for key in list(env):
            if key.startswith('APRL_') and ('PROBE' in key or 'AUDIT' in key or 'DIAGNOSTIC' in key):
                env.pop(key, None)
        print('START', case['name'], flush=True)
        launch(cmd+[str(Path(plan['workspace'])/'build/cluster-tools/binding_probe'), '16', '8'], directory, env, 120, 'binding.log')
        save(directory/'BINDING.json', binding(directory/'binding.log', 1, 16, 8))
        runtime = launch(cmd+[case['binary'], str(directory/'run.yaml')], directory, env, case['limit_seconds'])
        text = (directory/'aprl.log').read_text()
        require(all(s in text for s in ['Configured stop reached:', 'phase=solve-complete', 'phase=mpi-finalize-complete']), 'Incomplete normal stop')
        require(not re.search('Simulation failure|OpenFabrics device|Segmentation fault|DUE TO TIME LIMIT', text), 'Smoke simulation error')
        closed = verify(directory)
        rows,sums,ranks,steps = report_discrete.load(directory/'ledger', directory/'ledger/energy-ledger.h5')
        require(steps==plan['expected_steps'] and ranks==16, 'Changed smoke geometry/duration')
        _,_,ledger = report_discrete.calculate(rows,sums,steps)
        substep_files = []
        if case['substeps']:
            for path in sorted((directory/'detectors').glob('*-substeps.h5')):
                with h5py.File(path, 'r') as f:
                    g=f['particle_substep_plane'];require(int(g['complete'][()])==1,'Unclosed substeps')
                    a=g['records'][:]
                    require(len(a)==int(g['committed_records'][()])==plan['expected_particles'],'Missing substep crossing')
                    require(len(a)==len(np.unique(a['source_id'])),'Duplicate substep source')
                substep_files.append(path.name)
            require(len(substep_files)==3,'Expected all three original particle planes')
        result=dict(program_integrity_pass=True,physical_acceptance=False,steps=steps,
            runtime_seconds=runtime['seconds'],closed_scientific_files=closed,
            substep_files=substep_files,strict_ledger=ledger)
        save(directory/'RESULT.json',result);results[case['name']]=result
        print('COMPLETE',case['name'],flush=True)
    default = bitwise(root/'old',root/'default')
    diagnostics = bitwise(root/'default',root/'substeps')
    require(all(sha(path)==digest for path,digest in frozen.items()),'Smoke inputs changed during execution')
    save(root/'COMPLETE.json',dict(program_integrity_pass=True,physical_acceptance=False,
        cases=results,default_science_bitwise=default,substep_switch_legacy_science_bitwise=diagnostics))
    print('ALL COMPLETE',flush=True)


if __name__=='__main__':
    main(sys.argv[1])
