#!/usr/bin/env python3
"""Run and verify a frozen free-drift observer replay, within Slurm only."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np
sys.path.insert(0, str(Path(__file__).parent/'transport'))
from prepare import require, save, sha
from run import launch, binding
import review


def main():
    root = Path(sys.argv[1]).resolve()
    plan = json.loads((root/'plan.json').read_text())
    frozen = json.loads((root/'frozen.json').read_text())
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NUM_NODES') == '1', 'Need one Slurm node')
    require(os.environ.get('SLURM_NTASKS') == '16' and os.environ.get('SLURM_CPUS_PER_TASK') == '8', 'Need 16 x 8')
    require(all(sha(k) == v for k, v in frozen.items()), 'Frozen inputs changed')
    exe = root/'execution'; exe.mkdir()
    env = os.environ.copy()
    env.pop('APRL_GRID_WORK_AUDIT', None)
    env.pop('APRL_FIELD_STAGE_AUDIT', None)
    env['LD_LIBRARY_PATH'] = plan['build']+'/runtime-libs:'+env.get('LD_LIBRARY_PATH', '')
    def srun(ranks=16, threads=8):
        return ['srun', '--mpi='+env.get('SRUN_MPI', 'pmix_v3'), '--exact', '--nodes=1',
                '--ntasks='+str(ranks), '--cpus-per-task='+str(threads), '--cpu-bind=cores',
                '--kill-on-bad-exit=1', '--wait=30']
    try:
        launch(srun()+[plan['binding_probe'], '16', '8'], exe, env, 45, 'binding.log')
        save(exe/'BINDING.json', binding(exe/'binding.log', 1))
        single = env.copy(); single['OMP_NUM_THREADS'] = '1'
        launch(['ctest', '--test-dir', plan['build'], '-L', 'required', '--output-on-failure', '-j', '1'], exe, single, 300, 'required.log')
        analytic = exe/'observer-analytic'; analytic.mkdir()
        launch(srun(1, 1)+[plan['observer_checks'],str(analytic)], exe, single, 45, 'observer-checks.log')
        launch(srun(1, 1)+[plan['probe']], exe, single, plan['limits_s']['probe'], 'probe.jsonl')
        probe = [json.loads(x) for x in (exe/'probe.jsonl').read_text().splitlines() if x.startswith('{')]
        rows = [x for x in probe if x['type'] == 'result']
        require(len(rows) == 20 and all(x['pass'] for x in rows), 'Missing/failed independent probe checks')
        save(exe/'PROBE_REVIEW.json', {'checks_pass': True, 'results': rows})
        case = Path(plan['case'])
        env['APRL_GRID_WORK_AUDIT'] = str(case/'mesh-work')
        if plan.get('field_stages'): env['APRL_FIELD_STAGE_AUDIT']=str(case/'field-stages')
        launch(srun()+[plan['binary'], str(case/'run.yaml')], case, env, plan['limits_s']['case'])
        item = {'steps': 257, 'ledger': True, 'energy_planes': False, 'binary': 'new', 'magnet': 'free-drift'}
        checked = review.inspect(case, item, plan, export_csv=False)
        save(case/'CASE_REVIEW.json', checked)
        comparison = review.compare(Path(plan['baseline']), case)
        require(comparison['checks_pass'], 'Observer changed raw physical outputs')
        save(root/'PHYSICS_BITWISE_REVIEW.json', comparison)
        # Timers are intentionally affected; every scientific ledger component must agree.
        with h5py.File(Path(plan['baseline'])/'ledger/energy-ledger.h5') as f, h5py.File(case/'ledger/energy-ledger.h5') as g:
            a, b = f['energy_ledger/records'][:], g['energy_ledger/records'][:]
        require(a.shape == b.shape and a.dtype == b.dtype, 'Ledger schema changed')
        ignored = [k for k in a.dtype.names if 'seconds' in k]
        agreement = {k: np.ascontiguousarray(a[k]).tobytes() == np.ascontiguousarray(b[k]).tobytes()
                     for k in a.dtype.names if k not in ignored}
        require(all(agreement.values()), 'Observer changed energy ledger: '+repr([k for k,v in agreement.items() if not v]))
        save(root/'LEDGER_BITWISE_REVIEW.json', {'checks_pass': True, 'fields': agreement, 'ignored_timers': ignored})
        require(all(sha(k) == v for k, v in frozen.items()), 'Frozen inputs changed during run')
        save(root/'COMPLETION.json', {'checks_pass': True, 'job_id': os.environ['SLURM_JOB_ID'],
              'binary_sha256': plan['binary_sha256'], 'scope': 'Observer validation; not energy conservation acceptance'})
    except BaseException as error:
        save(root/'FAILED_OR_INCOMPLETE.json', {'checks_pass': False, 'error': str(error)})
        raise


if __name__ == '__main__':
    main()
