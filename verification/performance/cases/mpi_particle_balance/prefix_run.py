"""Execute an explicit duration-bounded prefix, never a completed FEL claim."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

from common import require, save, sha
from binding import binding
from prefix_closed import duration_status, verify_prefix_case
from analyze_closed import compare_cases


def launch(command, directory, env, limit, filename='aprl.log'):
    full = ['timeout', '--signal=TERM', '--kill-after=10', str(limit)] + command
    record = dict(command=full, cwd=str(directory), limit_s=limit, started_unix=time.time(),
        environment={name:env.get(name) for name in ['OMP_NUM_THREADS', 'OMP_DYNAMIC',
            'OMP_PLACES', 'OMP_PROC_BIND', 'OMPI_MCA_btl', 'APRL_PROFILE_INTERVAL', 'LD_LIBRARY_PATH']})
    process = None
    started = time.monotonic()
    try:
        with (directory / filename).open('x') as stream:
            process = subprocess.Popen(full, cwd=directory, env=env, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            record['exit_code'] = process.wait(timeout=limit + 15)
    finally:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        record['seconds'] = time.monotonic() - started
        save(directory / (filename + '.launch.json'), record)
    return record


def main(root):
    root = Path(root).resolve()
    plan = json.loads((root / 'plan.json').read_text())
    require(plan.get('accepted_for_execution') is True and plan.get('bounded_duration_prefix') is True,
            'Prefix campaign is not explicitly frozen and accepted')
    require(os.environ.get('SLURM_JOB_ID'), 'Use a bounded allocation')
    frozen = json.loads((root / 'FROZEN_INPUTS.json').read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()), 'Frozen prefix inputs changed')
    build = json.loads((Path(plan['build_record']) / 'COMPLETE.json').read_text())
    require(build.get('passed') is True and sha(plan['binary']) == build['binary_sha256'], 'Unaccepted prefix program')
    env = {key:value for key, value in os.environ.items() if not key.startswith('APRL_')}
    env['APRL_PROFILE_INTERVAL'] = str(plan.get('profile_interval', 0))
    require(env.get('OMP_NUM_THREADS') == '1', 'Prefix uses single-thread MPI only')
    results = {}
    for case in plan['cases']:
        require(case.get('bounded_duration_prefix') is True, 'A physical-stop case entered the prefix runner')
        directory = Path(case['directory'])
        nodes, ranks = case['nodes'], case['ranks']
        require(case['threads'] == 1 and ranks % nodes == 0, 'Invalid prefix MPI layout')
        command = ['srun', '--mpi=pmix_v3', '--exclusive', '--exact', '--nodes=' + str(nodes),
            '--ntasks=' + str(ranks), '--ntasks-per-node=' + str(ranks // nodes),
            '--cpus-per-task=1', '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']
        print('PREFIX START', case['name'], flush=True)
        probe = launch(command + [str(Path(plan['workspace']) / 'build/cluster-tools/binding_probe'),
                                 str(ranks), '1'], directory, env, 120, 'binding.log')
        require(probe['exit_code'] == 0 and probe['seconds'] < probe['limit_s'], 'Binding probe failed')
        save(directory / 'BINDING.json', binding(directory / 'binding.log', nodes, ranks, 1))
        runtime = launch(command + [case['binary'], str(directory / 'run.yaml')],
                         directory, env, case['limit_seconds'])
        # Proof of all-rank close must pass before any scientific HDF5 reads.
        termination = duration_status(directory, case, runtime)
        result = verify_prefix_case(case, runtime)
        result['runtime_seconds'] = runtime['seconds']
        save(directory / 'RESULT.json', result)
        save(directory / 'PREFIX_COMPLETE.json', dict(prefix_complete=True,
            closed_complete_steps=case['expected_steps'], physical_complete=False,
            program_integrity_pass=True, physical_acceptance=False,
            job=os.environ['SLURM_JOB_ID'], runtime=runtime, duration_exit=termination,
            binary_sha256=build['binary_sha256'], bounded_duration_prefix=True))
        results[case['name']] = result
        print('CLOSED PREFIX', case['name'], flush=True)
    comparisons = compare_cases(plan)
    require(all(sha(path) == digest for path, digest in frozen.items()), 'Frozen prefix inputs changed during run')
    save(root / 'PREFIX_COMPLETE.json', dict(prefix_complete=True, program_integrity_pass=True,
        physical_complete=False, physical_acceptance=False, bounded_duration_prefix=True,
        job=os.environ['SLURM_JOB_ID'], binary_sha256=build['binary_sha256'],
        cases=results, comparisons=comparisons,
        verification_scope=plan['verification_scope']))
    print('ALL PREFIXES CLOSED; PHYSICAL STOP NOT REACHED', flush=True)


if __name__ == '__main__':
    root = Path(sys.argv[1]).resolve()
    try:
        main(root)
    except Exception as error:
        save(root / 'FAILED.json', dict(passed=False, physical_complete=False,
             error=str(error), traceback=traceback.format_exc(), job=os.environ.get('SLURM_JOB_ID')))
        raise
