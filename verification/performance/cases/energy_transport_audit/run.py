#!/usr/bin/env python3
"""Run one frozen energy-transport stage inside Slurm; never submit jobs."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import json
from pathlib import Path
import re
import signal
import subprocess
import time
from prepare import require, save, sha
import review


def launch(command, directory, env, limit, filename='aprl.log'):
    full = ['timeout', '--signal=TERM', '--kill-after=10', str(limit)] + command
    record = {'command': full, 'cwd': str(directory), 'timeout_s': limit, 'started_unix': time.time(),
              'environment': {name: env.get(name) for name in ('OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_PLACES',
                  'OMP_PROC_BIND', 'OMPI_MCA_btl', 'APRL_PROFILE_INTERVAL', 'LD_LIBRARY_PATH')}}
    start, process = time.monotonic(), None
    try:
        with (directory/filename).open('x') as stream:
            process = subprocess.Popen(full, cwd=directory, env=env, stdin=subprocess.DEVNULL,
                stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            record['exit_code'] = process.wait(timeout=limit + 15)
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
        record['seconds'] = time.monotonic() - start
        save(directory/(filename + '.launch.json'), record)
    require(record['exit_code'] == 0 and record['seconds'] < limit, 'Launch failed/interrupted: ' + str(directory))
    return record


def binding(path, nodes):
    text = path.read_text()
    rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)', text)
    require('[binding] result=PASS failures=0' in text and len(rows) == 128, 'Incomplete binding check')
    require(all(row[-1] == 'PASS' for row in rows), 'Binding failures')
    require({(int(r[0]), int(r[1])) for r in rows} == {(r,t) for r in range(16) for t in range(8)}, 'Missing rank/thread')
    require(len({r[2] for r in rows}) == nodes, 'Wrong node count')
    require(len({(r[2],r[4],r[5]) for r in rows}) == 128, 'Physical core overlap')
    return {'pass': True, 'nodes': sorted({r[2] for r in rows}), 'ranks': 16, 'threads': 8, 'physical_cores': 128}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', choices=('science', 'performance', 'mpi'), required=True)
    args = parser.parse_args(); root = args.root.resolve(); stage = args.stage
    plan = json.loads((root/'plan.json').read_text()); nodes = 2 if stage == 'mpi' else 1
    require(plan.get('frozen') is True, 'Plan must be finalized against the actual candidate binary before execution')
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == nodes,
            'Require real Slurm allocation matching stage')
    require(int(os.environ.get('SLURM_NTASKS', '0')) == 16 and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 8,
            'Require fixed 16 MPI x 8 thread allocation')
    frozen = json.loads((root/'frozen.json').read_text())
    require(all(sha(path) == digest for path,digest in frozen.items()), 'Frozen campaign changed')
    if stage != 'science':
        gate = json.loads((root/'SCIENCE_REVIEW.json').read_text())
        require(gate['checks_pass'] and gate['binary_sha256'] == plan['binary_sha256'], 'Same-binary science gate required')
    execution = root/stage/'execution'; execution.mkdir(exist_ok=True)
    require(not (execution/'STARTED.json').exists(), 'Stage already started; use a new campaign')
    budget = plan['stage_budget_s'][stage]; start = time.monotonic(); deadline = start + budget
    save(execution/'STARTED.json', {'job_id': os.environ['SLURM_JOB_ID'], 'stage': stage,
        'started_unix': time.time(), 'active_budget_s': budget, 'binary_sha256': plan['binary_sha256']})
    def expired(signum, frame):
        raise TimeoutError('Stage active budget exhausted')
    def terminated(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGALRM, expired); signal.signal(signal.SIGTERM, terminated); signal.signal(signal.SIGINT, terminated)
    signal.setitimer(signal.ITIMER_REAL, budget)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='8', OMP_DYNAMIC='FALSE', OMP_PLACES='cores', OMP_PROC_BIND='close',
        APRL_PROFILE_INTERVAL='0', HDF5_USE_FILE_LOCKING='FALSE', **plan['mpi_environment'])
    def launcher(directory):
        return ['srun', '--mpi='+env.get('SRUN_MPI','pmix_v3'), '--exact', '--nodes='+str(nodes),
            '--ntasks=16', '--ntasks-per-node='+str(16//nodes), '--cpus-per-task=8', '--cpu-bind=cores',
            '--kill-on-bad-exit=1', '--wait=30', '--chdir='+str(directory)]
    completed = []
    try:
        runtime = {}
        for name, path in plan['binaries'].items():
            checked_env = env.copy(); checked_env['LD_LIBRARY_PATH'] = str(Path(path).parent/'runtime-libs') + ':' + env.get('LD_LIBRARY_PATH','')
            checked = subprocess.run(['ldd','-r',path], env=checked_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
            runtime[name] = {'returncode': checked.returncode, 'output': checked.stdout}
            require(checked.returncode == 0 and 'not found' not in checked.stdout and 'undefined symbol' not in checked.stdout,
                'Missing runtime dependencies: '+name)
        save(execution/'RUNTIME_LINKS.json', runtime)
        launch(launcher(execution)+[plan['binding_probe'],'16','8'], execution, env, 45, 'binding.log')
        save(execution/'BINDING.json', binding(execution/'binding.log', nodes))
        for item in plan['cases'][stage]:
            require(deadline-time.monotonic() > item['timeout_s']+25, 'Insufficient budget to start next bounded case')
            directory = Path(item['directory']); binary = Path(plan['binaries'][item['binary']])
            print('CASE_START '+item['name']+' directory='+str(directory), flush=True)
            case_env = env.copy(); case_env['LD_LIBRARY_PATH'] = str(binary.parent/'runtime-libs')+':'+env.get('LD_LIBRARY_PATH','')
            launch(launcher(directory)+[str(binary),str(directory/'run.yaml')], directory, case_env, item['timeout_s'])
            report = review.inspect(directory, item, plan)
            save(directory/'CASE_REVIEW.json', report)
            completed.append(item['name'])
            print('CASE_COMPLETE '+item['name']+' preparation_s='+str(report['reference_preparation_seconds'])+
                  ' step_s='+str(report['steps_40_240_seconds_per_step']), flush=True)
        result = review.stage_review(root, stage, plan)
        require(result['checks_pass'], 'Stage numerical/schema/lifecycle review failed')
        require(all(sha(path) == digest for path,digest in frozen.items()), 'Frozen files changed during execution')
        print('STAGE_REVIEW '+str(root/(stage.upper()+'_REVIEW.json')), flush=True)
    except BaseException as error:
        save(execution/'FAILED_OR_INCOMPLETE.json', {'checks_pass': False, 'error': str(error),
            'completed_cases': completed, 'elapsed_s': time.monotonic()-start})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
