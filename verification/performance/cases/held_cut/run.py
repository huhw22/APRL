#!/usr/bin/env python3
"""Bounded, frozen-input held-inflow/causal-cut verification on gnode1."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import time

WORK = Path('/public/home/wfollowm/APRL_Workspace')
SOURCE = WORK / 'source/aprl'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')


def prepare(args):
    root = args.root.resolve()
    require(not root.exists(), 'Use a fresh result directory')
    require(1 <= args.threads <= 16, 'Physical CPU budget exceeded')
    binary = args.binary.resolve()
    checks = args.checks.resolve() if args.checks else None
    binding = WORK / ('build/cluster-tools/binding_probe' if args.environment == 'slurm'
                      else 'build/node-survey/binding_probe')
    for path in (binary, binding) + ((checks,) if checks else ()):
        require(path.is_file() and os.access(path, os.X_OK), 'Missing executable '+str(path))
    cases = []
    for specification in args.case or []:
        parts = specification.split(':')
        require(len(parts) >= 2 and parts[0].replace('-', '').replace('_', '').isalnum(),
                'Use name:argument:argument ...')
        require(parts[0] not in [x['name'] for x in cases], 'Duplicate case')
        cases.append({'name': parts[0], 'arguments': parts[1:]})
    root.mkdir(parents=True)
    files = sorted((SOURCE/'src').glob('*'))
    files += [SOURCE/'CMakeLists.txt', SOURCE/'verification/CMakeLists.txt']
    files += list((SOURCE/'verification').glob('held_cut*'))
    files += [SOURCE/'verification/moving_window_remap.h', SOURCE/'verification/moving_window_remap.cpp']
    files += list(Path(__file__).resolve().parent.glob('*'))
    files = [p for p in files if p.is_file()]
    save(root/'source-files.json', {str(p): sha(p) for p in files})
    for p in files:
        dest = root/'source'/p.relative_to(SOURCE)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    shutil.copy2(__file__, root/'run.py')
    shutil.copy2(Path(__file__).with_name('run.sbatch'), root/'run.sbatch')
    production = WORK/'results/elegant_field_only_20261027/production'
    for name in ('run.yaml', 'geometry-preflight.json', 'aprl.log'):
        shutil.copy2(production/name, root/('production-'+name))
    moment_sample = WORK/'results/held_cut_20261028/source-moment-sample.json'
    if moment_sample.exists():
        shutil.copy2(moment_sample, root/'source-moment-sample.json')
    (root/'source.diff').write_bytes(subprocess.check_output(
        ['git','-C',str(SOURCE),'diff','--binary','HEAD']))
    (root/'source-status.txt').write_text(subprocess.check_output(
        ['git','-C',str(SOURCE),'status','--short'], text=True))
    plan = {'format':'aprl-held-cut-campaign-v1', 'cases':cases,
            'binary':str(binary), 'checks':str(checks) if checks else None,
            'binding':str(binding), 'threads':args.threads, 'ranks':1,
            'environment':args.environment,
            'case_timeout_s':600 if args.environment == 'slurm' else 180,
            'total_timeout_s':1500 if args.environment == 'slurm' else 7200,
            'source_base':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),
            'scope':'Production constants, shortened generated bunch and run; no production switch, no speed claim'}
    save(root/'plan.json', plan)
    frozen = [p for p in root.rglob('*') if p.is_file()]
    frozen += [binary, binding] + ([checks] if checks else [])
    frozen += [p for p in (WORK/'build/cowan-accepted-r2-20261025/runtime-libs').glob('*') if p.is_file()]
    save(root/'frozen.json', {str(p):sha(p) for p in frozen})
    print(root/'plan.json')


def idle():
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            exe = Path(os.fsdecode((p/'cmdline').read_bytes().split(b'\0')[0])).name
            require(not (exe in ('aprl','mpirun','mpiexec','orted') or exe.startswith('aprl_')),
                    'Existing compute process '+p.name+' '+exe)
        except (OSError, ValueError):
            pass
    def read():
        x = list(map(int,Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(x), x[3]
    a,b = read()
    time.sleep(2)
    c,d = read()
    fraction = (d-b)/(c-a)
    require(fraction >= .95, 'Host busy: idle fraction '+str(fraction))
    return fraction


def launch(command, directory, env, limit, log, record):
    started = time.monotonic()
    command = ['timeout','--signal=TERM','--kill-after=5',str(limit)] + command
    result = {'command':command, 'started_unix':time.time(), 'timeout_s':limit}
    proc = None
    try:
        with log.open('x') as stream:
            proc = subprocess.Popen(command,cwd=directory,env=env,stdin=subprocess.DEVNULL,
                                    stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            result['pid'] = proc.pid
            result['exit_code'] = proc.wait(timeout=limit+10)
            require(result['exit_code'] == 0,'Failed: '+str(log))
    except BaseException as error:
        result['error'] = str(error)
        raise
    finally:
        if proc and proc.poll() is None:
            os.killpg(proc.pid,signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL)
                proc.wait(timeout=5)
        result['seconds'] = time.monotonic()-started
        save(record,result)


def run(args):
    root = args.root.resolve()
    require(not (root/'STARTED.json').exists(),'Already started; use fresh directory')
    for name,digest in json.loads((root/'frozen.json').read_text()).items():
        require(sha(name) == digest,'Frozen input changed: '+name)
    plan = json.loads((root/'plan.json').read_text())
    scheduled = plan.get('environment') == 'slurm'
    if scheduled:
        require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES','0')) == 1,
                'Requires one allocated Slurm node')
        require(int(os.environ.get('SLURM_CPUS_PER_TASK','0')) == 128,
                'Reserve exclusive 128 CPUs; steps use the recorded smaller budget')
    else:
        require(socket.gethostname().split('.')[0] == 'gnode1' and
                os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1' and
                not os.environ.get('SLURM_JOB_ID'),'Requires authorized gnode1 CPU session')
    require(plan['ranks'] == 1 and 1 <= plan['threads'] <= 16 and
            plan['case_timeout_s'] == (600 if scheduled else 180) and
            plan['total_timeout_s'] == (1500 if scheduled else 7200),'Budget changed')
    with (WORK/'results/node_survey/.cpu-calibration.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        info = {'host':socket.gethostname(),'unix_time':time.time(),
                'idle_fraction':None if scheduled else idle(),
                'slurm_job_id':os.environ.get('SLURM_JOB_ID'),
                'allocated_cpus':128 if scheduled else 16}
        save(root/'STARTED.json',info)
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS=str(plan['threads']),OMP_DYNAMIC='FALSE',OMP_PLACES='cores',
                   OMP_PROC_BIND='close',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
                   NUMEXPR_NUM_THREADS='1',APRL_PROFILE_INTERVAL='0')
        launcher = (['srun','--mpi='+os.environ.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1',
                     '--ntasks=1','--ntasks-per-node=1','--cpus-per-task='+str(plan['threads']),
                     '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'] if scheduled else
                    ['mpiexec','-n','1','--map-by','slot:PE='+str(plan['threads']),
                     '--bind-to','core','--nooversubscribe'])
        deadline = time.monotonic()+plan['total_timeout_s']
        try:
            launch(launcher+[plan['binding'],'1',str(plan['threads'])],root,env,45,
                   root/'binding.log',root/'binding-launch.json')
            require('[binding] result=PASS failures=0' in (root/'binding.log').read_text(),
                    'CPU binding failed')
            if plan['checks']:
                directory = root/'checks'
                directory.mkdir()
                launch(launcher+[plan['checks'],str(directory/'metrics.json')],directory,env,180,
                       directory/'probe.log',directory/'launch.json')
            cases = []
            for case in plan['cases']:
                require(time.monotonic()+plan['case_timeout_s']+10 < deadline,'Campaign deadline reached')
                directory = root/case['name']
                directory.mkdir()
                launch(launcher+[plan['binary'],str(directory)]+case['arguments'],directory,env,plan['case_timeout_s'],
                       directory/'probe.log',directory/'launch.json')
                metric = directory/'metrics.json'
                require(metric.is_file(),'Missing metrics '+str(directory))
                parsed = json.loads(metric.read_text())
                require(parsed.get('completed') is True,'Incomplete case '+case['name'])
                cases.append({'name':case['name'],'metrics_sha256':sha(metric)})
                print('FINISHED '+case['name'],flush=True)
            save(root/'COMPLETED.json',{'status':'EXECUTION_COMPLETED_SCIENCE_REVIEW_REQUIRED',
                                      'unix_time':time.time(),'cases':cases})
        except BaseException as error:
            save(root/'FAILED.json',{'unix_time':time.time(),'error':str(error)})
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='action',required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--checks',type=Path)
    p.add_argument('--threads',type=int,default=16)
    p.add_argument('--environment',choices=('gnode1','slurm'),default='gnode1')
    p.add_argument('--case',action='append')
    p = sub.add_parser('run')
    p.add_argument('--root',type=Path,required=True)
    options = parser.parse_args()
    (prepare if options.action == 'prepare' else run)(options)
