#!/usr/bin/env python3
"""Run immutable generic-particle cards with a fixed 16-physical-core budget on gnode1."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
LAYOUTS = ((16, 1), (8, 2), (4, 4), (2, 8))
TIMEOUT = 180


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def launch(binary, ranks, threads, directory):
    cmd = ['timeout', '--kill-after=10', str(TIMEOUT), 'mpirun', '--host', 'gnode1:16',
           '-np', str(ranks), '--map-by', 'slot:PE=' + str(threads), '--bind-to',
           'core', '--nooversubscribe', '--report-bindings', str(binary),
           str(directory / 'run.yaml')]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
               OMP_PROC_BIND='close', APRL_PROFILE_INTERVAL='0')
    save(directory / 'launch.json', {'command': cmd, 'cwd': str(directory),
         'environment': {key: env.get(key) for key in
                         ('OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_PLACES', 'OMP_PROC_BIND',
                          'APRL_PROFILE_INTERVAL', 'LD_LIBRARY_PATH', 'SLURM_JOB_ID',
                          'OMPI_MCA_pml', 'OMPI_MCA_btl', 'UCX_TLS')}})
    start = time.monotonic()
    with (directory / 'run.log').open('w') as stream:
        result = subprocess.run(cmd, cwd=directory, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    return result.returncode, time.monotonic() - start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('balanced', 'field-heavy', 'particle-heavy', 'holdout'), required=True)
    parser.add_argument('--mode', choices=('pilot', 'matrix'), required=True)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if platform.node() != 'gnode1' or os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') != '1':
        parser.error('explicitly authorized gnode1 CPU tests only')
    if args.repeats < 1 or args.repeats > 5:
        parser.error('repeats must be 1..5')
    case_path = WORKSPACE / 'tests/cpu_calibration/cards' / (args.case + '.yaml')
    generic_path = WORKSPACE / 'source/aprl/verification/performance/cases/generic_particles/run.py'
    inspector_path = WORKSPACE / 'source/aprl/verification/performance/run.py'
    generic = load_module('aprl_generic_16', generic_path)
    inspector = load_module('aprl_inspect_16', inspector_path)
    payload, case = generic.read_case(case_path)
    if case['mesh_cells'][2] < 32:
        parser.error('longitudinal mesh too small for 16 ranks')
    binaries = {label: WORKSPACE / 'build' / label / 'aprl'
                for label in ('baseline', 'candidate')}
    manifest = WORKSPACE / 'tests/cpu_calibration/prepared-manifest.sha256'
    subprocess.run(['sha256sum', '--check', '--strict', '--quiet', str(manifest)],
                   cwd=WORKSPACE, check=True)
    root_parent = WORKSPACE / 'results/cpu_calibration'
    root_parent.mkdir(exist_ok=True)
    import tempfile
    root = Path(tempfile.mkdtemp(prefix=f'{args.case}-{args.mode}-gnode1-{datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")}-', dir=root_parent))
    card_hash = hashlib.sha256(payload).hexdigest()
    (root / 'case.yaml').write_bytes(payload)
    rows = ([('baseline', 16, 1, 1), ('candidate', 16, 1, 1)] if args.mode == 'pilot'
            else [(label, ranks, threads, rep)
                  for ranks, threads in LAYOUTS
                  for label, count in (('baseline', 1), ('candidate', args.repeats))
                  for rep in range(1, count + 1)])
    report = {'status': 'running', 'pass': False, 'host': platform.node(),
              'case_name': args.case, 'mode': args.mode, 'card_sha256': card_hash,
              'case': {k: v for k, v in case.items() if k != 'physical_card'},
              'cpu_budget': 16, 'physical_cores': 16, 'step_timeout_s': TIMEOUT,
              'bindings_preflight': str(WORKSPACE / 'results/node_survey/cpu-preflight-gnode1-20261023T183105Z-cqz0Cb/PASS'),
              'binaries': {k: {'path': str(v), 'sha256': sha(v)} for k, v in binaries.items()},
              'sources': {str(p): sha(p) for p in (case_path, generic_path, inspector_path,
                         WORKSPACE / 'tests/cpu_calibration/run.py',
                         WORKSPACE / 'tests/cpu_calibration/run.sh', manifest)},
              'comparison_policy': {'same_rank': 'strict bitwise; gating',
                                    'cross_rank': 'not compared for gating',
                                    'rtol': generic.RTOL, 'scaled_atol': generic.SCALED_ATOL},
              'runs': [], 'comparisons': [], 'environment': {'loadavg_start': os.getloadavg(),
                  'python': sys.executable, 'mpi': subprocess.run(['mpirun', '--version'],
                  capture_output=True, text=True).stdout.splitlines()[0],
                  'slurm_job_id': os.environ.get('SLURM_JOB_ID')},
              'notes': ['Each run uses the exact same card bytes. Between cases, work sizes differ.',
                        'A workspace lock prevents this script running twice, but does not reserve the node.',
                        'Reported wall includes initialization; loop wall is the measured simulation loop.']}
    save(root / 'report.json', report)
    print('Result directory:', root, flush=True)
    references = {}
    try:
        for label, ranks, threads, repeat in rows:
            name = f'rep{repeat:02d}-{label}-r{ranks}-t{threads}'
            directory = root / name
            directory.mkdir()
            (directory / 'run.yaml').write_bytes(payload)
            run = {'name': name, 'directory': str(directory), 'label': label,
                   'ranks': ranks, 'threads': threads, 'repeat': repeat,
                   'profile_interval': 0, 'status': 'running', 'valid': False}
            report['runs'].append(run)
            save(root / 'report.json', report)
            if sha(directory / 'run.yaml') != card_hash:
                raise ValueError(name + ': card changed before launch')
            print('Running', name, flush=True)
            rc, elapsed = launch(binaries[label], ranks, threads, directory)
            run.update(returncode=rc, process_wall_s=elapsed)
            save(root / 'report.json', report)
            if rc:
                raise ValueError(name + f': process exit {rc}; see run.log')
            if sha(directory / 'run.yaml') != card_hash:
                raise ValueError(name + ': card changed during run')
            summary, arrays = generic.inspect(inspector, run, case)
            run.update(summary, valid=True, status='passed')
            if label == 'baseline':
                references[ranks] = (name, arrays)
            else:
                reference_name, reference = references[ranks]
                comparison = {'reference': reference_name, 'candidate': name,
                              **generic.compare(inspector, reference, arrays, True)}
                report['comparisons'].append(comparison)
                if not comparison['pass']:
                    raise ValueError(name + ': strict same-rank comparison failed')
            save(root / 'report.json', report)
        if sha(case_path) != card_hash or sha(root / 'case.yaml') != card_hash:
            raise ValueError('case card changed')
        subprocess.run(['sha256sum', '--check', '--strict', '--quiet', str(manifest)],
                       cwd=WORKSPACE, check=True)
        report['pass'] = True
        report['status'] = 'passed'
        (root / 'PASS').write_text('All completed runs validated and same-rank fields strictly bitwise equal.\n')
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = str(exc)
        print('FAIL:', exc, flush=True)
    report['environment']['loadavg_end'] = os.getloadavg()
    report['timing_summary'] = []
    for label, ranks, threads in sorted({(row['label'], row['ranks'], row['threads'])
                                          for row in report['runs'] if row['valid']}):
        selection = [r for r in report['runs'] if r['valid'] and
                     (r['label'], r['ranks'], r['threads']) == (label, ranks, threads)]
        report['timing_summary'].append({'label': label, 'ranks': ranks, 'threads': threads,
            'repeats': len(selection), 'reported_wall_s_median': statistics.median(r['reported_wall_s'] for r in selection),
            'reported_wall_s_min': min(r['reported_wall_s'] for r in selection),
            'reported_wall_s_max': max(r['reported_wall_s'] for r in selection),
            'loop_wall_s_median': statistics.median(r['loop_wall_s'] for r in selection)})
    save(root / 'report.json', report)
    (root / 'EXIT_CODE').write_text('0\n' if report['pass'] else '1\n')
    print('CPU calibration:', report['status'], root / 'report.json', flush=True)
    return 0 if report['pass'] else 1


if __name__ == '__main__':
    sys.exit(main())
