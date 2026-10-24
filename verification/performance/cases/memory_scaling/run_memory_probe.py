#!/usr/bin/env python3
"""Bounded, serialized gnode1 array-throughput and first-touch diagnostic."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.node() != 'gnode1' or os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') != '1':
        parser.error('Explicitly authorized gnode1 CPU diagnostic only')
    if os.environ.get('SLURM_JOB_ID'):
        parser.error('Use the scheduled workflow inside Slurm')
    w = args.workspace.resolve()
    lock = (w / 'results/node_survey/.cpu-calibration.lock').open('a')
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    if os.getloadavg()[0] > 2:
        raise RuntimeError('Host is busy; do not start this independent diagnostic')
    topology = subprocess.check_output(['lscpu', '-p=CPU,CORE,SOCKET,NODE'], text=True)
    physical = [tuple(map(int, line.split(','))) for line in topology.splitlines() if not line.startswith('#')]
    selected = [row for row in physical if row[0] < 16]
    assert len(selected) == len({(r[2], r[1]) for r in selected}) == 16
    args.output.mkdir(parents=True, exist_ok=False)
    binary = w / 'build/memory-scaling/memory_probe'
    source = w / 'source/aprl/verification/performance/cases/memory_scaling/memory_probe.cpp'
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    report = dict(status='running', pass_=False, host=platform.node(),
                  started_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  binary=str(binary), binary_sha256=digest(binary), source_sha256=digest(source),
                  topology=physical, load_start=os.getloadavg(), results=[],
                  notes=['Effective array bytes/s uses two reads and one write per element; not hardware-counter DRAM bandwidth.',
                         'Anonymous NUMA pages include allocator/runtime overhead; these are page-placement counts, not remote-access counts.',
                         'All threads use distinct physical cores. OMP_PROC_BIND=spread changes the available NUMA domains with thread count.',
                         'Same arrays and operation compare serial with parallel first-touch at a fixed thread count.',
                         'This memory diagnostic does not measure or change APRL physics or its production runtime.'])
    destination = args.output / 'report.json'
    def save():
        temporary = destination.with_suffix('.tmp')
        temporary.write_text(json.dumps(report, indent=2) + '\n')
        temporary.replace(destination)
    rows = [(size, threads, touch) for size in (64,512,4096,32768,262144)
            for threads in (1,2,4,8,16) for touch in ('serial','parallel')]
    random.Random(20261025).shuffle(rows)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='16', OMP_DYNAMIC='FALSE', OMP_PROC_BIND='spread',
               OMP_PLACES=','.join('{'+str(i)+'}' for i in range(16)))
    start = time.monotonic()
    save()
    try:
        for size, threads, touch in rows:
            if time.monotonic()-start > 180:
                raise RuntimeError('Three-minute diagnostic deadline reached')
            command = ['timeout','--kill-after=2','10','taskset','-c','0-15',str(binary),str(size),str(threads),touch]
            result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=15)
            if result.returncode:
                raise RuntimeError(f'Probe failed: {command}: {result.stderr}')
            row = json.loads(result.stdout)
            assert row['correct'] and row['threads'] == threads
            assert len(set(row['cpus'])) == threads
            report['results'].append(row)
            save()
        assert digest(binary) == report['binary_sha256'] and digest(source) == report['source_sha256']
        report.update(status='passed', pass_=True)
    except Exception as error:
        report.update(status='failed', error=str(error))
    report.update(elapsed_s=time.monotonic()-start, load_end=os.getloadavg())
    save()
    print(json.dumps(dict(status=report['status'], samples=len(report['results']), report=str(destination))))
    return 0 if report['pass_'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
