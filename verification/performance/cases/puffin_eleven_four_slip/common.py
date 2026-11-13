"""Frozen, independent-source Puffin eleven-segment campaign helpers."""
import hashlib
import json
import os
import signal
import subprocess
import time
from pathlib import Path

W = Path('/public/home/wfollowm/APRL_Workspace')
C = 299792458.0
E = 1.602176634e-19


def read(path): return json.loads(Path(path).read_text())
def save(path, value): Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
def need(value, message):
    if not value: raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''): h.update(block)
    return h.hexdigest()


def frozen(path):
    for name, expected in read(path).items(): need(sha(name) == expected, 'Frozen identity changed: ' + name)


def srun(nodes, ranks):
    return ['srun', '--mpi=pmix_v3', '--exact', '--exclusive', f'--nodes={nodes}',
            f'--ntasks={ranks}', f'--ntasks-per-node={ranks//nodes}', '--cpus-per-task=1',
            '--cpu-bind=cores', '--kill-on-bad-exit=1', '--wait=30']


def launch(command, cwd, log, limit, progress=False):
    """Hard wall limit, 600s no-output watchdog and MPI-free process cleanup."""
    log = Path(log); start = time.monotonic(); last = start; size = -1
    with log.open('x') as stream:
        proc = subprocess.Popen(command, cwd=cwd, stdin=subprocess.DEVNULL, stdout=stream,
                                stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while proc.poll() is None:
                time.sleep(5)
                now = time.monotonic(); current = log.stat().st_size
                if current != size: last = now; size = current
                if now-start > limit or (progress and now-last > 600):
                    raise TimeoutError('Wall/no-output watchdog: ' + str(log))
            need(proc.returncode == 0, 'Command failed: ' + str(log))
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try: proc.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait(timeout=10)
            save(str(log)+'.launch.json', dict(command=command, exit_code=proc.returncode,
                 seconds=time.monotonic()-start, job=os.environ.get('SLURM_JOB_ID')))


def dump_at(directory, kind, step):
    import h5py
    choices=[]
    for path in Path(directory).glob(f'main_{kind}_*.h5'):
        with h5py.File(path, 'r') as f:
            if int(f['runInfo'].attrs['iCsteps']) == step: choices.append(path)
    need(bool(choices), 'Missing closed physical step: ' + str(step))
    return min(choices, key=lambda p: int(p.stem.rsplit('_', 1)[1]))
