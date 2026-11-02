"""Preserved bounded launch and binding helpers from energy_transport_audit/run.py."""
import os,re,signal,subprocess,time
from prepare import require,save

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
