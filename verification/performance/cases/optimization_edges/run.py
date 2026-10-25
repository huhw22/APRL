#!/usr/bin/env python3
"""Small 4 MPI x 4 OpenMP scientific equivalence gates; never a speed benchmark."""
import argparse
import copy
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
HELPER = HERE / 'science_helpers.py'
if not HELPER.is_file():
    HELPER = HERE.parent / 'full_scale_cowan/run.py'
spec = importlib.util.spec_from_file_location('edge_science', HELPER)
science = importlib.util.module_from_spec(spec)
spec.loader.exec_module(science)
require, digest, save = science.require, science.digest, science.save
CONTROL = WORKSPACE / 'build/cowan-accepted-r2-20261025/aprl'
CONTROL_REVISION = '77f1ef95267d037dcd629584e42aedc1b04668f3'
CONTROL_SHA = 'dc537da20eacc15648cc59360a1f61bdc9ef1176950ac0cdfda21714e2ba58bb'
CASE_NAMES = ('held-crossrank', 'held-retirement', 'cpml-exit', 'outer-exit')


def compiler_flags(binary, revision):
    result = {}
    paths = []
    for target in ('direct_eb_core', 'aprl'):
        path = binary.parent / 'CMakeFiles' / (target + '.dir') / 'flags.make'
        source = path.read_text()
        if target == 'direct_eb_core':
            matches = re.findall(r'APRL_SOURCE_REVISION=\\"([0-9a-f]{40})\\"', source)
            require(matches == [revision], 'Compiled revision does not match core flags')
        result[target] = {
            key: re.sub(r'APRL_SOURCE_REVISION=\\"[0-9a-f]{40}\\"',
                        'APRL_SOURCE_REVISION=<identity>', value)
            for key, value in re.findall(r'^(CXX_[A-Z]+) = (.*)$', source, re.M)}
        require('CXX_FLAGS' in result[target] and '-fopenmp' in result[target]['CXX_FLAGS'],
                'Missing real OpenMP compiler flags')
        paths.append(path)
    return result, paths


def required_evidence(directory, binary):
    directory = Path(directory).resolve()
    require((directory / 'PASS').is_file() and
            (directory / 'EXIT_CODE').read_text().strip() == '0' and
            not (directory / 'FAIL').exists(), 'Required checks did not pass: ' + str(directory))
    hashes = directory / 'build.sha256'
    matches = [line.split('  ', 1)[0] for line in hashes.read_text().splitlines()
               if line.endswith('  ' + str(binary))]
    require(matches == [digest(binary)], 'Required checks are not tied to this binary')
    return [directory / name for name in ('PASS', 'EXIT_CODE', 'build.sha256')]


def boundary_card(cpml):
    # Eight coincident, freely drifting macro-particles start 0.05 um inside
    # the selected +x face. The early hit is checked in the actual run log.
    return {
        'units': {'length': 'micrometer', 'time': 'femtosecond'},
        'runtime': {'mode': 'throughput', 'resource_monitor': {'enabled': False}},
        'output': {'overwrite': False, 'manifest': 'run-manifest.yaml'},
        'mesh': {'field_solver': 'cowan-z', 'current_filter': 'none',
                 'cells': [12, 12, 80], 'cell_size': [0.5, 0.5, 0.1],
                 'center': [0, 0, 0], 'duration': 30, 'boost_gamma': 2,
                 'particle_steps_per_undulator_period': 1,
                 'ballistic_outside_transverse': False, 'self_field_feedback': False},
        'boundary': ({'type': 'cpml', 'cells': [2, 2, 4], 'polynomial_order': 3,
                      'target_reflection': 1e-8, 'kappa_max': 8, 'alpha_fraction': 0}
                     if cpml else {'type': 'pec'}),
        'initial_self_field': {'enabled': False},
        'beam': {'reference': {'initial_center_z': -0.5}, 'input': {
            'type': 'generated-gaussian', 'electrons': 8, 'macroparticles': 8,
            'gamma': 4, 'direction': [0.2, 0, 0.9797958971132712],
            'center': [1.95 if cpml else 2.95, 0, 0],
            'sigma_position': [0, 0, 0], 'sigma_proper_velocity': [0, 0, 0],
            'random_seed': 260126}},
        'sources': {'incident_waves': [], 'magnetic_elements': []},
        'detectors': {'enabled': True, 'directory': 'detectors',
                      'field_planes': [{'name': 'field', 'z': 1, 'rhythm': 1,
                                        'buffer_samples': 4, 'compression': 0}],
                      'particle_planes': [{'name': 'after-exit', 'z': 1,
                                           'buffer_records': 16, 'compression': 0}]},
        'stop': {'mode': 'reference-center-z', 'z': 2},
        'trajectory': {'enabled': False},
        'energy_ledger': {'enabled': True, 'directory': 'ledger',
                          'filename': 'energy-ledger.h5', 'sample_interval_steps': 1,
                          'buffer_records': 8, 'compression': 0,
                          'warning_relative_tolerance': 0.05}}


def prepare(args):
    root = args.root.resolve()
    require(not root.exists(), 'Use a fresh result directory')
    variants = json.loads(args.variants.read_text())
    require(isinstance(variants, dict) and 2 <= len(variants) <= 5 and
            'control' in variants,
            'Provide accepted control plus one to four implementation candidates')
    frozen = []
    metadata = {}
    real_flags = {}
    for name, item in variants.items():
        require(re.fullmatch(r'[a-z0-9-]+', name), 'Invalid variant name')
        binary = Path(item['binary']).resolve()
        require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable')
        require(str(binary).startswith(str(WORKSPACE / 'build') + '/'), 'Use workspace builds')
        require(re.fullmatch(r'[0-9a-f]{40}', item['revision']), 'Use full compiled revision')
        if name == 'control':
            require(binary == CONTROL and item['revision'] == CONTROL_REVISION and
                    digest(binary) == CONTROL_SHA,
                    'Control must be the accepted r2 executable')
        else:
            frozen.extend(required_evidence(item['required_directory'], binary))
        metadata[name] = science.build_metadata(binary)
        real_flags[name], flags_paths = compiler_flags(binary, item['revision'])
        frozen.extend(flags_paths)
        item['binary'] = str(binary)
        item['sha256'] = digest(binary)
        frozen.extend([binary, Path(metadata[name]['cache'])])
        library_dir = binary.parent / 'runtime-libs'
        if library_dir.is_dir():
            frozen.extend(path for path in library_dir.iterdir() if path.is_file())
    require(all(m['selected'] == metadata['control']['selected'] for m in metadata.values()),
            'Compiler, MPI, HDF5 or OpenMP flags differ')
    require(all(flags == real_flags['control'] for flags in real_flags.values()),
            'Actual generated compiler flags differ beyond source revision')
    require(metadata['control']['selected']['APRL_ENABLE_OPENMP'] == 'ON', 'OpenMP must be enabled')
    root.mkdir(parents=True)
    scripts = root / 'scripts'
    scripts.mkdir()
    for name in ('run.py', 'held.yaml.in', 'README.txt'):
        shutil.copy2(HERE / name, scripts / name)
    shutil.copy2(HELPER, scripts / 'science_helpers.py')
    inputs = root / 'inputs'
    inputs.mkdir()
    held_source = CONTROL.parent / 'verification/work/hybrid-held/particles.h5'
    retirement_source = WORKSPACE / 'results/elegant_retirement_20261024/compact-observed-window/preflight-r2/functional'
    shutil.copy2(held_source, inputs / 'held.h5')
    shutil.copy2(retirement_source / 'particles.h5', inputs / 'retirement.h5')
    held = yaml.safe_load((HERE / 'held.yaml.in').read_text())
    held['beam']['input']['file'] = str(inputs / 'held.h5')
    retirement = yaml.safe_load((retirement_source / 'run.yaml').read_text())
    retirement['beam']['input']['file'] = str(inputs / 'retirement.h5')
    require(retirement['mesh']['cells'] == [16, 16, 4096] and
            retirement['particle_retirement']['enabled'] and
            retirement['mesh']['ballistic_outside_transverse'] is False,
            'Retirement reference fixture changed')
    cards = {'held-crossrank': held, 'held-retirement': retirement,
             'cpml-exit': boundary_card(True), 'outer-exit': boundary_card(False)}
    # Check inputs without running a generator or simulator.
    for filename, count in [('held.h5', 10410), ('retirement.h5', 512)]:
        with h5py.File(inputs / filename) as file:
            require(len(file['particles/records']) == count, 'Input particle cohort changed')
    for case, configuration in cards.items():
        configuration['output']['overwrite'] = False
        for name in variants:
            directory = root / case / name
            directory.mkdir(parents=True)
            redirected = science.redirect(configuration, directory)
            require(science.physics(redirected) == science.physics(configuration),
                    'Case physics differs between variants')
            (directory / 'run.yaml').write_text(yaml.safe_dump(redirected, sort_keys=False))
    probe = WORKSPACE / 'build/node-survey/binding_probe'
    require(probe.is_file(), 'Missing accepted binding probe')
    frozen.append(probe)
    plan = {'format': 'aprl-edge-acceptance-v1', 'variants': variants, 'toolchains': metadata,
            'compiler_flags_without_revision': real_flags,
            'cases': list(CASE_NAMES), 'ranks': 4, 'threads_per_rank': 4,
            'host': 'gnode1', 'timeout_per_run_s': 180, 'timeout_total_s': 1800,
            'binding_probe': str(probe),
            'reference_revision': '856e172',
            'scope': 'Small correctness fixtures only; no performance or radiation claim',
            'uncovered': ['Exact or tolerance-near CPML/outer/retirement/interface event ties']}
    save(root / 'plan.json', plan)
    frozen.extend(path for path in root.rglob('*') if path.is_file())
    save(root / 'frozen.json', {str(path): digest(path) for path in sorted(set(frozen))})
    print(root / 'plan.json')


def verify(root):
    for filename, sha in json.loads((root / 'frozen.json').read_text()).items():
        require(digest(Path(filename)) == sha, 'Frozen input changed: ' + filename)
    return json.loads((root / 'plan.json').read_text())


def environment(binary):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='4', OMP_DYNAMIC='FALSE', OMP_PLACES='cores',
               OMP_PROC_BIND='close', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', APRL_PROFILE_INTERVAL='0')
    library_dir = Path(binary).parent / 'runtime-libs'
    if library_dir.is_dir():
        env['LD_LIBRARY_PATH'] = str(library_dir) + ':' + env.get('LD_LIBRARY_PATH', '')
    return env


def idle_host_snapshot():
    """Same bounded idle policy as the accepted fixed-MPI harness."""
    other = []
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            name = Path(os.fsdecode((entry / 'cmdline').read_bytes().split(b'\0')[0])).name
            if name in ('aprl', 'mpirun', 'mpiexec', 'orted'):
                other.append({'pid': int(entry.name), 'executable': name})
        except (OSError, ValueError):
            continue
    def counters():
        values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(values), values[3]
    before, idle_before = counters()
    time.sleep(2)
    after, idle_after = counters()
    require(after > before, 'Could not measure current host idleness')
    fraction = (idle_after - idle_before) / (after - before)
    snapshot = {'loadavg': os.getloadavg(), 'sample_seconds': 2,
                'cpu_idle_fraction': fraction, 'observed_compute_processes': other}
    require(not other and fraction >= .95, 'gnode1 has active work: ' + json.dumps(snapshot))
    return snapshot


def launch(command, directory, env, deadline, timeout=180):
    budget = min(timeout, int(deadline - time.monotonic()))
    require(budget > 0, 'Total campaign deadline expired')
    command = ['timeout', '--signal=TERM', '--kill-after=10', str(budget)] + command
    save(directory / 'command.json', {'argv': command, 'cwd': str(directory),
                                     'environment': {k: env.get(k) for k in (
                                         'OMP_NUM_THREADS', 'OMP_DYNAMIC', 'OMP_PLACES',
                                         'OMP_PROC_BIND', 'LD_LIBRARY_PATH')}})
    with (directory / 'aprl.log').open('x') as output:
        process = subprocess.Popen(command, cwd=directory, env=env,
                                   stdout=output, stderr=subprocess.STDOUT)
        save(directory / 'LAUNCHED.json', {'pid': process.pid, 'host': socket.gethostname(),
                                           'started_unix_s': time.time()})
        code = process.wait(timeout=budget + 20)
    (directory / 'EXIT_CODE').write_text(str(code) + '\n')
    require(code == 0, 'Simulation/binding failed or timed out: ' + str(directory))


def check_case(root, case, name, entry):
    directory = root / case / name
    require((directory / 'EXIT_CODE').read_text().strip() == '0', 'Missing successful exit')
    log = (directory / 'aprl.log').read_text(errors='replace')
    for marker in ('Compute threads per MPI rank: 4;',
                   'Direct SI E/B simulation active;', 'Configured stop reached:'):
        require(marker in log, 'Missing real simulator marker: ' + marker)
    require('Simulation failure' not in log and science.EXHAUSTED not in log,
            'Simulation did not stop physically')
    manifest = yaml.safe_load((directory / 'run-manifest.yaml').read_text())
    require(manifest['mpi_size'] == 4 and manifest['source_revision'] == entry['revision'],
            'Wrong launched rank count or compiled revision')
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    require(yaml.safe_load(manifest['configuration_yaml']) == card, 'Launched card differs')
    outputs = science.inspect_outputs(directory, True)
    with h5py.File(directory / 'ledger/energy-ledger.h5') as file:
        records = file['energy_ledger/records'][:]
    require(len(records) > 2, 'No actual scientific steps')
    coverage = {}
    if case == 'held-crossrank':
        hit = re.search(r'cross-rank activation pushes=(\d+)', log)
        require(hit and int(hit[1]) > 0, 'Held fixture did not cross a rank at release')
        coverage['cross_rank_activation_pushes'] = int(hit[1])
        for detector in ('activation', 'downstream'):
            with h5py.File(directory / 'detectors' / (detector + '.h5')) as file:
                require(len(file['particle_plane/records']) == 10410, 'Held particle cohort incomplete')
    elif case == 'held-retirement':
        require('Experimental particle retirement: entries=512 ' in log and
                re.search(r'Experimental particle retirement:.*completed=512\b', log),
                'Retirement entry or taper completion not exercised')
        require(records[-1]['active_macroparticles'] == 0 and
                records[-1]['removed_macroparticles'] == 512, 'Retirement cohort incomplete')
        empty = np.flatnonzero(records['active_macroparticles'] == 0)
        require(len(empty) and len(records) - int(empty[0]) >= 20,
                'Missing post-retirement field propagation')
        require(manifest['resolved_beam_reference']['startup_mode'] == 'held', 'Held startup missing')
        coverage['retired_particles'] = 512
    else:
        label = 'CPML physical-domain exits:' if case == 'cpml-exit' else 'Direct outer-domain exits:'
        require(re.search(re.escape(label) + r'[^\n]*x\+=8\b', log),
                'Selected boundary did not remove all eight particles')
        require(records[-1]['active_macroparticles'] == 0 and
                records[-1]['removed_macroparticles'] == 8, 'Boundary cohort incomplete')
        coverage['boundary_exits'] = 8
    (directory / 'STARTED_AND_COMPLETED').write_text('Actual simulator, manifest and scientific steps verified.\n')
    result = {'pass': True, 'case': case, 'variant': name, 'coverage': coverage,
              'outputs': outputs, 'ledger_samples': len(records)}
    save(directory / 'acceptance.json', result)
    return result


def run(args):
    root = args.root.resolve()
    plan = verify(root)
    require(socket.gethostname().split('.')[0] == 'gnode1' and
            os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') == '1' and
            not os.environ.get('SLURM_JOB_ID'), 'Run only on the authorized gnode1 CPU test host')
    lock_path = WORKSPACE / 'results/node_survey/.cpu-calibration.lock'
    if args.lock_fd is None:
        lock_file = lock_path.open('a')
        lock_fd = lock_file.fileno()
    else:
        lock_fd = args.lock_fd
        require(Path(os.readlink('/proc/self/fd/' + str(lock_fd))).resolve() == lock_path.resolve(),
                'Inherited lock descriptor points to another file')
    fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    save(root / 'initial-resources.json', idle_host_snapshot())
    (root / 'RUN_STARTED').mkdir()
    deadline = time.monotonic() + plan['timeout_total_s']
    for name, entry in plan['variants'].items():
        env = environment(entry['binary'])
        linked = subprocess.run(['ldd', '-r', entry['binary']], env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
        (root / ('ldd-' + name + '.txt')).write_text(linked.stdout)
        require(linked.returncode == 0 and 'not found' not in linked.stdout and
                'undefined symbol' not in linked.stdout, 'Runtime dependency failure: ' + name)
    mpi = ['mpirun', '--host', 'gnode1:16', '-np', '4', '--map-by', 'slot:PE=4',
           '--bind-to', 'core', '--nooversubscribe', '--report-bindings']
    binding = root / 'binding'
    binding.mkdir()
    launch(mpi + [plan['binding_probe'], '4', '4'], binding,
           environment(plan['variants']['control']['binary']), deadline, 60)
    require('[binding] result=PASS failures=0' in (binding / 'aprl.log').read_text(), 'Binding failed')
    reports, comparisons = [], []
    names = ['control'] + [name for name in plan['variants'] if name != 'control']
    for case in plan['cases']:
        for name in names:
            verify(root)
            entry = plan['variants'][name]
            directory = root / case / name
            launch(mpi + [entry['binary'], str(directory / 'run.yaml')], directory,
                   environment(entry['binary']), deadline, plan['timeout_per_run_s'])
            reports.append(check_case(root, case, name, entry))
            if name != 'control':
                result = science.compare_outputs(root / case / 'control', directory)
                result.update(case=case, variant=name, mode='same 4 MPI x 4 OpenMP scientific bitwise')
                comparisons.append(result)
                save(root / 'comparisons.json', comparisons)
    save(root / 'acceptance.json', {'pass': True, 'reports': reports, 'comparisons': comparisons,
                                  'uncovered': plan['uncovered'], 'scope': plan['scope']})
    (root / 'PASS').write_text(
        f"All {len(plan['cases'])} small cases matched accepted r2 for "
        f"{len(names) - 1} candidate(s) at fixed 4 MPI x 4 OpenMP "
        f"({len(reports)} completed runs).\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'run', 'verify'))
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--variants', type=Path)
    parser.add_argument('--lock-fd', type=int)
    args = parser.parse_args()
    if args.mode == 'prepare':
        require(args.variants is not None, 'prepare requires --variants')
        prepare(args)
    elif args.mode == 'run':
        try:
            run(args)
        except BaseException as error:
            if args.root.is_dir():
                (args.root / 'FAIL').write_text(type(error).__name__ + ': ' + str(error) + '\n')
            raise
    else:
        verify(args.root.resolve())
        print('Frozen inputs verified; no computation executed')


if __name__ == '__main__':
    main()
