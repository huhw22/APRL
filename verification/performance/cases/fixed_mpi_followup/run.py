#!/usr/bin/env python3
"""Fixed 4 MPI x 4 threads, same-physics implementation comparison."""
import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import time

import yaml

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
CONTROL_SHA = 'dc537da20eacc15648cc59360a1f61bdc9ef1176950ac0cdfda21714e2ba58bb'
CONTROL_REVISION = '77f1ef95267d037dcd629584e42aedc1b04668f3'
CONTROL_BUILD = Path('build/cowan-accepted-r2-20261025')
REQUIRED_TEST_COUNT = 72
CONTROL_RUNTIME_LIBS = {
    'libfftw3.so.3': 'b63413c26b5ae054d25df7899e5937a342e3d676254738f2fb891c8260296150',
    'libfftw3_threads.so.3': 'e606bf0c0fec66837c4ab36530701b5bedfd6fcbbb3a03bec7b2e77d4ee228ea',
}
CASE_NAMES = ('grid016', 'grid128', 'particles262144')
LAYOUT = (4, 4, 'compact')
REPEATS = 5
SEED = 20261025
CACHE_KEYS = (
    'CMAKE_BUILD_TYPE', 'CMAKE_CXX_COMPILER', 'CMAKE_C_COMPILER',
    'CMAKE_CXX_FLAGS', 'CMAKE_CXX_FLAGS_RELEASE', 'CMAKE_C_FLAGS',
    'CMAKE_C_FLAGS_RELEASE', 'MPI_CXX_COMPILER', 'HDF5_C_COMPILER_EXECUTABLE',
    'HDF5_C_LIBRARY_hdf5', 'APRL_REQUIRE_PARALLEL_HDF5', 'APRL_ENABLE_OPENMP',
    'OpenMP_CXX_FLAGS', 'OpenMP_CXX_LIB_NAMES', 'OpenMP_gomp_LIBRARY',
    'yaml-cpp_DIR')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def helper_path():
    here = Path(__file__).resolve().parent
    frozen = here / 'memory_driver.py'
    return frozen if frozen.exists() else here.parent / 'memory_scaling/run.py'


memory = module('fixed_mpi_memory_helpers', helper_path())
save, digest = memory.save, memory.digest
# The inherited binding routine must exercise this single, fixed layout only.
memory.LAYOUTS = [LAYOUT]


def build_metadata(binary):
    cache = binary.parent / 'CMakeCache.txt'
    if not cache.is_file():
        raise ValueError('Binary must have its build CMakeCache.txt beside it: ' + str(binary))
    values = {}
    for line in cache.read_text().splitlines():
        match = re.match(r'([^#/:][^:]*):[^=]+=(.*)$', line)
        if match:
            values[match[1]] = match[2]
    missing = [key for key in CACHE_KEYS if key not in values]
    if missing:
        raise ValueError('Incomplete compiler/dependency provenance: ' + str(missing))
    return dict(cache=str(cache), cache_sha256=digest(cache),
                selected={key: values[key] for key in CACHE_KEYS})


def compiled_revision_evidence(binary):
    """Read the revision recorded in this Makefile build's core compile flags."""
    flags = binary.parent / 'CMakeFiles/direct_eb_core.dir/flags.make'
    if not flags.is_file():
        raise ValueError('Compiled revision flags are missing: ' + str(flags))
    text = flags.read_text()
    matches = re.findall(r'-DAPRL_SOURCE_REVISION=\\"([0-9a-f]{40})\\"(?=\s|$)', text)
    if len(matches) != 1 or text.count('-DAPRL_SOURCE_REVISION=') != 1:
        raise ValueError('Expected one full APRL_SOURCE_REVISION in ' + str(flags))
    return dict(flags=str(flags), flags_sha256=digest(flags), compiled_revision=matches[0])


def required_evidence(pass_file, trial):
    """Bind a completed 72-test suite to the exact executable and build cache."""
    root = pass_file.parent
    exit_file = root / 'EXIT_CODE'
    hashes_file = root / 'build.sha256'
    driver_file = root / 'driver.log'
    tests_file = root / 'Testing/Temporary/LastTest.log'
    cache_copy = root / 'CMakeCache.txt'
    evidence_files = (pass_file, exit_file, hashes_file, driver_file, tests_file, cache_copy)
    for file in evidence_files:
        if not file.is_file():
            raise ValueError('Required-regression evidence is missing: ' + str(file))
    if exit_file.read_text().strip() != '0':
        raise ValueError('Required-regression EXIT_CODE is not zero')
    if (root / 'FAIL').exists():
        raise ValueError('Required-regression directory contains a FAIL marker')
    records = []
    for line in hashes_file.read_text().splitlines():
        match = re.fullmatch(r'([0-9a-f]{64}) [ *](.+)', line)
        if not match or not Path(match[2]).is_absolute():
            raise ValueError('Malformed or non-absolute required build.sha256 record')
        records.append((match[1], match[2]))
    verified = {}
    for target in (trial, trial.parent / 'CMakeCache.txt'):
        matches = [sha for sha, path in records if path == str(target)]
        if len(matches) != 1:
            raise ValueError('Required build.sha256 must contain exactly one absolute record for ' + str(target))
        if matches[0] != digest(target):
            raise ValueError('Required-regression SHA256 mismatch: ' + str(target))
        verified[str(target)] = matches[0]
    if digest(cache_copy) != verified[str(trial.parent / 'CMakeCache.txt')]:
        raise ValueError('Required-regression saved build cache differs from the tested build')
    driver = driver_file.read_text(errors='replace')
    summaries = re.findall(r'^100% tests passed, 0 tests failed out of ([0-9]+)$', driver, re.M)
    if summaries != [str(REQUIRED_TEST_COUNT)]:
        raise ValueError('Required-regression driver must report all 72 tests passed')
    passed = re.findall(r'^\s*([0-9]+)/([0-9]+) Test\s+#\s*[0-9]+: (required\.[^ ]+) +\.+ +Passed +[0-9.]+ sec$', driver, re.M)
    if len(passed) != REQUIRED_TEST_COUNT or {int(row[0]) for row in passed} != set(range(1, REQUIRED_TEST_COUNT + 1)) or \
            any(int(row[1]) != REQUIRED_TEST_COUNT for row in passed) or len({row[2] for row in passed}) != REQUIRED_TEST_COUNT:
        raise ValueError('Required-regression driver does not contain 72 unique successful test records')
    tests = tests_file.read_text(errors='replace')
    logged = re.findall(r'^([0-9]+)/([0-9]+) Test: (required\.[^\n]+)$', tests, re.M)
    if logged != passed or len(re.findall(r'^Test (?:Passed\.|Pass Reason:)$', tests, re.M)) != REQUIRED_TEST_COUNT or \
            re.search(r'^Test Failed', tests, re.M) or not re.search(r'^End testing:', tests, re.M):
        raise ValueError('Required-regression detailed test log is incomplete or disagrees with the driver')
    names = {row[2] for row in logged}
    if not {'required.runtime.duration-close', 'required.runtime.duration-barrier-watchdog'}.issubset(names):
        raise ValueError('Required-regression duration-close/watchdog checks are missing')
    # LastTestsFailed.log can survive an earlier CTest invocation. The complete,
    # matching current driver and LastTest.log above are authoritative.
    files = {str(file): digest(file) for file in evidence_files}
    return dict(pass_path=str(pass_file), sha256=files[str(pass_file)],
                trial_path=str(trial), trial_sha256=verified[str(trial)],
                cache_sha256=verified[str(trial.parent / 'CMakeCache.txt')],
                tests_passed=REQUIRED_TEST_COUNT, files=files)


def paired_order(case_names, repeat):
    cases = list(case_names)
    random.Random(SEED + repeat).shuffle(cases)
    for index, case in enumerate(cases):
        variants = ('control', 'trial') if (index + repeat) % 2 else ('trial', 'control')
        for variant in variants:
            yield case, variant


def prepare(args):
    workspace, root = args.workspace.resolve(), args.root.resolve()
    control, trial = args.control.resolve(), args.trial.resolve()
    if control == trial:
        raise ValueError('Control and trial must be separate executables')
    if digest(control) != CONTROL_SHA:
        raise ValueError('Control must be the accepted Cowan executable (dc537da / compiled 77f1ef9)')
    if not re.fullmatch(r'[0-9a-f]{40}', args.trial_revision):
        raise ValueError('--trial-revision must be the full compiled source revision')
    trial_revision = compiled_revision_evidence(trial)
    if args.trial_revision != trial_revision['compiled_revision']:
        raise ValueError('--trial-revision differs from APRL_SOURCE_REVISION in ' + trial_revision['flags'])
    metadata = {name: build_metadata(binary) for name, binary in (('control', control), ('trial', trial))}
    if metadata['control']['selected'] != metadata['trial']['selected']:
        differences = {key: [metadata['control']['selected'][key], metadata['trial']['selected'][key]]
                       for key in CACHE_KEYS
                       if metadata['control']['selected'][key] != metadata['trial']['selected'][key]}
        raise ValueError('Compiler/flags/MPI/HDF5/OpenMP toolchain differs: ' + json.dumps(differences))
    if metadata['control']['selected']['APRL_ENABLE_OPENMP'] != 'ON' or \
            metadata['control']['selected']['APRL_REQUIRE_PARALLEL_HDF5'] != 'ON':
        raise ValueError('OpenMP and required parallel HDF5 must remain enabled')
    runtime_directory = control.parent / 'runtime-libs'
    runtime_files = {str(runtime_directory / name): expected
                     for name, expected in CONTROL_RUNTIME_LIBS.items()}
    for file, expected_sha in runtime_files.items():
        if not Path(file).is_file() or digest(file) != expected_sha:
            raise ValueError('Accepted FFTW runtime library is missing or changed: ' + file)
    runtime = dict(directory=str(runtime_directory), files=runtime_files)
    required_pass = args.required_pass.resolve()
    required = required_evidence(required_pass, trial)
    source = args.source_campaign.resolve()
    source_plan = json.loads((source / 'plan.json').read_text())
    source_frozen = json.loads((source / 'frozen.json').read_text())
    if not (source / 'PASS').exists():
        raise ValueError('Source cards must come from the completed workload campaign')
    expected = {item['name']: item for item in source_plan['cases']}
    for name in CASE_NAMES:
        file = source / 'cards' / (name + '.yaml')
        if digest(file) != expected[name]['sha256'] or digest(file) != source_frozen[str(file)]:
            raise ValueError('Source card changed: ' + str(file))
    probe = workspace / 'build/node-survey/binding_probe'
    if not probe.is_file():
        raise ValueError('Prepared binding probe missing')
    root.mkdir(parents=True, exist_ok=False)
    scripts, cards = root / 'scripts', root / 'cards'
    scripts.mkdir(); cards.mkdir()
    here = Path(__file__).resolve().parent
    for name in ('run.py', 'run.sh', 'README.md'):
        shutil.copy2(here / name, scripts / name)
    shutil.copy2(helper_path(), scripts / 'memory_driver.py')
    for name in ('generic_inspector.py', 'basic_inspector.py'):
        original = source / 'scripts' / name
        if digest(original) != source_frozen[str(original)]:
            raise ValueError('Source inspector changed: ' + str(original))
        shutil.copy2(original, scripts / name)
    cases = []
    for name in CASE_NAMES:
        file = cards / (name + '.yaml')
        shutil.copy2(source / 'cards' / file.name, file)
        record = dict(expected[name]); record.update(card=str(file), sha256=digest(file))
        cases.append(record)
    binaries = {name: dict(path=str(binary), sha256=digest(binary),
                          compiled_revision=revision, build=metadata[name])
                for name, binary, revision in (('control', control, CONTROL_REVISION),
                                                ('trial', trial, args.trial_revision))}
    plan = dict(format='aprl-fixed-mpi-followup-plan-v1', workspace=str(workspace), root=str(root),
                prepared_at=datetime.now(timezone.utc).isoformat(),
                source_campaign=str(source), binaries=binaries, runtime_libraries=runtime,
                trial_compiled_revision=trial_revision,
                required_regression=required,
                binding_probe=str(probe), binding_probe_sha256=digest(probe),
                layout=dict(ranks=4, threads=4, placement='compact', cpus=16),
                expected_host='gnode1', physical_core_budget=16, smt_workers=False,
                cases=cases, repeats=REPEATS, order_seed=SEED,
                run_timeout_s=180, total_timeout_s=7200,
                primary_throughput_interval=[20, 80], profile_interval=20,
                comparison='Every control/trial output is bitwise compared to the first completed 4x4 control of that case; no cross-rank comparison.',
                limits=['MPI ranks and threads are fixed at 4x4 for every measured simulation.',
                        'Each case uses byte-identical input cards for both executables.',
                        'Different cases change work and are not implementation speedups.',
                        'The accepted production executable is read-only and never rebuilt here.',
                        'Only separate profiling trials sample NUMA pages; no sampling during throughput.',
                        'NUMA snapshots include runtime pages and are not hardware remote-access counters.',
                        'Communication stages include waiting; nested stage times and rank maxima are not additive.'])
    save(root / 'plan.json', plan)
    frozen = [root / 'plan.json', *sorted(scripts.iterdir()), *sorted(cards.iterdir()),
              control, trial, probe, Path(trial_revision['flags']),
              *[Path(file) for file in required['files']],
              *[Path(file) for file in runtime_files],
              *[Path(value['cache']) for value in metadata.values()]]
    save(root / 'frozen.json', {str(file): digest(file) for file in frozen})
    print(json.dumps(dict(prepared=str(root), layout='4 MPI x 4 threads',
                         pilot_command=f'bash {scripts / "run.sh"} --root {root} --stage pilot',
                         all_command=f'bash {scripts / "run.sh"} --root {root} --stage all'), indent=2))


class Campaign(memory.Campaign):
    def __init__(self, root):
        super().__init__(root)
        self.report['format'] = 'aprl-fixed-mpi-followup-result-v1'
        if self.plan['layout'] != dict(ranks=4, threads=4, placement='compact', cpus=16):
            raise ValueError('The reviewed fixed layout has changed')

    def binding(self):
        super().binding()
        record = self.report['binding'][0]
        if len(record['numa_domains']) != 8:
            raise ValueError('4x4 workers did not cover all eight NUMA domains')
        for rank in range(4):
            workers = [worker for worker in record['workers'] if worker['rank'] == rank]
            if len(workers) != 4 or len({worker['numa'] for worker in workers}) != 2:
                raise ValueError('Each four-thread rank must cover exactly two NUMA domains')

    def profile_summary(self, row, summary):
        stages = [item for item in summary['performance']
                  if item.get('stage') in memory.STAGES and 20 < int(item['step']) <= 80]
        workload = [item for item in summary['performance']
                    if item.get('phase') == 'workload' and 20 < int(item['step']) <= 80]
        if {item['stage'] for item in stages} != set(memory.STAGES) or not workload:
            raise ValueError('Missing profiling stages or particle statistics')
        typed = []
        for stage in memory.STAGES:
            records = [item for item in stages if item['stage'] == stage]
            nsteps = sum(int(item['steps']) for item in records)
            if nsteps != 60:
                raise ValueError('Profile stage does not cover exactly steps 20-80')
            typed.append(dict(stage=stage, steps=nsteps,
                              rank_mean_seconds_per_step=sum(float(item['mean_s']) for item in records) / nsteps,
                              rank_max_seconds_per_step=sum(float(item['max_s']) for item in records) / nsteps,
                              maximum_rank_ids=sorted({int(item['max_rank']) for item in records})))
        load = [dict(step=int(item['step']), particles_total=int(item['particles_total']),
                     particles_min=int(item['particles_min']), particles_max=int(item['particles_max']),
                     particles_mean=int(item['particles_total']) / 4,
                     active_total=int(item['active_total']), frozen_total=int(item['frozen_total']))
                for item in workload]
        if any(item['particles_total'] != self.cases[row['case']]['particles'] for item in load):
            raise ValueError('Profiling particle total changed')
        return dict(stages=typed, workload=load,
                    particle_max_to_mean_max=max(item['particles_max'] / item['particles_mean'] for item in load),
                    any_rank_empty=any(item['particles_min'] == 0 for item in load),
                    interpretation='Communication includes waits and is nested. Rank maxima are not additive. NUMA maps are one process-page snapshot.')

    def one(self, case_name, variant, phase, repeat=0):
        name = f'{phase}-{case_name}-{variant}-r4-t4-rep{repeat}'
        if any(row['name'] == name for row in self.report['runs']):
            raise ValueError('Refusing duplicate run: ' + name)
        case, binary = self.cases[case_name], self.plan['binaries'][variant]
        if digest(binary['path']) != binary['sha256']:
            raise ValueError('Executable changed before launch')
        directory = self.root / name; directory.mkdir()
        shutil.copy2(case['card'], directory / 'run.yaml')
        interval = 20 if phase == 'profile' else 0
        row = dict(name=name, case=case_name, directory=str(directory), label=variant,
                   ranks=4, threads=4, placement='compact', phase=phase, repeat=repeat,
                   cpus=16, profile_interval=interval, valid=False, status='running',
                   binary_sha256=binary['sha256'], compiled_revision=binary['compiled_revision'],
                   configured_grid_cells=case['grid_cells'], configured_particles=case['particles'],
                   loadavg_before=os.getloadavg())
        self.report['runs'].append(row); self.checkpoint(); print('Running', name, flush=True)
        code, wall = self.execute(memory.mpi_command(LAYOUT, binary['path'], [directory / 'run.yaml']),
                                  directory, 4, interval)
        row.update(returncode=code, process_wall_s=wall, loadavg_after=os.getloadavg()); self.checkpoint()
        if code:
            raise ValueError(name + f': process exit {code}; see run.log')
        if digest(directory / 'run.yaml') != case['sha256']:
            raise ValueError('Run card changed')
        manifest = yaml.safe_load((directory / 'run-manifest.yaml').read_text())
        if manifest['source_revision'] != binary['compiled_revision'] or manifest['mpi_size'] != 4:
            raise ValueError('Actual executable revision or MPI size differs')
        if manifest['resolved_mesh']['cells'] != case['cells']:
            raise ValueError('Actual mesh differs from frozen case')
        summary, arrays = self.generic.inspect(self.inspector, row, self.case_info[case_name])
        if summary['reported_threads'] != 4:
            raise ValueError('Executable did not report four effective threads')
        row.update(summary)
        progress = {}
        for line in (directory / 'run.log').read_text(errors='replace').splitlines():
            if not line.startswith('[resource] phase=progress '):
                continue
            fields = dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s,]+)', line))
            step = int(fields['step'])
            if step in progress:
                raise ValueError('Duplicate progress step')
            progress[step] = float(fields['elapsed_s'])
        counts = self.report.setdefault('case_steps', {})
        if case_name not in counts:
            if variant != 'control' or phase != 'pilot' or summary['steps'] < 81:
                raise ValueError('First endpoint reference must be a completed control pilot')
            counts[case_name] = summary['steps']
        if summary['steps'] != counts[case_name] or not {20, 80}.issubset(progress):
            raise ValueError('Physical endpoint/common progress markers differ')
        steady = progress[80] - progress[20]
        if not math.isfinite(steady) or steady <= 0:
            raise ValueError('Invalid steady timing')
        row.update(steady_steps=60, steady_wall_s=steady, steady_seconds_per_step=steady / 60,
                   steady_core_seconds=steady * 16,
                   steady_grid_cell_updates_per_s=case['grid_cells'] * 60 / steady,
                   steady_particle_updates_per_s=case['particles'] * 60 / steady,
                   peak_rss_rank_max_mib=float(summary['final_resource'][-1]['peak_rss_rank_max_mib']),
                   peak_rss_total_mib=float(summary['final_resource'][-1]['peak_rss_total_mib']))
        if phase == 'profile':
            row['profile_summary'] = self.profile_summary(row, summary)
            save(directory / 'profile-summary.json', row['profile_summary'])
        reference = self.report['references'].get(case_name)
        if reference:
            self.compare(row, arrays, reference, True)
        elif variant == 'control' and phase == 'pilot':
            self.report['references'][case_name] = name
        else:
            raise ValueError('No same-MPI control reference')
        row.update(valid=True, status='passed'); self.checkpoint()

    def pilot(self):
        required = self.plan['required_regression']
        actual = required_evidence(Path(required['pass_path']), Path(self.plan['binaries']['trial']['path']))
        if actual != required:
            raise ValueError('Required-regression evidence changed or disappeared')
        if (self.root / 'PILOT_PASS').exists():
            return
        self.binding()
        for case in CASE_NAMES:
            for variant in ('control', 'trial'):
                self.one(case, variant, 'pilot')
        self.report.update(pilot_pass=True, status='pilot_passed'); self.checkpoint()
        (self.root / 'PILOT_PASS').write_text('Fixed 4x4 binding and all three control/trial pilots passed bitwise checks.\n')

    def matrix(self):
        if (self.root / 'PASS').exists():
            return
        if not (self.root / 'PILOT_PASS').exists():
            raise ValueError('Pilot must pass first')
        self.report['status'] = 'running_matrix'; self.checkpoint()
        for case, variant in paired_order(CASE_NAMES, 0):
            self.one(case, variant, 'warmup')
        for repeat in range(1, REPEATS + 1):
            for case, variant in paired_order(CASE_NAMES, repeat):
                self.one(case, variant, 'throughput', repeat)
        self.report['status'] = 'running_profiles'; self.checkpoint()
        for case, variant in paired_order(CASE_NAMES, REPEATS + 1):
            self.one(case, variant, 'profile')
        self.summarize(); self.report.update(status='passed', pass_=True); self.checkpoint()
        (self.root / 'PASS').write_text('All fixed-4x4 control/trial runs and strict bitwise comparisons passed.\n')

    def summarize(self):
        groups = {}
        for row in self.report['runs']:
            if row.get('valid') and row['phase'] == 'throughput':
                groups.setdefault((row['case'], row['label']), []).append(row)
        metrics = ('steady_seconds_per_step', 'steady_core_seconds', 'loop_wall_s',
                   'reported_wall_s', 'process_wall_s', 'peak_rss_rank_max_mib', 'peak_rss_total_mib')
        records = []
        for (case, variant), runs in sorted(groups.items()):
            record = dict(case=case, variant=variant, ranks=4, threads=4, cpus=16, repeats=len(runs))
            for metric in metrics:
                for label, function in (('median', statistics.median), ('min', min), ('max', max)):
                    record[metric + '_' + label] = function(row[metric] for row in runs)
            records.append(record)
        self.report['timing_summary'] = records
        self.write_csv('summary.csv', records)
        comparisons = []
        for case in CASE_NAMES:
            control = {row['repeat']: row for row in groups.get((case, 'control'), [])}
            trial = {row['repeat']: row for row in groups.get((case, 'trial'), [])}
            common = sorted(control.keys() & trial.keys())
            if not common:
                continue
            key = 'steady_seconds_per_step'
            ratios = [control[repeat][key] / trial[repeat][key] for repeat in common]
            c = statistics.median(control[repeat][key] for repeat in common)
            t = statistics.median(trial[repeat][key] for repeat in common)
            comparisons.append(dict(case=case, paired_repeats=len(common),
                                    control_seconds_per_step_median=c, trial_seconds_per_step_median=t,
                                    control_over_trial_medians=c / t, trial_time_reduction_fraction=1 - t / c,
                                    paired_speedup_median=statistics.median(ratios),
                                    paired_speedup_min=min(ratios), paired_speedup_max=max(ratios)))
        self.report['implementation_comparison'] = comparisons
        self.write_csv('comparison.csv', comparisons)
        profiles = []
        for row in self.report['runs']:
            if row.get('valid') and row['phase'] == 'profile':
                for stage in row['profile_summary']['stages']:
                    profiles.append(dict(case=row['case'], variant=row['label'], ranks=4, threads=4,
                                         stage=stage['stage'],
                                         rank_mean_seconds_per_step=stage['rank_mean_seconds_per_step'],
                                         rank_max_seconds_per_step=stage['rank_max_seconds_per_step']))
        self.write_csv('profile-summary.csv', profiles)

    def write_csv(self, filename, rows):
        if rows:
            with (self.root / filename).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)


def idle_host_snapshot():
    """One bounded observation; residual load average alone is not current work."""
    other = []
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            argv = (entry / 'cmdline').read_bytes().split(b'\0')
            name = Path(os.fsdecode(argv[0])).name
            if name in ('aprl', 'mpirun', 'mpiexec', 'orted'):
                other.append(dict(pid=int(entry.name), executable=name))
        except (OSError, ValueError):
            continue
    def counters():
        values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:9]))
        return sum(values), values[3]
    total_before, idle_before = counters()
    time.sleep(2)
    total_after, idle_after = counters()
    delta = total_after - total_before
    if delta <= 0:
        raise ValueError('Could not measure current host idleness')
    fraction = (idle_after - idle_before) / delta
    snapshot = dict(loadavg=os.getloadavg(), sample_seconds=2, cpu_idle_fraction=fraction,
                    observed_compute_processes=other,
                    policy='No other APRL/MPI processes and at least 95% aggregate CPU idle; no waiting or system changes.')
    if other or fraction < .95:
        raise ValueError('gnode1 has active work: ' + json.dumps(snapshot))
    return snapshot


def self_test(workspace):
    """Checks the immutable design and ordering; never launches MPI or APRL."""
    assert LAYOUT == (4, 4, 'compact') and memory.LAYOUTS == [LAYOUT]
    expected = {(case, variant) for case in CASE_NAMES for variant in ('control', 'trial')}
    for repeat in range(REPEATS + 2):
        order = list(paired_order(CASE_NAMES, repeat))
        assert len(order) == 6 and set(order) == expected
        for index in range(0, len(order), 2):
            assert order[index][0] == order[index + 1][0]
            assert order[index][1] != order[index + 1][1]
    command = memory.mpi_command(LAYOUT, 'aprl', ['run.yaml'])
    assert command[command.index('-np') + 1] == '4'
    assert command[command.index('--map-by') + 1] == 'slot:PE=4'
    import tempfile
    with tempfile.TemporaryDirectory(prefix='aprl-fixed-mpi-metadata-') as temporary:
        temporary = Path(temporary)
        trial_build = temporary / 'trial-build'; trial_build.mkdir()
        for filename in ('aprl', 'CMakeCache.txt'):
            shutil.copy2(workspace / CONTROL_BUILD / filename, trial_build / filename)
        flags_relative = Path('CMakeFiles/direct_eb_core.dir/flags.make')
        (trial_build / flags_relative).parent.mkdir(parents=True)
        shutil.copy2(workspace / CONTROL_BUILD / flags_relative, trial_build / flags_relative)
        required_dir = temporary / 'metadata-only-required'; required_dir.mkdir()
        evidence = required_dir / 'PASS'
        evidence.write_text('Metadata-only fixture, not a scientific test result.\n')
        (required_dir / 'EXIT_CODE').write_text('0\n')
        (required_dir / 'build.sha256').write_text(''.join(
            digest(trial_build / name) + '  ' + str(trial_build / name) + '\n'
            for name in ('aprl', 'CMakeCache.txt')))
        shutil.copy2(trial_build / 'CMakeCache.txt', required_dir / 'CMakeCache.txt')
        test_names = ['required.runtime.duration-close', 'required.runtime.duration-barrier-watchdog']
        test_names += ['required.metadata-fixture.' + str(index) for index in range(REQUIRED_TEST_COUNT - 2)]
        (required_dir / 'driver.log').write_text(''.join(
            f'{index:2}/72 Test #{index:2}: {name} .... Passed 0.01 sec\n'
            for index, name in enumerate(test_names, 1)) + '100% tests passed, 0 tests failed out of 72\n')
        detailed = required_dir / 'Testing/Temporary/LastTest.log'
        detailed.parent.mkdir(parents=True)
        detailed.write_text(''.join(
            f'{index}/72 Test: {name}\n' + ('Test Pass Reason:\nMatched required expression\n' if index <= 5 else 'Test Passed.\n')
            for index, name in enumerate(test_names, 1)) + 'End testing: metadata fixture\n')
        args = argparse.Namespace(workspace=workspace, root=temporary / 'prepared',
                                  control=workspace / CONTROL_BUILD / 'aprl',
                                  trial=trial_build / 'aprl', trial_revision=CONTROL_REVISION,
                                  required_pass=evidence,
                                  source_campaign=workspace / 'results/memory_scaling/campaign-20261025')
        prepare(args)
        memory.verify_frozen(args.root)
        campaign = Campaign(args.root)
        assert len(campaign.cases) == 3 and campaign.plan['layout']['ranks'] == 4
        for name in CASE_NAMES:
            original = args.source_campaign / 'cards' / (name + '.yaml')
            copied = args.root / 'cards' / original.name
            assert original.read_bytes() == copied.read_bytes()
        altered = args.root / 'cards/grid016.yaml'
        altered.write_text(altered.read_text() + '\n')
        try:
            memory.verify_frozen(args.root)
        except ValueError:
            pass
        else:
            raise AssertionError('Changed frozen card was accepted')
        args.root = temporary / 'wrong-compiled-revision'
        args.trial_revision = '0' * 40
        try:
            prepare(args)
        except ValueError as error:
            assert 'differs from APRL_SOURCE_REVISION' in str(error)
            assert not args.root.exists()
        else:
            raise AssertionError('Incorrect compiled revision was accepted')
        args.trial_revision = CONTROL_REVISION
        evidence_hashes = required_dir / 'build.sha256'
        good_hashes = evidence_hashes.read_text()
        evidence_hashes.write_text('0' * 64 + '  ' + str(trial_build / 'aprl') + '\n')
        args.root = temporary / 'wrong-required-binary'
        try:
            prepare(args)
        except ValueError as error:
            assert 'SHA256 mismatch' in str(error)
        else:
            raise AssertionError('Required evidence for a different binary was accepted')
        evidence_hashes.write_text(good_hashes)
        exit_file = required_dir / 'EXIT_CODE'
        exit_file.write_text('8\n')
        args.root = temporary / 'failed-required-suite'
        try:
            prepare(args)
        except ValueError as error:
            assert 'EXIT_CODE is not zero' in str(error)
        else:
            raise AssertionError('Failed required suite was accepted')
        exit_file.write_text('0\n')
        saved_cache = required_dir / 'CMakeCache.txt'
        original_cache = saved_cache.read_bytes()
        saved_cache.write_bytes(original_cache + b'\n')
        try:
            required_evidence(evidence, trial_build / 'aprl')
        except ValueError as error:
            assert 'saved build cache differs' in str(error)
        else:
            raise AssertionError('Required evidence with a different cache was accepted')
        saved_cache.write_bytes(original_cache)
        driver = required_dir / 'driver.log'
        original_driver = driver.read_text()
        driver.write_text(original_driver.replace('out of 72', 'out of 70'))
        try:
            required_evidence(evidence, trial_build / 'aprl')
        except ValueError as error:
            assert 'all 72 tests passed' in str(error)
        else:
            raise AssertionError('A partial required suite was accepted')
        driver.write_text(original_driver)
        original_details = detailed.read_text()
        detailed.write_text(original_details.replace('Test Passed.', 'Test Failed.', 1))
        try:
            required_evidence(evidence, trial_build / 'aprl')
        except ValueError as error:
            assert 'detailed test log' in str(error)
        else:
            raise AssertionError('Conflicting required test logs were accepted')
        detailed.write_text(original_details)
        stale_failure = detailed.parent / 'LastTestsFailed.log'
        stale_failure.write_text('1:required.obsolete-failure\n')
        required_evidence(evidence, trial_build / 'aprl')
        cache = trial_build / 'CMakeCache.txt'
        cache.write_text(cache.read_text().replace('CMAKE_CXX_FLAGS:STRING=', 'CMAKE_CXX_FLAGS:STRING=-ffast-math '))
        args.root = temporary / 'different-toolchain'
        try:
            prepare(args)
        except ValueError as error:
            assert 'toolchain differs' in str(error)
        else:
            raise AssertionError('Changed compiler flags were accepted')
    print('Self-test PASS: fixed 4x4 command, paired ordering, exact card copying, frozen-card, compiled revision, required-binary/cache SHA, 72-test evidence/exit-code, stale failure handling, and toolchain mismatch rejection; no compute launched.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--stage', choices=('prepare', 'pilot', 'matrix', 'all'))
    parser.add_argument('--control', type=Path, default=WORKSPACE / CONTROL_BUILD / 'aprl')
    parser.add_argument('--trial', type=Path)
    parser.add_argument('--trial-revision')
    parser.add_argument('--required-pass', type=Path)
    parser.add_argument('--source-campaign', type=Path,
                        default=WORKSPACE / 'results/memory_scaling/campaign-20261025')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test(args.workspace.resolve()); return 0
    if args.root is None or args.stage is None:
        parser.error('--root and --stage are required')
    if args.stage == 'prepare':
        if args.trial is None or args.trial_revision is None or args.required_pass is None:
            parser.error('prepare requires --trial, --trial-revision and --required-pass')
        prepare(args); return 0
    if platform.node() != 'gnode1' or os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') != '1' or os.environ.get('SLURM_JOB_ID'):
        parser.error('Execution requires explicitly authorized gnode1 CPU-only environment without Slurm variables')
    root = args.root.resolve(); memory.verify_frozen(root)
    campaign = Campaign(root)
    if (root / 'PASS').exists() or (args.stage == 'pilot' and (root / 'PILOT_PASS').exists()):
        print('Requested stage already completed; no duplicate launches.', flush=True); return 0
    if campaign.report['status'] == 'failed':
        raise ValueError('Failed campaign cannot be silently resumed; prepare a new directory')
    if not campaign.report['runs']:
        campaign.report['initial_resources'] = idle_host_snapshot()
        campaign.checkpoint()
    try:
        campaign.report['status'] = 'running'; campaign.checkpoint()
        if args.stage in ('pilot', 'all'):
            campaign.pilot()
        if args.stage in ('matrix', 'all'):
            campaign.matrix()
        memory.verify_frozen(root)
        return 0
    except Exception as error:
        campaign.report.update(status='failed', pass_=False, error=str(error))
        campaign.summarize(); campaign.checkpoint()
        (root / 'FAILED').write_text(str(error) + '\n'); print('FAILED:', error, flush=True)
        return 1
    finally:
        campaign.checkpoint()


if __name__ == '__main__':
    raise SystemExit(main())
