#!/usr/bin/env python3
"""Prepare and inspect a bounded, fixed-256-rank production-size followup optimization pair."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import h5py
import numpy as np
import yaml

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
STEPS = 600
EXHAUSTED = 'mesh.duration was exhausted before the configured physical stop'
SPECS = {
    'control': ('cowan-accepted-r2-20261025', 'dc537da20eacc15648cc59360a1f61bdc9ef1176950ac0cdfda21714e2ba58bb',
                '77f1ef95267d037dcd629584e42aedc1b04668f3'),
    'trial': ('compute-followup-20261026', '5edd096eec8d02fef03613a4f74731d90eba5574cb6c109fe26c7fb4f16a9483',
              '8adce330f646701123ef5b312b914b37152af874')}
RUNTIME_LIBRARIES = {
    'libfftw3.so.3': 'b63413c26b5ae054d25df7899e5937a342e3d676254738f2fb891c8260296150',
    'libfftw3_threads.so.3': 'e606bf0c0fec66837c4ab36530701b5bedfd6fcbbb3a03bec7b2e77d4ee228ea',
}
CACHE_KEYS = ('CMAKE_BUILD_TYPE', 'CMAKE_CXX_COMPILER', 'CMAKE_C_COMPILER',
              'CMAKE_CXX_FLAGS', 'CMAKE_CXX_FLAGS_RELEASE', 'CMAKE_C_FLAGS',
              'CMAKE_C_FLAGS_RELEASE', 'MPI_CXX_COMPILER', 'HDF5_C_COMPILER_EXECUTABLE',
              'HDF5_C_LIBRARY_hdf5', 'APRL_REQUIRE_PARALLEL_HDF5', 'APRL_ENABLE_OPENMP',
              'OpenMP_CXX_FLAGS', 'OpenMP_CXX_LIB_NAMES', 'OpenMP_gomp_LIBRARY', 'yaml-cpp_DIR')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def read_yaml(path):
    return yaml.safe_load(path.read_text())


def physics(card, duration=False, input_path=False):
    card = copy.deepcopy(card)
    del card['output']['manifest']
    del card['detectors']['directory']
    del card['energy_ledger']['directory']
    if duration:
        del card['mesh']['duration']
    if input_path:
        del card['beam']['input']['file']
    return card


def redirect(card, case):
    card = copy.deepcopy(card)
    card['output']['manifest'] = str(case / 'run-manifest.yaml')
    card['detectors']['directory'] = str(case / 'detectors')
    card['energy_ledger']['directory'] = str(case / 'ledger')
    return card


def build_metadata(binary):
    cache = binary.parent / 'CMakeCache.txt'
    values = {}
    for line in cache.read_text().splitlines():
        match = re.match(r'([^#/:][^:]*):[^=]+=(.*)$', line)
        if match:
            values[match[1]] = match[2]
    require(all(key in values for key in CACHE_KEYS), 'Incomplete toolchain provenance')
    return {'cache': str(cache), 'cache_sha256': digest(cache),
            'selected': {key: values[key] for key in CACHE_KEYS}}



def compile_flags(binary, revision):
    flags = binary.parent / 'CMakeFiles/direct_eb_core.dir/flags.make'
    contents = flags.read_text()
    pattern = r'-DAPRL_SOURCE_REVISION=\\"([0-9a-f]{40})\\"(?=\s|$)'
    matches = re.findall(pattern, contents)
    require(matches == [revision] and contents.count('-DAPRL_SOURCE_REVISION=') == 1,
            'Actual compile flags do not match expected full source revision: ' + str(binary))
    return {'path': str(flags), 'sha256': digest(flags),
            'normalized': re.sub(pattern, '-DAPRL_SOURCE_REVISION=VERIFIED', contents)}


def required_evidence(workspace, trial):
    required = workspace / 'results/code_optimization_20261026/combined-required'
    names = ('PASS', 'EXIT_CODE', 'build.sha256', 'CMakeCache.txt',
             'driver.log', 'Testing/Temporary/LastTest.log')
    files = [required / name for name in names]
    require(all(path.is_file() for path in files), 'Trial required evidence is missing')
    require((required / 'EXIT_CODE').read_text().strip() == '0' and
            not (required / 'FAIL').exists(), 'Trial required regression did not pass')
    records = (required / 'build.sha256').read_text().splitlines()
    for path in (trial, trial.parent / 'CMakeCache.txt'):
        require(records.count(digest(path) + '  ' + str(path)) == 1,
                'Required evidence does not identify exact executable/build cache')
    require(digest(required / 'CMakeCache.txt') == digest(trial.parent / 'CMakeCache.txt'),
            'Required evidence build cache changed')
    require(re.findall(r'^100% tests passed, 0 tests failed out of ([0-9]+)$',
                       (required / 'driver.log').read_text(), re.M) == ['72'],
            'Expected all 72 required tests to pass')
    acceptance = workspace / 'results/code_optimization_20261026/ACCEPTANCE_REVIEW.json'
    accepted = json.loads(acceptance.read_text())
    require(accepted['status'] == 'PASS' and accepted['binary'] == str(trial) and
            accepted['binary_sha256'] == digest(trial), 'Combined acceptance evidence differs')
    for path, sha in accepted['evidence_sha256'].items():
        require(digest(path) == sha, 'Combined acceptance evidence changed: ' + path)
    files += [acceptance] + [Path(path) for path in accepted['evidence_sha256']]
    return {str(path): digest(path) for path in files}


def prepare(root, workspace):
    reference = workspace / 'results/elegant_retirement_20261024/compact-observed-window/production'
    preflight = reference.parent / 'preflight-r2'
    original = read_yaml(reference / 'run.yaml')
    manifest = read_yaml(reference / 'run-manifest.yaml')
    original_plan = json.loads((reference / 'run-plan.json').read_text())
    require(original['units'] == {'length': 'm', 'time': 's'}, 'SI reference required')
    require(original['mesh']['field_solver'] == 'cowan-z', 'Cowan dt=dz/c required')
    require(manifest['resolved_mesh']['cells'] == [160, 160, 46592], 'Unexpected production geometry')
    require(original['runtime']['resource_monitor']['progress_interval_steps'] == 50,
            'Preserve production progress cadence')
    require((preflight / 'PASS').is_file(), 'Production geometry/functionality preflight missing')
    variants = {}
    for name, (build, sha, revision) in SPECS.items():
        binary = workspace / 'build' / build / 'aprl'
        require(digest(binary) == sha, name + ' binary differs from verified implementation')
        variants[name] = {'binary': str(binary), 'binary_sha256': sha,
                          'source_revision': revision, 'toolchain': build_metadata(binary),
                          'compile_flags': compile_flags(binary, revision)}
    require(variants['control']['toolchain']['selected'] == variants['trial']['toolchain']['selected'],
            'Control/trial compiler, flags or dependencies differ')
    require(variants['control']['compile_flags']['normalized'] ==
            variants['trial']['compile_flags']['normalized'],
            'Actual compile flags differ beyond source revision')
    required = required_evidence(workspace, Path(variants['trial']['binary']))
    runtime_dir = workspace / 'build/cowan-accepted-r2-20261025/runtime-libs'
    runtime = {str(runtime_dir / name): sha for name, sha in RUNTIME_LIBRARIES.items()}
    for path, sha in runtime.items():
        require(digest(path) == sha, 'Accepted runtime library differs: ' + path)
    require(not root.exists(), 'Use a new result directory')
    root.mkdir(parents=True)
    shutil.copytree(Path(__file__).resolve().parent, root / 'scripts',
                    ignore=shutil.ignore_patterns('__pycache__'))
    refs = root / 'reference'
    refs.mkdir()
    for name in ('run.yaml', 'run-manifest.yaml', 'run-plan.json'):
        shutil.copy2(reference / name, refs / name)
    for name in ('run.yaml', 'run-manifest.yaml', 'particles.h5'):
        shutil.copy2(preflight / 'functional' / name, refs / ('preflight-' + name))
    small = read_yaml(refs / 'preflight-run.yaml')
    dt = original['mesh']['cell_size'][2] / 299792458.0
    gamma = original['mesh']['boost_gamma']
    beta = math.sqrt(1 - 1 / gamma**2)
    longitudinal = manifest['resolved_mesh']['longitudinal']
    detector_geometry = []
    for detector in original['detectors']['field_planes']:
        first = longitudinal['translation_m'] + (
            detector['z'] - original['beam']['reference']['activation_plane_z']) / gamma
        last = first - beta * 299792458.0 * STEPS * dt
        lo, hi = longitudinal['physical_bounds_m']
        require(lo < last <= first < hi, 'Field plane leaves preserved physical window during timing sample')
        require(math.isclose(detector['rhythm'], dt / gamma, rel_tol=1e-12),
                'Field cadence must remain one sample per field step for this checker')
        detector_geometry.append({'name': detector['name'], 'box_z_at_step_0_m': first,
                                  'box_z_at_step_600_m': last, 'expected_samples': STEPS + 1})
    for name in SPECS:
        case = root / name
        case.mkdir()
        card = redirect(original, case)
        card['mesh']['duration'] = (STEPS - 0.5) * dt
        require(physics(card, duration=True) == physics(original, duration=True),
                'Unexpected production card change')
        (case / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
        mini = root / 'preflight' / name
        mini.mkdir(parents=True)
        small_card = redirect(small, mini)
        small_card['beam']['input']['file'] = str(refs / 'preflight-particles.h5')
        require(physics(small_card, input_path=True) == physics(small, input_path=True),
                'Unexpected preflight card change')
        (mini / 'run.yaml').write_text(yaml.safe_dump(small_card, sort_keys=False))
    plan = {'format': 'aprl-full-scale-followup-pair-v1', 'prepared_at': datetime.now(timezone.utc).isoformat(),
            'workspace': str(workspace), 'reference': str(reference), 'variants': variants,
            'source_revision': subprocess.check_output(['git', '-C', str(workspace / 'source/aprl'),
                                                       'rev-parse', 'HEAD'], text=True).strip(),
            'steps': STEPS, 'dt_box_s': dt, 'expected_mesh': manifest['resolved_mesh'],
            'field_sample_geometry': detector_geometry,
            'nodes': 2, 'ranks': 256, 'threads_per_rank': 1, 'exclusive': True,
            'run_timeout_s': 3000, 'preflight_timeout_s': 600, 'total_limit_s': 7200,
            'runtime_libraries': {'directory': str(runtime_dir), 'files': runtime},
            'intervals': [[50, 100], [100, 400], [400, 600]],
            'changes_from_production': ['mesh.duration: bounded timing stop', 'output paths'],
            'expected_termination': '600 complete steps, duration exhausted, closed incomplete diagnostics',
            'limitations': ['One paired early-production sample, not repeated throughput statistics.',
                            'Control already includes accepted Cowan optimization; only incremental followup gains are measured.',
                            'Not a complete production run or validation of late retirement/radiation physics.',
                            'Same MPI layout only; no claim about existing cross-rank differences.',
                            'No subsequent job submission or automatic production replacement.']}
    for key in ('input_file', 'input_sha256', 'input_size_bytes', 'input_mtime_ns'):
        plan[key] = original_plan[key]
    input_stat = Path(plan['input_file']).stat()
    require((input_stat.st_size, input_stat.st_mtime_ns) ==
            (plan['input_size_bytes'], plan['input_mtime_ns']), 'Production input identity changed')
    plan['required_evidence'] = required
    save(root / 'plan.json', plan)
    frozen = sorted(p for p in root.rglob('*') if p.is_file())
    frozen += [Path(v['binary']) for v in variants.values()]
    frozen += [Path(v['toolchain']['cache']) for v in variants.values()]
    frozen += [Path(v['compile_flags']['path']) for v in variants.values()]
    frozen += [Path(path) for path in runtime]
    frozen += [workspace / 'build/cluster-tools/binding_probe']
    frozen += [Path(p) for p in plan['required_evidence']]
    (root / 'frozen.sha256').write_text(''.join(f'{digest(p)}  {p}\n' for p in frozen))
    print(root)


def verify(root):
    for line in (root / 'frozen.sha256').read_text().splitlines():
        sha, file = line.split('  ', 1)
        require(digest(file) == sha, 'Frozen file changed: ' + file)
    plan = json.loads((root / 'plan.json').read_text())
    file = Path(plan['input_file'])
    stat = file.stat()
    require((stat.st_size, stat.st_mtime_ns) == (plan['input_size_bytes'], plan['input_mtime_ns']),
            'Production input stat changed')
    require(digest(file) == plan['input_sha256'], 'Production input SHA256 changed')
    save(root / 'input-verified.json', {'pass': True, 'input_sha256': plan['input_sha256'],
                                     'checked_at': datetime.now(timezone.utc).isoformat()})


def finite(array, label):
    if array.dtype.names:
        for name in array.dtype.names:
            finite(array[name], label + '.' + name)
    elif array.dtype.kind in 'fc':
        require(np.isfinite(array).all(), 'Nonfinite scientific values: ' + label)


def field_cadence_report(times, dt_box_s, boost_gamma, expected_samples,
                         lab_origin_time_s=0.0):
    """Check the whole time grid with a bound derived from FP64 rounding.

    The detector writes lab_origin + t_box/gamma + a constant flight time.
    At a ~19 ns absolute time, adjacent 4.5 as samples cannot have an exactly
    constant FP64 difference. Anchor every sample to the first timestamp so
    a cumulative cadence error cannot hide inside independent pairwise tests.
    """
    times = np.asarray(times)
    require(times.dtype == np.dtype('float64') and times.ndim == 1 and
            len(times) == expected_samples and expected_samples >= 2,
            'Unexpected field time-grid shape or precision')
    require(np.isfinite(times).all() and np.all(np.diff(times) > 0),
            'Nonfinite or nonmonotonic field time grid')
    require(math.isfinite(dt_box_s) and dt_box_s > 0 and
            math.isfinite(boost_gamma) and boost_gamma >= 1 and
            math.isfinite(lab_origin_time_s), 'Invalid cadence reference')
    wide = np.longdouble
    k = np.arange(expected_samples, dtype=wide)
    unit_roundoff = wide(np.finfo(np.float64).eps) / 2
    require(k[-1] * unit_roundoff < 1, 'Time-grid rounding bound is undefined')
    dt, gamma = wide(dt_box_s), wide(boost_gamma)
    exact_offsets = k * dt / gamma
    measured_offsets = times.astype(wide) - wide(times[0])
    errors = np.abs(measured_offsets - exact_offsets)
    # Each persisted endpoint has at most half an ULP of rounding error.
    ulps = np.abs(np.spacing(times)).astype(wide)
    endpoint_bound = (ulps + ulps[0]) / 2
    # t_box += dt is performed once per step. gamma_n bounds the accumulated
    # FP64 addition error; division and the small lab-origin addition follow
    # boostframe.cpp::labTFromLabZBoxT. The constant flight term cancels in
    # anchored differences. The factor two also covers the printed origin's
    # insignificant decimal rounding when supplied from the simulation log.
    accumulation_bound = (k * unit_roundoff / (1 - k * unit_roundoff)) * k * dt / gamma
    conversion_bound = 2 * unit_roundoff * (
        2 * abs(wide(lab_origin_time_s)) + np.abs(exact_offsets) + accumulation_bound)
    bound = endpoint_bound + accumulation_bound + conversion_bound
    require(np.all(errors <= bound), 'Field sampling cadence changed beyond FP64 time-grid rounding')
    return {'method': 'anchored sample grid with FP64 endpoint ULP and accumulated-time bounds',
            'samples': int(expected_samples), 'expected_cadence_s': float(dt / gamma),
            'first_time_s': float(times[0]), 'last_time_s': float(times[-1]),
            'maximum_absolute_anchor_error_s': float(np.max(errors)),
            'maximum_allowed_anchor_error_s': float(np.max(bound)),
            'maximum_error_fraction_of_bound': float(np.max(errors / bound)),
            'timestamp_ulp_min_s': float(np.min(ulps)),
            'timestamp_ulp_max_s': float(np.max(ulps)),
            'maximum_accumulated_box_time_bound_lab_s': float(np.max(accumulation_bound)),
            'maximum_conversion_roundoff_bound_s': float(np.max(conversion_bound))}


def slices(dataset):
    if not dataset.shape:
        yield ()
    else:
        # Bound field slabs to about 16 MiB, particle slabs to 65536 records.
        row_bytes = max(1, dataset.dtype.itemsize * math.prod(dataset.shape[1:]))
        count = max(1, min(65536, 16 * 1024 * 1024 // row_bytes))
        for start in range(0, dataset.shape[0], count):
            yield slice(start, min(start + count, dataset.shape[0]))


def dataset_names(file):
    names = []
    file.visititems(lambda name, obj: names.append(name) if isinstance(obj, h5py.Dataset) else None)
    return names


def output_files(case):
    card = read_yaml(case / 'run.yaml')
    names = [Path('detectors') / (d['name'] + '.h5')
             for key in ('field_planes', 'particle_planes') for d in card['detectors'][key]]
    names.append(Path('ledger') / card['energy_ledger']['filename'])
    require(set((case / 'detectors').glob('*.h5')) == {case / n for n in names if n.parts[0] == 'detectors'},
            'Missing or unexpected detector output')
    return names


def inspect_outputs(case, completed):
    result = []
    for relative in output_files(case):
        with h5py.File(case / relative, 'r') as file:
            group_name = next(n for n in ('field_plane', 'particle_plane', 'energy_ledger') if n in file)
            g = file[group_name]
            field = group_name == 'field_plane'
            count = int(g['committed_samples' if field else 'committed_records'][()])
            data = g['fields' if field else 'records']
            require(int(g['complete'][()]) == int(completed) and len(data) == count,
                    'Output transaction/complete marker mismatch: ' + str(relative))
            for name in dataset_names(file):
                for selection in slices(file[name]):
                    finite(np.asarray(file[name][selection]), str(relative) + '/' + name)
            if field:
                times = g['time_s'][:]
                require(len(times) == count and np.all(np.diff(times) > 0), 'Invalid field time axis')
            result.append({'file': str(relative), 'committed': count, 'complete': int(completed)})
    return result


def bitwise(left, right, label):
    require(left.dtype == right.dtype and left.shape == right.shape, 'Array layout differs: ' + label)
    if left.dtype.names:
        # Compound padding bytes are not scientific values; compare every named field.
        for name in left.dtype.names:
            bitwise(left[name], right[name], label + '.' + name)
    else:
        require(left.tobytes() == right.tobytes(), 'Bitwise scientific difference: ' + label)


def compare_outputs(left, right):
    files = output_files(left)
    require(files == output_files(right), 'Output file sets differ')
    compared = []
    for relative in files:
        with h5py.File(left / relative, 'r') as a, h5py.File(right / relative, 'r') as b:
            names = dataset_names(a)
            require(names == dataset_names(b), 'Dataset sets differ: ' + str(relative))
            object_names = ['']
            a.visit(object_names.append)
            other_names = ['']
            b.visit(other_names.append)
            require(object_names == other_names, 'HDF5 object sets differ: ' + str(relative))
            ignored = {'run_id', 'run_manifest', 'configuration_digest_fnv1a64', 'source_revision'}
            for object_name in object_names:
                ao, bo = (a, b) if not object_name else (a[object_name], b[object_name])
                keys = set(ao.attrs) - ignored
                require(keys == set(bo.attrs) - ignored, 'Scientific attribute sets differ')
                for key in keys:
                    bitwise(np.asarray(ao.attrs[key]), np.asarray(bo.attrs[key]),
                            str(relative) + '/' + object_name + '@' + key)
            for name in names:
                require(a[name].shape == b[name].shape and a[name].dtype == b[name].dtype,
                        'Dataset layout differs: ' + name)
                for selection in slices(a[name]):
                    bitwise(np.asarray(a[name][selection]), np.asarray(b[name][selection]), str(relative) + '/' + name)
                compared.append(str(relative) + '/' + name)
    return {'pass': True, 'mode': 'same-256-rank bitwise scientific datasets, attributes and complete markers',
            'ignored_identity_attributes': sorted(ignored), 'datasets': compared}


def check_manifest(root, case, variant, small=False):
    plan = json.loads((root / 'plan.json').read_text())
    manifest = read_yaml(case / 'run-manifest.yaml')
    reference = read_yaml(root / 'reference' / ('preflight-run-manifest.yaml' if small else 'run-manifest.yaml'))
    require(manifest['source_revision'] == plan['variants'][variant]['source_revision'], 'Wrong executable revision')
    require(manifest['mpi_size'] == 256, 'Wrong MPI layout')
    for key in ('resolved_mesh', 'resolved_beam_reference', 'resolved_initial_self_field', 'resolved_runtime',
                'build', 'mpi_library'):
        require(manifest[key] == reference[key], 'Manifest differs from preserved reference: ' + key)
    card = read_yaml(case / 'run.yaml')
    expected = read_yaml(root / 'reference' / ('preflight-run.yaml' if small else 'run.yaml'))
    require(physics(card, duration=not small, input_path=small) ==
            physics(expected, duration=not small, input_path=small), 'Unexpected card changes')
    require(yaml.safe_load(manifest['configuration_yaml']) == card, 'Manifest does not describe frozen card')
    require(digest(plan['variants'][variant]['binary']) == plan['variants'][variant]['binary_sha256'],
            'Executable changed during run')
    return manifest


def check_preflight(root):
    reports = {}
    for variant in SPECS:
        case = root / 'preflight' / variant
        require((case / 'EXIT_CODE').read_text().strip() == '0', 'Small run failed: ' + variant)
        check_manifest(root, case, variant, small=True)
        log = (case / 'aprl.log').read_text(errors='replace')
        require('Configured stop reached:' in log and 'Simulation failure' not in log, 'Small run did not finish')
        require(re.search(r'Compute threads per MPI rank:\s*1;', log),
                'Small run did not use exactly one compute thread per rank')
        activation = re.search(r'Frozen lab-plane injection: activations=(\d+).*?still frozen=(\d+)', log)
        retirement = re.search(r'Experimental particle retirement: entries=(\d+).*?completed=(\d+), residual Q at taper exit=([^ ]+)', log)
        require(activation and tuple(map(int, activation.groups())) == (512, 0), 'Incomplete small held activation')
        require(retirement and tuple(map(int, retirement.group(1, 2))) == (512, 512), 'Incomplete small retirement')
        residual = float(retirement.group(3))
        require(math.isfinite(residual) and abs(residual) <= 64 * np.finfo(float).eps * 512 * 1.602176634e-19,
                'Retirement residual charge exceeds unchanged roundoff bound')
        reports[variant] = {'outputs': inspect_outputs(case, True), 'residual_charge_C': residual}
        with h5py.File(root / 'reference/preflight-particles.h5') as f:
            source = f['particles/records'][:]
        for name in ('activation', 'retirement-entry'):
            with h5py.File(case / 'detectors' / (name + '.h5')) as f:
                records = f['particle_plane/records'][:]
            records = records[np.argsort(records['source_id'])]
            require(len(records) == 512 and np.array_equal(records['source_id'], source['source_id']),
                    'Small detector cohort differs')
            if name == 'activation':
                require(np.max(np.abs(records['position_m'][:, :2] - source['plane_position_m'])) <= 1e-12,
                        'Held transverse activation positions differ')
        with h5py.File(case / 'ledger/energy-ledger.h5') as f:
            records = f['energy_ledger/records'][:]
            empty = np.flatnonzero(records['active_macroparticles'] == 0)
            require(len(empty) > 0 and len(records) - empty[0] >= 20 and
                    int(records[-1]['removed_macroparticles']) == 512, 'Missing post-retirement propagation')
    reports['comparison'] = compare_outputs(root / 'preflight/control', root / 'preflight/trial')
    save(root / 'preflight-result.json', {'pass': True, 'reports': reports})
    (root / 'PREFLIGHT_PASS').write_text('Both 256-rank small held/retirement runs complete and bitwise equal.\n')


def check_duration_exit(log, code, steps=STEPS):
    """Accept only the expected bounded stop after root closed all output.

    MPI_Abort may truncate a peer's error prefix after the readiness barrier.
    Keep incomplete prefixes as evidence, but never discard a fully formed
    rank error or let a following expected error conceal its message.
    """
    require(code == 1 and EXHAUSTED in log,
            'Expected clean duration-exhaustion exit code 1')
    expected_ranks = []
    truncated = []
    for position in re.finditer(r'Simulation failure', log):
        tail = log[position.start():]
        segment = re.split(r'\n|(?=Simulation failure)|(?=srun:)', tail[len('Simulation failure'):],
                           maxsplit=1)[0]
        segment = 'Simulation failure' + segment
        full = re.match(r'^Simulation failure on MPI rank ([0-9]+):[ \t]*', segment)
        if full:
            message = segment[full.end():]
            require(message.startswith(EXHAUSTED),
                    'Unexpected complete rank failure: ' + segment)
            expected_ranks.append(int(full.group(1)))
        else:
            require(re.fullmatch(r'Simulation failure(?: on MPI rank(?: [0-9]*)?)?[ \t]*', segment),
                    'Unexpected malformed simulation failure: ' + segment)
            truncated.append(segment)
    require(expected_ranks, 'Missing a complete expected duration failure message')
    for phase in ('detector-close-complete', 'ledger-close-complete',
                  'final-resource-report-complete', 'duration-exhaustion-readiness-sync-enter'):
        require(re.search(r'\[teardown\] rank=0\b[^\n]* phase=' + phase + r'\b', log),
                'Missing root output-close/readiness marker: ' + phase)
    require(re.search(r'\[resource\] phase=final step=' + str(steps) + r'\b', log),
            'Missing final complete-step report')
    return {'expected_complete_rank_messages': len(expected_ranks),
            'truncated_rank_prefixes': truncated,
            'root_outputs_closed_before_readiness_barrier': True,
            'readiness_complete_marker_observed': bool(re.search(
                r'\[teardown\] rank=0\b[^\n]* phase=duration-exhaustion-readiness-sync-complete\b', log))}


def check_run(root, variant):
    plan = json.loads((root / 'plan.json').read_text())
    case = root / variant
    require((root / 'PREFLIGHT_PASS').is_file(), 'Paired small preflight has not passed')
    code = int((case / 'EXIT_CODE').read_text())
    log = (case / 'aprl.log').read_text(errors='replace')
    exit_report = check_duration_exit(log, code)
    require(re.search(r'Compute threads per MPI rank:\s*1;', log),
            'Full run did not use exactly one compute thread per rank')
    points = {int(step): float(elapsed) for step, elapsed in re.findall(
        r'\[resource\] phase=progress step=(\d+) elapsed_s=([\d.eE+-]+)', log)}
    require(sorted(points) == list(range(50, 601, 50)) and
            all(math.isfinite(v) and v > 0 for v in points.values()) and
            all(points[b] > points[a] for a, b in zip(sorted(points), sorted(points)[1:])),
            'Incomplete or nonmonotonic timing sequence')
    check_manifest(root, case, variant)
    dt = re.search(r'\[resource\] phase=pre-run .*?dt_s=([\d.eE+-]+)', log)
    require(dt and math.isclose(float(dt.group(1)), plan['dt_box_s'], rel_tol=1e-9, abs_tol=0),
            'Logged field timestep differs')
    outputs = inspect_outputs(case, False)
    card = read_yaml(case / 'run.yaml')
    origin = re.search(r'origin_lab_t=([\d.eE+-]+)', log)
    require(origin is not None, 'Missing Lorentz lab time origin')
    cadence_reports = []
    for detector in card['detectors']['field_planes']:
        with h5py.File(case / 'detectors' / (detector['name'] + '.h5')) as f:
            g = f['field_plane']
            require(g['fields'].shape == (601, 160, 160), 'Unexpected field sample shape')
            cadence_reports.append({'detector': detector['name'], **field_cadence_report(
                g['time_s'][:], plan['dt_box_s'], card['mesh']['boost_gamma'], 601,
                float(origin.group(1)))})
    with h5py.File(case / 'ledger/energy-ledger.h5') as f:
        records = f['energy_ledger/records'][:]
        require(np.array_equal(records['step'], np.arange(0, 601, 10)), 'Ledger cadence changed')
        require(math.isclose(float(records[-1]['time_box_s']), 600 * plan['dt_box_s'], rel_tol=1e-12),
                'Ledger final physical time differs')
    intervals = [{'start_step': a, 'end_step': b, 'seconds_per_step': (points[b] - points[a]) / (b - a)}
                 for a, b in plan['intervals']]
    report = {'pass': True, 'variant': variant, 'simulation_complete': False,
              'expected_duration_exhaustion': True, 'steps': STEPS, 'exit_code': code,
              'duration_exit_validation': exit_report,
              'timing_intervals': intervals, 'progress_elapsed_s': points, 'outputs': outputs,
              'field_cadence_validation': cadence_reports,
              'limitations': plan['limitations']}
    save(case / 'timing-result.json', report)
    (case / 'TIMING_PASS').write_text('600-step bounded timing run checked; not a complete physical simulation.\n')
    print(json.dumps(intervals, indent=2))


def compare(root):
    for variant in SPECS:
        require((root / variant / 'TIMING_PASS').is_file(), 'Run not checked: ' + variant)
    comparison = compare_outputs(root / 'control', root / 'trial')
    reports = {v: json.loads((root / v / 'timing-result.json').read_text()) for v in SPECS}
    intervals = []
    for a, b in zip(reports['control']['timing_intervals'], reports['trial']['timing_intervals']):
        require((a['start_step'], a['end_step']) == (b['start_step'], b['end_step']), 'Timing windows differ')
        old, new = a['seconds_per_step'], b['seconds_per_step']
        require(old > 0 and new > 0, 'Nonpositive interval time')
        intervals.append({'start_step': a['start_step'], 'end_step': a['end_step'],
                          'control_s_per_step': old, 'trial_s_per_step': new,
                          'speedup': old / new, 'time_reduction_fraction': 1 - new / old})
    plan = json.loads((root / 'plan.json').read_text())
    trial_nodes = (root / 'allocation-nodes.txt').read_text().splitlines()
    require(len(trial_nodes) == plan['nodes'] and len(set(trial_nodes)) == plan['nodes'],
            'Missing or invalid allocation node provenance')
    allocation_comparison = {'control_nodes': trial_nodes, 'trial_nodes': trial_nodes,
                             'same_node_set': True, 'separate_allocations': False}
    if 'reused_control' in plan:
        trial_nodes = (root / 'allocation-nodes.txt').read_text().splitlines()
        require(len(trial_nodes) == plan['nodes'] and len(set(trial_nodes)) == plan['nodes'],
                'Missing or invalid trial allocation node provenance')
        control_nodes = plan['reused_control']['nodes']
        allocation_comparison = {'control_nodes': control_nodes, 'trial_nodes': trial_nodes,
                                 'same_node_set': set(control_nodes) == set(trial_nodes),
                                 'separate_allocations': True}
    # Recheck all frozen identities and the shared input after both runs and
    # the scientific comparison; changed inputs must not receive PASS.
    verify(root)
    save(root / 'comparison.json', {'pass': True, 'scientific_comparison': comparison,
                                    'timing_intervals': intervals, 'limitations': plan['limitations'],
                                    'allocation_comparison': allocation_comparison})
    (root / 'PASS').write_text('Bounded same-256-rank comparison complete; partial production timing only.\n')
    print(json.dumps(intervals, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'verify', 'check-preflight', 'check-run', 'compare'))
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--workspace', default=WORKSPACE, type=Path)
    parser.add_argument('--variant', choices=tuple(SPECS))
    args = parser.parse_args()
    root = args.root.resolve()
    if args.mode == 'prepare':
        prepare(root, args.workspace.resolve())
    elif args.mode == 'verify':
        verify(root)
    elif args.mode == 'check-preflight':
        check_preflight(root)
    elif args.mode == 'check-run':
        if args.variant is None:
            parser.error('check-run requires --variant')
        check_run(root, args.variant)
    else:
        compare(root)


if __name__ == '__main__':
    main()
