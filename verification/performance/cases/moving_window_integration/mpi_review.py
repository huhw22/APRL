#!/usr/bin/env python3
"""Read-only cross-layout review of frozen moving-window integration outputs.

Uses one CPU, streams every field sample and trajectory record, preserves E/B
units, and computes per-ID kinetic-energy differences without subtracting two
rounded total energies. No simulation or external job is started.
"""
import os
if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'

import argparse
import hashlib
import json
from pathlib import Path
import re
import struct

import h5py
import numpy as np
import yaml

LD = np.longdouble
REST_EV = LD('510998.9499961642')
IDENTITY = ('source_id', 'particle_id', 'charge_C', 'mass_kg', 'weight')


def require(test, message):
    if not test:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def finite(data):
    if data.dtype.names:
        return all(finite(data[name]) for name in data.dtype.names)
    return bool(np.isfinite(data).all()) if data.dtype.kind in 'fc' else True


def array_error(reference, candidate):
    require(reference.shape == candidate.shape, 'Compared array shapes differ')
    first = np.asarray(reference, dtype=LD)
    difference = np.asarray(candidate, dtype=LD) - first
    norm, delta = np.sum(first * first), np.sum(difference * difference)
    return {'bitwise_equal': reference.dtype == candidate.dtype and reference.tobytes() == candidate.tobytes(),
            'max_abs': float(np.max(abs(difference))) if difference.size else 0.,
            'rms_abs': float(np.sqrt(delta / difference.size)) if difference.size else 0.,
            'relative_l2': float(np.sqrt(delta / norm)) if norm else None}


def delta_energy(reference_u, candidate_u):
    a, b = np.asarray(reference_u, dtype=LD), np.asarray(candidate_u, dtype=LD)
    return np.sum((b-a)*(b+a), axis=1) / (
        np.sqrt(1+np.sum(a*a, axis=1)) + np.sqrt(1+np.sum(b*b, axis=1))) * REST_EV


def weighted_energy(reference, candidate):
    weights = abs(np.asarray(reference['charge_C'], dtype=LD))
    weights /= weights.sum()
    delta = delta_energy(reference['proper_velocity'], candidate['proper_velocity'])
    return {'mean_eV': float(np.sum(weights*delta)),
            'rms_eV': float(np.sqrt(np.sum(weights*delta*delta))),
            'maximum_abs_eV': float(np.max(abs(delta)))}


def read_plane(root, name):
    with h5py.File(root/'moving/detectors'/f'{name}.h5') as stream:
        records = stream['particle_plane/records'][:]
    return records[np.argsort(records['source_id'])]


def execution_evidence(root, ranks, threads):
    text = (root/'binding.log').read_text()
    rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)', text)
    require(len(rows) == ranks*threads and all(row[-1] == 'PASS' for row in rows), 'Binding records incomplete')
    hosts = {}
    for rank in range(ranks):
        owned = [row for row in rows if int(row[0]) == rank]
        require({int(row[1]) for row in owned} == set(range(threads)), 'Binding thread IDs differ')
        require(len({row[4:6] for row in owned}) == threads, 'Bound threads overlap physical cores')
        require(len({row[2] for row in owned}) == 1, 'One rank appears on multiple hosts')
        hosts[str(rank)] = owned[0][2]
    if (root/'JOB_EXIT_CODE').exists():
        require((root/'JOB_EXIT_CODE').read_text().strip() == '0', 'Nonzero persisted case exit code')
    return {'actual_rank_hosts': hosts, 'threads_per_rank': threads, 'binding_records': len(rows),
            'binding_sha256': sha(root/'binding.log'),
            'persisted_case_exit_zero': (root/'JOB_EXIT_CODE').exists(),
            'scope': 'Case completion/launcher and actual binding checked. Final scheduler accounting is independently checked by the task owner.'}


def inspect_trajectory(path, release):
    """Bounded-memory scan: no full multi-million-record trajectory sort."""
    count = len(release)
    with h5py.File(path) as stream:
        group = stream['trajectory']
        records = group['records']
        first = np.zeros(count, dtype=records.dtype)
        terminal = np.zeros(count, dtype=records.dtype)
        first_time, last_time = np.full(count, np.inf), np.full(count, -np.inf)
        periodic_counts, terminal_counts = np.zeros(count, dtype=np.int64), np.zeros(count, dtype=np.int64)
        for start in range(0, len(records), 65536):
            chunk = records[start:start+65536]
            require(finite(chunk), 'Nonfinite trajectory records: '+str(path))
            require(np.all((chunk['source_id'] >= 1) & (chunk['source_id'] <= count)), 'Unexpected trajectory ID')
            index = chunk['source_id'].astype(np.int64)-1
            for key in ('particle_id', 'charge_C', 'weight'):
                require(np.array_equal(chunk[key], release[key][index]), 'Trajectory identity/weight changed: '+key)
            regular = chunk[chunk['event_type'] == 0]
            if len(regular):
                regular = regular[np.lexsort((regular['time_s'], regular['source_id']))]
                rid = regular['source_id'].astype(np.int64)-1
                same = rid[1:] == rid[:-1]
                require(np.all(np.diff(regular['time_s'])[same] > 0), 'Duplicate periodic ID/time')
                unique, starts, lengths = np.unique(rid, return_index=True, return_counts=True)
                ends = starts+lengths-1
                require(np.all(regular['time_s'][starts] > last_time[unique]), 'Per-rank trajectory time not increasing')
                new = periodic_counts[unique] == 0
                first[unique[new]] = regular[starts[new]]
                first_time[unique[new]] = regular['time_s'][starts[new]]
                last_time[unique] = regular['time_s'][ends]
                periodic_counts[unique] += lengths
            events = chunk[chunk['event_type'] != 0]
            if len(events):
                require(np.all(events['event_type'] == 2) and np.all(events['boundary_face'] == 5),
                        'Unexpected terminal event type/face')
                eid = events['source_id'].astype(np.int64)-1
                require(len(np.unique(eid)) == len(eid) and np.all(terminal_counts[eid] == 0), 'Duplicate terminal event')
                terminal[eid] = events
                terminal_counts[eid] += 1
        summary = {'records': len(records), 'rank': int(group.attrs['mpi_rank']),
                   'mpi_size': int(group.attrs['mpi_size']),
                   'periodic_records': int(periodic_counts.sum()),
                   'terminal_records': int(terminal_counts.sum()),
                   'source_ids_seen': int(np.count_nonzero(periodic_counts+terminal_counts))}
    return summary, first, terminal, first_time, last_time, periodic_counts, terminal_counts


def inspect_case(root):
    plan = json.loads((root/'plan.json').read_text())
    count, ranks = plan['source']['records'], plan['ranks']
    planes = {name: read_plane(root, name) for name in ('release', 'aux', 'exit')}
    ids = np.arange(1, count+1, dtype=np.uint64)
    for name, records in planes.items():
        require(np.array_equal(records['source_id'], ids), 'Missing/duplicate detector IDs: '+name)
        require(np.all(records['crossing_direction'] == 1), 'Backward detector crossing')
        for key in IDENTITY:
            require(np.array_equal(records[key], planes['release'][key]), 'Plane identity/weight changed: '+key)
    files, trajectories = {}, []
    for path in sorted((root/'moving').glob('*/*.h5')):
        with h5py.File(path) as stream:
            require(len(stream) == 1, 'Unexpected H5 groups: '+str(path))
            group = stream[next(iter(stream))]
            require(int(group['complete'][()]) == 1, 'Incomplete H5: '+str(path))
            item = {'complete': True, 'datasets': {}}
            for name, dataset in group.items():
                if not isinstance(dataset, h5py.Dataset):
                    continue
                item['datasets'][name] = {'shape': list(dataset.shape), 'dtype': str(dataset.dtype)}
                if path.parent.name == 'trajectory' and name == 'records':
                    continue  # Streamed and validated below, including identity.
                pieces = [dataset[()]] if not dataset.shape else (
                    dataset[k:k+8] for k in range(0, len(dataset), 8))
                require(all(finite(np.asarray(piece)) for piece in pieces), 'Nonfinite H5: '+str(path)+'/'+name)
            for committed, values in (('committed_records', 'records'), ('committed_samples', 'time_s')):
                if committed in group:
                    require(int(group[committed][()]) == len(group[values]), 'Truncated H5: '+str(path))
            if 'fields' in group:
                times = group['time_s'][:]
                require(len(group['fields']) == len(times) and np.all(np.diff(times) > 0), 'Invalid field time axis')
        if path.parent.name == 'trajectory':
            trajectory = inspect_trajectory(path, planes['release'])
            require(trajectory[0]['mpi_size'] == ranks, 'Trajectory MPI metadata differs')
            trajectories.append(trajectory)
            item.update(trajectory[0])
        files[str(path.relative_to(root))] = item
    require(len(files) == 5+ranks and len(trajectories) == ranks, 'Unexpected output H5 set')
    require(sorted(t[0]['rank'] for t in trajectories) == list(range(ranks)), 'Missing/duplicate trajectory rank')
    first = trajectories[0][1].copy()
    first_times = trajectories[0][3].copy()
    terminal = trajectories[0][2].copy()
    terminal_counts = trajectories[0][6].copy()
    periodic_counts = trajectories[0][5].copy()
    for trajectory in trajectories[1:]:
        replace = trajectory[3] < first_times
        first[replace], first_times[replace] = trajectory[1][replace], trajectory[3][replace]
        present = trajectory[6] > 0
        terminal[present] = trajectory[2][present]
        terminal_counts += trajectory[6]
        periodic_counts += trajectory[5]
    require(np.all(terminal_counts == 1) and np.all(periodic_counts > 0), 'Missing/duplicate whole-cohort trajectories')
    both = 0
    for left in range(ranks):
        for right in range(left+1, ranks):
            lower, upper = trajectories[left], trajectories[right]
            common = (lower[5] > 0) & (upper[5] > 0)
            # This moving -c case has monotonically increasing local particle z.
            require(np.all(lower[4][common] < upper[3][common]), 'Overlapping MPI ownership in trajectory time')
            both += int(np.count_nonzero(common))
    log = (root/'moving/aprl.log').read_text()
    release_match = re.search(r'Frozen lab-plane injection: activations=(\d+), cross-rank activation pushes=(\d+), still frozen=(\d+)', log)
    moving_match = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)', log)
    cleanup = re.search(r'Residual CPML-carrier cleanup:.*z\+=(\d+) \(Q=([^ ]+) C\)', log)
    require(release_match and moving_match and cleanup, 'Missing lifecycle counters')
    activations, cross_rank, held = map(int, release_match.groups())
    steps, exits, cropped = map(int, moving_match.groups())
    cleanup_count, cleanup_charge = int(cleanup[1]), float(cleanup[2])
    require(activations == count and held == 0 and exits == count and cropped+cleanup_count == count,
            'Release/exit/carrier lifecycle count does not close')
    csv = {}
    for path in sorted((root/'moving/reference').glob('window-rank-*.csv')):
        rows = np.genfromtxt(path, delimiter=',', names=True)
        require(len(rows) == steps+1 and np.array_equal(rows['step'][:-1], np.arange(1, steps+1)), 'Window CSV step sequence')
        require(rows[-1].tobytes() == rows[-2].tobytes(), 'Window CSV duplicate final summary differs')
        require(all(np.all(np.diff(rows[key]) >= 0) for key in ('particle_exits', 'cropped_carriers')), 'Window counters decreased')
        csv[path.name] = {key: float(rows[key][-1]) for key in rows.dtype.names}
    require(len(csv) == ranks, 'Missing window rank CSV')
    totals = {key: sum(row[key] for row in csv.values()) for key in (
        'inserted_field_J', 'cropped_field_J', 'particle_exits', 'exit_charge_C', 'cropped_carriers', 'cropped_charge_C')}
    charge = float(np.sum(np.asarray(planes['release']['charge_C'], dtype=LD)))
    # Charge counters use sequential double sums; the cleanup log prints 10 digits.
    charge_error = abs(totals['exit_charge_C']-charge)
    charge_bound = max(abs(charge)*8*count*np.finfo(float).eps, 1e-30)
    cleanup_error = abs(totals['cropped_charge_C']+cleanup_charge-charge)
    require(charge_error <= charge_bound and cleanup_error <= charge_bound+abs(cleanup_charge)*1e-9,
            'Lifecycle charge accounting does not close within summation/log precision')
    with h5py.File(root/'moving/ledger/energy-ledger.h5') as stream:
        energy = stream['energy_ledger/records'][:]
    require(np.array_equal(energy['step'], np.arange(steps+1)), 'Energy ledger steps')
    require(energy['active_macroparticles'][-1] == 0 and energy['removed_macroparticles'][-1] == count,
            'Final active/removed cohort mismatch')
    require(np.all(energy['closure_valid'] == 0), 'Frozen moving-volume ledger must not claim exact closure')
    launch = json.loads((root/'moving/aprl.log.launch.json').read_text())
    require(launch['exit_code'] == 0 and 'Configured stop reached:' in log and 'phase=mpi-finalize-complete' in log,
            'Run did not complete normally')
    with h5py.File(root/'inputs/held.h5') as stream:
        source = stream['particles/records'][:]
    source = source[np.argsort(source['source_id'])]
    require(np.array_equal(source['source_id'], ids), 'Frozen input cohort mismatch')
    geometry = re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+).*fixed_lab_particle_exit_m=([^ ]+)', log)
    require(geometry, 'Actual moving geometry missing')
    summary = {'hdf5': files, 'ranks': ranks, 'source_count': count, 'configured_stop_complete': True,
        'particle_planes_counts': {name: len(rows) for name, rows in planes.items()},
        'trajectory': {'periodic_records': int(periodic_counts.sum()), 'terminal_records': int(terminal_counts.sum()),
            'terminal_type': 'outer_domain_exit', 'terminal_face': 'upper_z',
            'IDs_on_multiple_rank_pairs': both, 'per_rank_time_intervals_disjoint': True},
        'source_release_xy': array_error(source['plane_position_m'], planes['release']['position_m'][:, :2]),
        'source_release_u': array_error(source['proper_velocity'], planes['release']['proper_velocity']),
        'lifecycle': {'steps': steps, 'activations': activations, 'still_held': held,
            'cross_rank_activation_pushes': cross_rank, 'active_particle_exits': exits,
            'cropped_carriers': cropped, 'outer_carrier_cleanup': cleanup_count, 'remaining_carriers_by_count_closure': 0,
            'exit_charge_error_C': charge_error, 'exit_charge_roundoff_bound_C': charge_bound,
            'crop_plus_cleanup_charge_error_C': cleanup_error, 'cleanup_log_charge_C': cleanup_charge},
        'window_csv_final': csv, 'window_cumulative_totals': totals, 'window_csv_final_summary_duplicates_last_step': True,
        'actual_geometry': {'reference_nz': int(geometry[1]), 'active_nz': int(geometry[2]),
            'initial_offset': int(geometry[3]), 'safe_exit_lab_z_m': float(geometry[4])},
        'closure_valid': False, 'launcher': launch,
        'execution_evidence': execution_evidence(root, ranks, plan['threads'])}
    return summary, planes, first, terminal, periodic_counts, source


def compare_fields(reference, candidate, cohort_end):
    with h5py.File(reference/'moving/detectors/exit-field.h5') as af, h5py.File(candidate/'moving/detectors/exit-field.h5') as bf:
        a, b = af['field_plane'], bf['field_plane']
        times = a['time_s'][:]
        require(np.array_equal(times, b['time_s'][:]), 'Field clocks differ; no resampling allowed')
        keys = ('plane_z_m', 'x_first_m', 'y_first_m', 'dx_m', 'dy_m', 'nx', 'ny', 'magnetic_time_stagger', 'spatial_sampling', 'frame')
        require(all(np.array_equal(a.attrs[key], b.attrs[key]) for key in keys), 'Detector geometry/native phase differs')
        shape = a['fields'].shape
        require(shape == b['fields'].shape, 'Field output shape differs')
        before = times <= cohort_end
        stats = {field: {period: {'n': LD(0), 'd': LD(0), 'nc': np.zeros(3, dtype=LD),
            'dc': np.zeros(3, dtype=LD), 'max': 0., 'values': 0, 'equal': True, 'trace': []}
            for period in ('all', 'through_cohort')} for field in ('electric_V_per_m', 'magnetic_T')}
        for index in range(shape[0]):
            x, y = a['fields'][index], b['fields'][index]
            for field, periods in stats.items():
                xx = np.asarray(x[field], dtype=LD)
                difference = np.asarray(y[field], dtype=LD)-xx
                norm, delta = np.sum(xx*xx, axis=(0, 1)), np.sum(difference*difference, axis=(0, 1))
                for period, item in periods.items():
                    if period == 'through_cohort' and not before[index]:
                        continue
                    item['n'] += norm.sum(); item['d'] += delta.sum()
                    item['nc'] += norm; item['dc'] += delta
                    item['max'] = max(item['max'], float(abs(difference).max()))
                    item['values'] += difference.size
                    item['equal'] &= x[field].tobytes() == y[field].tobytes()
                    item['trace'].append((index, float(norm.sum()), float(delta.sum())))
        result = {}
        for field, periods in stats.items():
            result[field] = {}
            for period, item in periods.items():
                floor = max(row[1] for row in item['trace'])*1e-12
                instant = [(i, np.sqrt(d/n)) for i, n, d in item['trace'] if n > floor]
                worst = max(instant, key=lambda pair: pair[1]) if instant else (None, None)
                result[field][period] = {'relative_l2': float(np.sqrt(item['d']/item['n'])),
                    'component_relative_l2': [float(np.sqrt(d/n)) if n else None for d, n in zip(item['dc'], item['nc'])],
                    'max_abs': item['max'], 'rms_abs': float(np.sqrt(item['d']/item['values'])),
                    'bitwise_equal': bool(item['equal']), 'sample_count': len(item['trace']),
                    'scalar_values': item['values'], 'max_instantaneous_relative_l2': worst[1],
                    'max_instantaneous_sample': worst[0], 'instantaneous_denominator_floor_relative_to_peak_norm_squared': 1e-12}
        attrs = {key: (a.attrs[key].decode() if isinstance(a.attrs[key], bytes) else
                 a.attrs[key].item() if isinstance(a.attrs[key], np.generic) else a.attrs[key]) for key in keys}
    return {'shape': list(shape), 'time_bitwise_equal': True, 'geometry_native_phase_equal': True,
        'time_range_lab_s': [float(times[0]), float(times[-1])],
        'lab_spacing_min_max_s': [float(np.diff(times).min()), float(np.diff(times).max())],
        'attributes': attrs, 'errors': result,
        'scope': 'All exported raw lab-plane pixels/times, E/B in separate SI units. Not a full-volume or radiation-only field.'}


def nested_diff(a, b, prefix=''):
    if isinstance(a, dict) and isinstance(b, dict):
        return sum((nested_diff(a.get(key), b.get(key), prefix+'.'+str(key)) for key in sorted(set(a)|set(b))), [])
    return [] if a == b else [{'path': prefix, 'reference': a, 'candidate': b}]


def function_body(text, signature):
    start = text.index(signature)
    brace = text.index('{', start)
    depth = 1
    for index in range(brace+1, len(text)):
        depth += (text[index] == '{')-(text[index] == '}')
        if depth == 0:
            return text[brace:index+1]
    raise RuntimeError('Unclosed function '+signature)


def frozen_identity(a, b):
    cards = [yaml.safe_load((root/'moving/run.yaml').read_text()) for root in (a, b)]
    diff = nested_diff(*cards)
    require({row['path'] for row in diff} <= {'.beam.input.file', '.mesh.moving_window.reference_directory'}, 'Physical card settings differ')
    manifests = {}
    for root in (a, b):
        frozen = json.loads((root/'frozen.json').read_text())
        failures = [path for path, digest in frozen.items() if not Path(path).is_file() or sha(path) != digest]
        require(not failures, 'Frozen input/source hash mismatch: '+str(failures))
        manifests[str(root)] = {'sha256': sha(root/'frozen.json'), 'verified_files': len(frozen)}
    names = sorted(set(str(p.relative_to(a/'source')) for p in (a/'source').rglob('*') if p.is_file()) |
                   set(str(p.relative_to(b/'source')) for p in (b/'source').rglob('*') if p.is_file()))
    source_hashes = {name: {label: sha(root/'source'/name) if (root/'source'/name).is_file() else None
        for label, root in (('reference', a), ('candidate', b))} for name in names}
    changed = [name for name, hashes in source_hashes.items() if hashes['reference'] != hashes['candidate']]
    preserved = ['eb_particles.cpp', 'eb_particles.h', 'eb_deposition.cpp', 'eb_deposition.h',
        'eb_sources.cpp', 'eb_sources.h', 'eb_cpml.cpp', 'eb_cpml.h', 'boostframe.cpp', 'boostframe.h',
        'eb_mpi.cpp', 'eb_mpi.h', 'blocked_cowan_reference.cpp', 'blocked_cowan_reference.h',
        'moving_window_fields.cpp', 'moving_window_fields.h']
    require(all('src/'+name not in changed for name in preserved), 'Core push/deposition/field communication/reference changed')
    fields = [(root/'source/src/eb_field.cpp').read_text() for root in (a, b)]
    functions = ('Double EBFieldGrid::sampleComponent(', 'void EBFieldGrid::advanceMagneticCowan(',
                 'void EBFieldGrid::advanceMagneticYee(', 'void EBFieldGrid::advanceElectric(')
    unchanged = {name: function_body(fields[0], name) == function_body(fields[1], name) for name in functions}
    require(all(unchanged.values()), 'Original field/interpolation formula changed')
    prepared = {}
    for name in ('rho.bin', 'checkpoints.bin', 'phi.bin', 'metadata.txt'):
        hashes = {label: sha(root/'moving/reference'/name) for label, root in (('reference', a), ('candidate', b))}
        equal = hashes['reference'] == hashes['candidate']
        if name != 'metadata.txt':
            require(equal, 'Prepared native scalar file changed: '+name)
        prepared[name] = dict(hashes, bitwise_equal=equal)
    metadata = []
    for root in (a, b):
        raw = (root/'moving/reference/metadata.txt').read_bytes()
        body, checksum = raw.rsplit(b'checksum ', 1)
        value = 14695981039346656037
        for byte in body:
            value = ((value ^ byte)*1099511628211) & ((1 << 64)-1)
        require(value == int(checksum.strip(), 16), 'Prepared metadata checksum mismatch')
        metadata.append({line.split()[0]: line.split()[1:] for line in raw.decode().splitlines()})
    metadata_difference = []
    for key in sorted(set(metadata[0]) | set(metadata[1])):
        x, y = metadata[0].get(key), metadata[1].get(key)
        if x == y:
            continue
        item = {'key': key, 'reference': x, 'candidate': y}
        if key.endswith('_bits') and x and y and len(x) == len(y):
            item['binary64_values'] = [{'reference': struct.unpack('>d', bytes.fromhex(v))[0],
                'candidate': struct.unpack('>d', bytes.fromhex(w))[0],
                'candidate_minus_reference_bit_pattern_integer': int(w, 16)-int(v, 16)} for v, w in zip(x, y)]
        metadata_difference.append(item)
    require({row['key'] for row in metadata_difference} <= {
        'beta_bits', 'translation_bits', 'magnetic_interpolation_bits', 'checksum'},
        'Unexpected prepared metadata contract change')
    return {'card_only_path_differences': diff, 'manifests': manifests, 'changed_source_files': changed,
        'source_hashes': source_hashes, 'preserved_core_files': preserved, 'original_function_bodies_identical': unchanged,
        'prepared_files': prepared, 'prepared_metadata_checksums_valid': True,
        'prepared_metadata_differences': metadata_difference,
        'prepared_metadata_scope': 'Every difference is reported. MPI partitioning can change the long-double mass/velocity reduction before its binary64 cast. This is not a claim that the metadata contracts are bitwise equal, nor an extra physical tolerance.',
        'scope': 'R3 adds moving-only neighbor values for the existing particle interpolator and its step lifetime. R1 to R2 separately fixed directory visibility. Original Cowan/Boris/deposition/native interpolation weights are retained.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--unit-root', type=Path)
    parser.add_argument('--prior-review', type=Path)
    parser.add_argument('--field-relative-tolerance', type=float, default=2e-12,
                        help='Existing MPI module roundoff comparison bound; not a physical radiation tolerance')
    args = parser.parse_args()
    a, b, output = args.reference.resolve(), args.candidate.resolve(), args.output.resolve()
    require(not output.exists(), 'Refuse to overwrite review; use a new output')
    print('Inspecting reference complete output', flush=True)
    sa, pa, ia, ta, na, source_a = inspect_case(a)
    print('Inspecting candidate complete output', flush=True)
    sb, pb, ib, tb, nb, source_b = inspect_case(b)
    require(source_a.tobytes() == source_b.tobytes(), 'Frozen fixture records differ')
    require(sa['actual_geometry'] == sb['actual_geometry'], 'Actual moving geometry differs')
    require(np.array_equal(na, nb), 'Per-ID periodic trajectory sample counts differ')
    comparison = {'planes': {}}
    for name in pa:
        for key in IDENTITY:
            require(np.array_equal(pa[name][key], pb[name][key]), 'Cross-layout identity changed')
        comparison['planes'][name] = {'delta_K': weighted_energy(pa[name], pb[name]),
            **{key: array_error(pa[name][key], pb[name][key]) for key in ('position_m', 'proper_velocity', 'time_s')}}
    for name, first, second in (('first_released_periodic', ia, ib), ('terminal_events', ta, tb)):
        comparison[name] = {key: array_error(first[key], second[key]) for key in
            ('source_id', 'particle_id', 'charge_C', 'weight', 'position_m', 'proper_velocity', 'time_s')}
    print('Comparing every exported field sample', flush=True)
    comparison['field_plane'] = compare_fields(a, b, min(pa['exit']['time_s'].max(), pb['exit']['time_s'].max()))
    field_pass = all(item['all']['relative_l2'] <= args.field_relative_tolerance
        for item in comparison['field_plane']['errors'].values())
    bitwise_particles = all(item['proper_velocity']['bitwise_equal'] for item in comparison['planes'].values())
    print('Verifying frozen files and unchanged formulas', flush=True)
    identity = frozen_identity(a, b)
    status = ('PASS_SCOPED_BITWISE_PARTICLE_PLANES_AND_ROUNDOFF_FIELDS' if bitwise_particles and field_pass else
              'FUNCTIONAL_PASS_SCIENTIFIC_DIFFERENCES_REPORTED')
    report = {'format': 'aprl-moving-window-mpi-review-v2', 'status': status,
        'reference': str(a), 'candidate': str(b), 'cases': {'reference': sa, 'candidate': sb},
        'comparison': comparison, 'frozen_identity': identity,
        'acceptance': {'complete_outputs_and_lifecycle': True, 'all_plane_proper_velocities_bitwise_equal': bitwise_particles,
            'field_relative_l2_within_existing_module_roundoff_bound': field_pass,
            'field_relative_l2_bound': args.field_relative_tolerance,
            'particle_energy_threshold': 'No threshold introduced; all per-ID results reported, bitwise equality stated only when true.'},
        'method': {'script': str(Path(__file__).resolve()), 'sha256': sha(__file__),
            'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
            'energy': 'long-double ((u_candidate-u_reference)*(u_candidate+u_reference)).sum/(gamma_candidate+gamma_reference)*mc2',
            'trajectory': 'Bounded chunks, all records finite and identity checked; per-rank ID time intervals ordered, one terminal event per ID.'},
        'limits': ['One short functional fixture; not long-duration/production/radiation-only acceptance.',
            'Reference and candidate use different MPI/thread layouts; elapsed-time ratios are not same-MPI implementation speedups.',
            'Frozen-reference versus full fixed-window approximation remains separate from MPI consistency.',
            'Stationary energy-ledger closure_valid=false; inserted/cropped energy is volume bookkeeping, not radiation.',
            'Formal trajectories omit held states. First periodic sample means first released recorded state, not t_box=0.',
            'Initial left-entry witness labels do not prove entry coverage; check actual geometry separately.']}
    if args.unit_root:
        units = {}
        for ranks in (1, 2, 4):
            path = args.unit_root/f'ranks-{ranks}.json'
            unit = json.loads(path.read_text())
            require(unit['status'] == 'PASS' and unit['ranks'] == ranks, 'Unit result incomplete')
            require(unit['maximum_E_relative_error'] <= unit['relative_tolerance'] and
                    unit['maximum_B_relative_error'] <= unit['relative_tolerance'], 'Unit field threshold failed')
            if ranks > 1:
                require(unit['negative_control_samples'] > 0 and unit['legacy_missing_plane_E_error'] > 1e-8 and
                        unit['legacy_missing_plane_B_error'] > 1e-8, 'Unit negative control absent')
            units[str(ranks)] = {'path': str(path.resolve()), 'sha256': sha(path), 'result': unit}
        report['native_neighbor_unit_chain'] = units
    if args.prior_review:
        prior = json.loads(args.prior_review.read_text())
        rms = prior['comparison']['planes']['exit']['weighted_rms_energy_difference_eV']
        report['prior_failed_cross_layout_comparison'] = {'path': str(args.prior_review.resolve()),
            'sha256': sha(args.prior_review), 'status': prior['status'], 'exit_delta_K_rms_eV': rms,
            'current_exit_delta_K_rms_eV': comparison['planes']['exit']['delta_K']['rms_eV'],
            'scope': 'The prior unresolved discrepancy remains recorded; it is not reclassified by relaxing its criteria.'}
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps({'output': str(output), 'sha256': sha(output), 'status': status,
        'plane_energy_differences': {name: item['delta_K'] for name, item in comparison['planes'].items()},
        'field_relative_l2': {name: item['all']['relative_l2'] for name, item in comparison['field_plane']['errors'].items()}}, indent=2))


if __name__ == '__main__':
    main()
