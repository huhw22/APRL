#!/usr/bin/env python3
"""Additional checkpoint branch and nonempty-boundary-state gates; no auto-submit."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import socket
import struct
import time

import h5py
import numpy as np
import yaml

import run as common

HERE = Path(__file__).resolve().parent


def checkpoint(directory, step):
    return common.checkpoint_path(directory, step)


def settings(directory, steps):
    return {'directory': str(directory / 'checkpoints'), 'steps': steps,
            'io_timeout_seconds': 60}


def core_equal(first, second):
    """Numerical state, including local ordering, without native struct padding."""
    checked = []
    with h5py.File(first) as a, h5py.File(second) as b:
        for name in ('fields', 'cpml', 'particles', 'carriers'):
            if name not in a['blocks']:
                common.require(name not in b['blocks'], 'Boundary block presence differs')
                continue
            def visit(path, obj):
                if isinstance(obj, h5py.Dataset):
                    target = 'blocks/' + name + '/' + path
                    common.require(target in b, 'Missing restart numerical block ' + target)
                    common.bitwise_equal(obj[()], b[target][()], target)
                    checked.append(target)
            a['blocks'][name].visititems(visit)
    common.require(checked, 'No core snapshot state compared')
    return checked


def branch_suffix(cold, branch, cold_snapshot, actual_parent, step, reset_fields=False):
    """Compare a branch against the matching changed-model cold run.

    The cold snapshot provides exact prefix counts; the actual parent archive
    provides provenance. For reset fields sampled at every Maxwell step, the
    branch deliberately includes the snapshot-time sample once.
    """
    reference, trial = common.output_files(cold), common.output_files(branch)
    common.require(reference.keys() == trial.keys(), 'Branch output files differ')
    prefixes, cold_time = common.snapshot_metadata(cold_snapshot)
    _, parent_time = common.snapshot_metadata(actual_parent)
    common.require(cold_time == parent_time, 'Cold and source snapshot times differ')
    checked = []
    for name in reference:
        ga, aa, attrs_a = common.inspect_file(reference[name])
        gb, ab, attrs_b = common.inspect_file(trial[name])
        common.require(ga == gb and aa.keys() == ab.keys(), 'Branch output layout differs')
        prefix = prefixes[str(reference[name].resolve())]
        parent = attrs_b.get('restart_parent', '')
        if isinstance(parent, bytes):
            parent = parent.decode()
        scope = attrs_b.get('output_scope', '')
        if isinstance(scope, bytes):
            scope = scope.decode()
        common.require(Path(parent).resolve() == actual_parent.resolve() and
                       int(attrs_b.get('restart_step', -1)) == step and
                       float(attrs_b.get('restart_time_box_s', -1)) == parent_time and
                       scope.startswith('restart-segment'), 'Branch provenance differs')
        for attr in set(attrs_a) | set(attrs_b):
            if attr in common.IDENTITY_ATTRIBUTES:
                continue
            common.require(attr in attrs_a and attr in attrs_b, 'Missing branch science attribute')
            common.bitwise_equal(attrs_a[attr], attrs_b[attr], name + '@' + attr)
        for key, actual in ab.items():
            expected = aa[key]
            if np.ndim(actual) == 0:
                if key.startswith('committed_'):
                    continue
            elif key == 'records' and actual.dtype.names and 'step' in actual.dtype.names:
                expected = expected[expected['step'] > step]
            else:
                begin = int(prefix.get('validation_records', 0) if key == 'validation_records'
                            else prefix['records'])
                if reset_fields and key in ('time_s', 'fields'):
                    common.require(begin > 0, 'Reset-field cold reference has no boundary sample')
                    begin -= 1
                expected = expected[begin:]
            common.bitwise_equal(expected, actual, name + '/' + key)
            checked.append(name + '/' + key)
    return checked


def cpml_carrier_card():
    # This is the existing optimization_edges CPML exit geometry, with small
    # diagnostic cadences. At step1 all real particles enter the +x CPML but
    # still have nonzero current-only carriers, which must survive restart.
    card = common.gaussian_card()
    card['mesh'].update(cells=[12, 12, 80], cell_size=[.5, .5, .1], duration=30,
                        self_field_feedback=False, particle_steps_per_undulator_period=1)
    card['boundary']['cells'] = [2, 2, 4]
    card['initial_self_field'] = {'enabled': False}
    card['beam'] = {'reference': {'initial_center_z': -.5}, 'input': {
        'type': 'generated-gaussian', 'electrons': 8, 'macroparticles': 8,
        'gamma': 4, 'direction': [.2, 0, .9797958971132712], 'center': [1.95, 0, 0],
        'sigma_position': [0, 0, 0], 'sigma_proper_velocity': [0, 0, 0],
        'random_seed': 260126}}
    card['detectors'] = {'enabled': True, 'directory': 'detectors',
        'field_planes': [{'name': 'field', 'z': 1, 'rhythm': .4,
                          'buffer_samples': 4, 'compression': 0}],
        'particle_planes': [{'name': 'after-exit', 'z': 1,
                             'buffer_records': 16, 'compression': 0}]}
    card['stop'] = {'mode': 'reference-center-z', 'z': 2}
    card['trajectory'] = {'enabled': False}
    card['energy_ledger']['sample_interval_steps'] = 1
    return card


def tfsf_card():
    card = yaml.safe_load((HERE.parent / 'cards/detector-stop.yaml').read_text())
    card['mesh']['field_solver'] = 'yee'
    card['mesh']['duration'] = 30
    card['stop'] = {'mode': 'reference-center-z', 'z': 2}
    card['sources']['incident_waves'] = [{
        'profile': 'plane', 'position': [0, 0, 0], 'direction': [0, 0, 1],
        'polarization': [1, 0, 0], 'wavelength': 1,
        'peak_electric_field_V_per_m': 1000,
        'envelope': {'type': 'gaussian', 'center_time': 0, 'duration': 100,
                     'carrier_phase_rad': .3}}]
    card['energy_ledger'] = common.gaussian_card()['energy_ledger']
    return card


def magnet_card():
    card = common.gaussian_card()
    card['beam']['reference']['initial_center_z'] = -5
    card['beam']['input']['sigma_position'][2] = .5
    card['sources']['magnetic_elements'] = [{
        'type': 'uniform-dipole', 'entrance_z': 0, 'length': 10,
        'field_T': .1, 'polarization_angle_rad': .3}]
    return card


def saved_lab_z(path, gamma):
    """Decode only the documented geometry prefix, then use the saved frame."""
    with h5py.File(path) as file:
        header = file['header'][:]
        endian = '<' if int(header[9]) == 1 else '>'
        block = file['blocks/state/geometry']
        data = block['data'][:int(block['rank_counts'][0])].tobytes()
        version = struct.unpack_from(endian + 'I', data)[0]
        common.require(version == 1, 'Unsupported simulation geometry state version')
        # uint32 version, size_t cells[3], double extent[3], double origins[3].
        offset = 4 + 3 * int(header[5]) + 3 * 8
        _, lab_z_origin, box_z_origin = struct.unpack_from(endian + '3d', data, offset)
        z = file['blocks/particles/data'][:]['z']
    _, time_box = common.snapshot_metadata(path)
    beta = math.sqrt(1 - 1 / (gamma * gamma))
    # lab_time_origin matters for absolute time but not this lab-z transform.
    return lab_z_origin + gamma * (z - box_z_origin + beta * common.C * time_box)


def prepare(args):
    root, source = args.root.resolve(), args.science_root.resolve()
    common.require(not root.exists(), 'Use a fresh result directory')
    common.require((source / 'PASS').is_file() and
                   json.loads((source / 'RESULT.json').read_text())['status'] == 'PASS',
                   'Supply a completed checkpoint science campaign')
    parent = json.loads((source / 'plan.json').read_text())
    binary = args.binary.resolve() if args.binary else Path(parent['binary'])
    common.require(binary.is_absolute() and binary.is_file() and os.access(binary, os.X_OK),
                   'Candidate executable is missing')
    plan = {key: parent[key] for key in ('launcher', 'ranks', 'threads', 'timeout_seconds')}
    plan.update(format='aprl-checkpoint-branches-v1', binary=str(binary),
                binary_sha256=common.digest(binary), science_root=str(source), cases={})
    root.mkdir(parents=True)
    held = copy.deepcopy(parent['cases']['held-retirement']['base'])
    held['particle_retirement']['length'] = 12
    gaussian = copy.deepcopy(parent['cases']['gaussian-cpml']['base'])
    gaussian['detectors']['field_planes'][0].update(z=61, rhythm=.1)
    gaussian['detectors']['field_planes'].append({
        'name': 'new-field', 'z': 62, 'rhythm': .1, 'buffer_samples': 4, 'compression': 0})
    gaussian['detectors']['particle_planes'].append({
        'name': 'new-particles', 'z': 70, 'buffer_records': 16, 'compression': 0})
    gaussian['trajectory'] = {'enabled': False}
    reference = copy.deepcopy(parent['cases']['ballistic-validation']['base'])
    reference['detectors']['field_planes'][0]['rhythm'] = .05
    reference['trajectory'] = {'enabled': False}
    reference_source = copy.deepcopy(reference)
    reference_source['detectors']['field_planes'][0].pop('particle_background')
    automatic = yaml.safe_load((HERE.parent / 'cards/auto-longitudinal.yaml').read_text())
    automatic['energy_ledger'] = copy.deepcopy(held['energy_ledger'])
    automatic['energy_ledger']['sample_interval_steps'] = 2
    automatic['detectors']['field_planes'] = [{
        'name': 'field', 'z': .5, 'rhythm': .07, 'buffer_samples': 4, 'compression': 0}]
    definitions = {
        'retirement-before-entry': (held, [7, 55], source / 'held-retirement/continuous', 7),
        'reset-moved-new-planes': (gaussian, [7, 55], source / 'gaussian-cpml/continuous', 7),
        'new-ballistic-reference': (reference, [1, 17], root / 'reference-source', 1),
        'automatic-longitudinal': (automatic, [1, 3], None, 1),
        'cpml-carriers': (cpml_carrier_card(), [1, 7, 12], None, 1),
        'tfsf-incident': (tfsf_card(), [7], None, 7),
        'magnet-tail-head': (magnet_card(), [17, 55], None, 17),
    }
    common.write_card(root / 'reference-source', reference_source,
                      settings(root / 'reference-source', [1]))
    frozen_inputs = [source / 'plan.json', source / 'RESULT.json', source / 'PASS',
                     source / 'inputs/held.h5']
    for name, (card, steps, parent_directory, step) in definitions.items():
        cold = root / name / 'cold'
        common.write_card(cold, card, settings(cold, steps))
        actual_parent = checkpoint(parent_directory or cold, step)
        if parent_directory and str(parent_directory).startswith(str(source) + '/'):
            frozen_inputs.append(actual_parent)
        branch = root / name / 'restart'
        restart = settings(branch, [steps[-1]] if steps[-1] > step else [])
        restart['restart_from'] = str(actual_parent)
        restart['reset_diagnostics'] = name in ('reset-moved-new-planes', 'new-ballistic-reference')
        restart['allow_retirement_change'] = name == 'retirement-before-entry'
        common.write_card(branch, card, restart)
        plan['cases'][name] = {'base': card, 'steps': steps, 'restart_step': step,
                               'actual_parent': str(actual_parent),
                               'reset_fields': restart['reset_diagnostics']}
    common.save(root / 'plan.json', plan)
    frozen_inputs += [p for p in root.rglob('*') if p.is_file()]
    frozen_inputs += [HERE / 'run_branches.py', HERE / 'run.py']
    common.save(root / 'frozen.json', {str(p): common.digest(p) for p in frozen_inputs})
    print(root / 'plan.json')


def run(args):
    root = args.root.resolve()
    common.require(not (root / 'RUN_STARTED.json').exists(), 'Campaign already started')
    host = socket.gethostname().split('.')[0]
    common.require(os.environ.get('SLURM_JOB_ID') or (host == 'gnode1' and args.allow_gnode1),
                   'Use a real Slurm allocation or authorized gnode1 session')
    for path, expected in json.loads((root / 'frozen.json').read_text()).items():
        common.require(common.digest(Path(path)) == expected, 'Frozen input changed: ' + path)
    plan = json.loads((root / 'plan.json').read_text())
    common.require(common.digest(Path(plan['binary'])) == plan['binary_sha256'], 'Binary changed')
    common.save(root / 'RUN_STARTED.json', {'host': host, 'started_unix': time.time(),
        'slurm_job_id': os.environ.get('SLURM_JOB_ID')})
    common.launch(root, plan, root / 'reference-source')
    results = []
    for name, item in plan['cases'].items():
        cold, branch = root / name / 'cold', root / name / 'restart'
        step = item['restart_step']
        common.launch(root, plan, cold)
        actual_parent, cold_snapshot = Path(item['actual_parent']), checkpoint(cold, step)
        same_before = core_equal(actual_parent, cold_snapshot)
        state = common.inspect_snapshot(actual_parent)
        if name == 'retirement-before-entry':
            common.require(state['particles'] == 64 and state['retirement_carriers'] == 0,
                           'Retirement branch checkpoint already lost particles')
            with h5py.File(cold / 'ledger/energy-ledger.h5') as f:
                records = f['energy_ledger/records'][:]
                common.require(np.all(records[records['step'] <= step]['removed_macroparticles'] == 0),
                               'Retirement entered before the requested branch')
        elif name == 'new-ballistic-reference':
            prefixes, _ = common.snapshot_metadata(cold_snapshot)
            reference_prefixes = [p for p in prefixes.values()
                                  if p['group'].lstrip('/') == 'ballistic_reference']
            common.require(len(reference_prefixes) == 1 and
                           all(p['records'] == 0 for p in reference_prefixes),
                           'New reference branch is after an entrance crossing')
        elif name == 'cpml-carriers':
            common.require(state['cpml_carriers'] > 0 and state['particles'] == 0,
                           'CPML fixture must contain nonzero carriers and no live particles')
        elif name == 'tfsf-incident':
            with h5py.File(actual_parent) as f:
                common.require(np.any(f['blocks/fields/ex/data'][:] != 0),
                               'TF/SF fixture has no electric wave to restore')
        elif name == 'magnet-tail-head':
            lab_z = saved_lab_z(actual_parent, item['base']['mesh']['boost_gamma'])
            inside = int(np.count_nonzero((lab_z >= 0) & (lab_z <= 10e-6)))
            downstream = int(np.count_nonzero(lab_z > 10e-6))
            common.require(inside > 0 and downstream > 0 and np.all(lab_z >= 0),
                           'Magnet snapshot must straddle its interaction region and downstream')
            with h5py.File(cold / 'ledger/energy-ledger.h5') as file:
                records = file['energy_ledger/records'][:]
                previous = records[records['step'] <= step]
                common.require(len(previous) > 0 and
                               np.any(previous['prescribed_work_J'] != 0),
                               'Magnet fixture did not accumulate nonzero prescribed-device work')
            state.update(in_magnet=inside, downstream_of_magnet=downstream,
                         minimum_lab_z_m=float(np.min(lab_z)),
                         maximum_lab_z_m=float(np.max(lab_z)))
        common.launch(root, plan, branch)
        compared = branch_suffix(cold, branch, cold_snapshot, actual_parent, step,
                                  reset_fields=item['reset_fields'])
        core_after = []
        if item['steps'][-1] > step:
            core_after = core_equal(checkpoint(cold, item['steps'][-1]),
                                    checkpoint(branch, item['steps'][-1]))
        if name == 'automatic-longitudinal':
            cold_manifest = yaml.safe_load((cold / 'run-manifest.yaml').read_text())
            resumed_manifest = yaml.safe_load((branch / 'run-manifest.yaml').read_text())
            common.require(item['base']['mesh']['cells'][2] == 'auto' and
                           isinstance(cold_manifest['resolved_mesh']['cells'][2], int) and
                           cold_manifest['resolved_mesh'] == resumed_manifest['resolved_mesh'],
                           'Automatic resolved geometry was not retained exactly')
        results.append({'case': name, 'snapshot_coverage': state, 'datasets': compared,
                        'core_before': same_before, 'core_after': core_after})
    # At step12 the deterministic reference centre is already beyond1.25um;
    # this new stop is nevertheless downstream of the original field plane1um,
    # so YAML validation succeeds and the restart-specific stop gate is tested.
    stopped = copy.deepcopy(plan['cases']['cpml-carriers']['base'])
    stopped['stop']['z'] = 1.25
    directory = root / 'already-at-stop'
    common.write_card(directory, stopped, {'restart_from': str(
        checkpoint(root / 'cpml-carriers/cold', 12)), 'io_timeout_seconds': 60})
    common.launch(root, plan, directory, rejection=r'already satisfies the requested physical stop')
    common.save(root / 'RESULT.json', {'status': 'PASS', 'cases': results,
                'already_at_stop_rejected': True,
                'limitations': ['Small same-rank correctness only; no production I/O throughput claim']})
    (root / 'PASS').write_text('Checkpoint branch and boundary-state gates passed.\n')
    print(root / 'RESULT.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    prepare_parser = commands.add_parser('prepare')
    prepare_parser.add_argument('--root', type=Path, required=True)
    prepare_parser.add_argument('--science-root', type=Path, required=True)
    prepare_parser.add_argument('--binary', type=Path)
    run_parser = commands.add_parser('run')
    run_parser.add_argument('--root', type=Path, required=True)
    run_parser.add_argument('--allow-gnode1', action='store_true')
    args = parser.parse_args()
    try:
        prepare(args) if args.action == 'prepare' else run(args)
    except Exception as error:
        if args.action == 'run' and args.root.exists():
            common.save(args.root / 'FAIL.json', {'error': str(error), 'time_unix': time.time()})
        raise


if __name__ == '__main__':
    main()
