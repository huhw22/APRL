"""Prepare a paired cold-uniform-sphere ownership regression; never launch."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re

import h5py
import numpy as np
import yaml


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def rewrite(value, original, destination):
    if isinstance(value, dict):
        return {key: rewrite(item, original, destination) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, original, destination) for item in value]
    return value.replace(str(original), str(destination)) if isinstance(value, str) else value


def main(workspace, output, binary, build_record, ranks=64):
    workspace, output, binary, build_record = [Path(value).resolve()
        for value in [workspace, output, binary, build_record]]
    original = workspace / ('results/energy_ledger_test_record_20261030/physical_benchmarks/'
                             'space_charge_pic_20261030_r2/uniform')
    base = yaml.safe_load((original / 'run.yaml').read_text())
    reference = json.loads((original / 'case.json').read_text())
    source = Path(base['beam']['input']['file'])
    source_sha = sha(source)
    assert source_sha == reference['input_sha256']
    with h5py.File(source, 'r') as handle:
        records = handle['particles/records'][:]
    assert len(records) == 32768
    assert np.array_equal(records['source_id'], np.arange(1, 32769, dtype=np.uint64))
    assert np.all(np.isfinite(records['proper_velocity']))
    gamma = float(base['mesh']['boost_gamma'])
    dz = float(base['mesh']['cell_size'][2])
    dt = dz / 299792458.
    initial_reference_z = float(reference['geometry']['reference_initial_lab_z_m'])
    reference_advance = math.sqrt(gamma - 1.) * math.sqrt(gamma + 1.) * dz
    stop_steps_real = (float(base['stop']['z']) - initial_reference_z) / reference_advance
    expected_steps = math.ceil(stop_steps_real)
    assert expected_steps == 1338 == reference['geometry']['expected_stop_step']
    with h5py.File(original / 'ledger/energy-ledger.h5', 'r') as handle:
        ledger = handle['energy_ledger']
        assert int(ledger['complete'][()]) == 1
        assert int(ledger['records'][-1]['step']) == expected_steps
    log = (original / 'aprl.log').read_text()
    assert re.search(r'phase=final step=1338\b', log)
    assert base['mesh']['cells'] == [160, 160, 3328]
    assert base['mesh']['moving_window']['cells'] == 1332
    assert ranks >= 1 and 1332 >= 8 * ranks
    balance = dict(enabled=True, interval_steps=25, minimum_cells=2,
                   minimum_improvement=0., frozen_weight=.2, field_cell_weight=400.)
    theory_paths = [original / 'case.json', original / 'reference_input/expected_events.npz']
    theory = {str(path):sha(path) for path in theory_paths}
    output.mkdir(parents=True, exist_ok=False)
    cases = []
    cards = {}
    canonical_hashes = []
    for mode in ['off', 'on']:
        name = 'cold-uniform-' + mode
        directory = output / name
        directory.mkdir()
        card = rewrite(copy.deepcopy(base), original, directory)
        # The physical source is inside the old output tree. Preserve it
        # explicitly rather than rewriting it to a nonexistent new input.
        card['beam'] = copy.deepcopy(base['beam'])
        card['mesh']['particle_coupling'] = 'conservative-source'
        card['mesh']['longitudinal_shape'] = 'linear-binomial3'
        card['mesh'].pop('mpi_particle_balance', None)
        if mode == 'on':
            card['mesh']['mpi_particle_balance'] = copy.deepcopy(balance)
        card['detectors']['substep_particle_planes'] = True
        for folder in ['detectors', 'ledger']:
            (directory / folder).mkdir()
        if card.get('trajectory', {}).get('enabled'):
            (directory / 'trajectory').mkdir()
        assert not (directory / 'reference').exists()
        assert card['beam']['input']['file'] == str(source)
        assert card['mesh']['cells'] == base['mesh']['cells']
        assert card['mesh']['duration'] == base['mesh']['duration']
        assert card['stop'] == base['stop']
        path = directory / 'run.yaml'
        path.write_text(yaml.safe_dump(card, sort_keys=False))
        digest = sha(path)
        cards[str(path)] = digest
        canonical = rewrite(copy.deepcopy(card), directory, '@CASE')
        canonical['mesh'].pop('mpi_particle_balance', None)
        canonical_sha = hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()
        canonical_hashes.append(canonical_sha)
        planes = [dict(name=plane['name'], z=float(plane['z']), expected_records=32768)
                  for plane in card['detectors']['particle_planes']]
        expected = []
        for plane in planes:
            for kind, suffix in [('particle_planes', '.h5'),
                                  ('substep_particle_planes', '-substeps.h5')]:
                expected.append(dict(kind=kind, name=plane['name'], z_m=plane['z'],
                    filename=plane['name'] + suffix, count=32768, count_field='records'))
        case = dict(name=name, family='cold-uniform', mode=mode, directory=str(directory),
            binary=str(binary), build_record=str(build_record), nodes=1, ranks=ranks, threads=1,
            expected_particles=32768, expected_steps=expected_steps,
            configured_duration_step_upper=math.ceil(float(base['mesh']['duration']) / dt),
            limit_seconds=1800, case_limit_seconds=1800, card_sha256=digest,
            normalized_physical_card_sha256=canonical_sha, input=str(source), input_sha256=source_sha,
            represented_electrons=float(base['beam']['input']['electrons']), shape='linear-binomial3',
            active_grid=[160, 160, 1332], reference_grid=[160, 160, 3328],
            particle_planes=planes, expected_detector_files=expected,
            require_all_particles_through_planes=True, strict_ledger=True,
            discrete_operator_audit=True, profile_interval=0,
            balance=copy.deepcopy(balance) if mode == 'on' else None,
            original_card=str(original / 'run.yaml'), original_card_sha256=sha(original / 'run.yaml'),
            reference_model='Unchanged cold uniform sphere; theoretical shell events and original input unchanged.')
        save(directory / 'run-plan.json', case)
        cases.append(case)
    assert len(set(canonical_hashes)) == 1
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(build_record), accepted_for_execution=False, cases=cases,
        nodes=1, ranks=ranks, threads=1, recommended_slurm_minutes=65, profile_interval=0,
        source_inputs={str(source):source_sha}, frozen_cards=cards,
        theory_reference=theory, balance_configuration=balance,
        stop_derivation=dict(initial_reference_z_m=initial_reference_z,
            advance_per_field_step_m=reference_advance, continuous_stop_step=stop_steps_real,
            expected_steps=expected_steps, original_closed_ledger_steps=expected_steps,
            formula='ceil((stop_z-initial_reference_lab_z)/(sqrt(gamma-1)*sqrt(gamma+1)*dz))'),
        common_baseline_changes=['Both versions use conservative-source instead of original Boris',
            'Both versions use matched linear-binomial3 source and native substep planes'],
        unchanged=['Source HDF5/source IDs/charge', 'Theoretical shell model and expected_events',
            'Nx/Ny160/reference3328/active1332', 'Spacing/dt/boost/duration/stop',
            'All original particle planes', 'CPML and original supply options', 'Strict energy-ledger cadence'],
        comparison='Only new-off versus new-on is an algorithm-matched pair. The historical Boris result is not claimed to use the same pusher.',
        scope='Cold spherical space-charge regression for exact integer ownership changes; physical acceptance must be assessed separately.',
        physical_acceptance=False, submitted=False)
    save(output / 'plan.json', plan)
    save(output / 'FROZEN_CARDS.json', cards)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
        accepted_for_execution=False, cases=2, source_unchanged=True,
        theory_unchanged=True, geometry_unchanged=True, expected_steps=expected_steps,
        source_ids_exact=True, reference_directories_absent=True))
    (output / 'README.md').write_text(
        '# Cold uniform sphere: paired MPI ownership regression\n\n'
        'Prepared only; no simulation is launched. Both versions use the same original32,768 '
        'source events and1e8 represented electrons, with the original analytical shell '
        'reference and six lab planes. Geometry remains160x160, reference3328, active1332, '
        'with unchanged spacing, timestep, duration and stop.\n\n'
        'The completed original ledger and an independent reference-worldline calculation '
        'both give1338 actual steps; duration alone permits1339. Both new versions use '
        'conservative-source and matched three-point coupling. Therefore the matched '
        'algorithm comparison is new-off versus new-on, not new versus historical Boris.\n\n'
        'Balanced uses interval25, min2, zero improvement threshold and field-cell weight400. '
        'The positive field baseline keeps nearly empty field regions from being combined '
        'into a huge slab. Single-node64MPIx1;1800s per case and65min total are proposed.\n')
    print(json.dumps(dict(complete=True, output=str(output), cases=2,
                          expected_steps=expected_steps, ranks=ranks), indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--out', '--output', dest='output', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/cold-sphere-r1')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-record', '--buildrecord', dest='build_record', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=64)
    args = parser.parse_args()
    main(args.workspace, args.output, args.binary, args.build_record, args.ranks)
