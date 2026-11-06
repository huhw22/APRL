"""Prepare nine single-thread controls for exact MPI integer-slab balancing.

No job is launched and no reference directory is precreated. Within each
physical family the parent, candidate-off and candidate-on cards have identical
input events, geometry, timestep and diagnostics. Only the candidate-on card
contains the balancing configuration. The shared three-point source shape and
native substep diagnostics are explicitly recorded as common baseline choices.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re

import h5py
import yaml


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
BALANCE = dict(enabled=True, interval_steps=5, minimum_cells=3,
               minimum_improvement=.1, frozen_weight=.2, field_cell_weight=0.)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def rewrite(value, source, destination):
    if isinstance(value, dict):
        return {key: rewrite(item, source, destination) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, source, destination) for item in value]
    return value.replace(str(source), str(destination)) if isinstance(value, str) else value


def detector_headers(directory, card, expected_particles):
    result = []
    for kind in ['particle_planes', 'field_planes', 'energy_planes']:
        for plane in card['detectors'].get(kind, []):
            filename = plane['name'] + ('-energy.h5' if kind == 'energy_planes' else '.h5')
            path = directory / 'detectors' / filename
            with h5py.File(path, 'r') as handle:
                group = handle[next(iter(handle))]
                assert int(group['complete'][()]) == 1, str(path)
                if kind == 'field_planes':
                    values, counter = 'time_s', 'committed_samples'
                else:
                    values, counter = 'records', 'committed_records'
                count = len(group[values])
                assert int(group[counter][()]) == count, str(path)
                if kind == 'particle_planes':
                    assert count == expected_particles, str(path)
                result.append(dict(kind=kind, name=plane['name'], z_m=float(plane['z']),
                    filename=filename, count=count, count_field=values,
                    closed_baseline_file=str(path), baseline_file_bytes=path.stat().st_size))
    return result


def main(workspace, output, binary, build_record, ranks=32):
    workspace, output, binary, build_record = [Path(value).resolve()
        for value in [workspace, output, binary, build_record]]
    assert ranks > 0
    parent_binary = workspace / 'build/binomial3-fel-20261106/aprl'
    parent_record_path = workspace / 'results/binomial3_fel_20261106/build-r1/COMPLETE.json'
    parent_record = json.loads(parent_record_path.read_text())
    parent_sha = sha(parent_binary)
    assert parent_record.get('passed') is True and parent_sha == parent_record['binary_sha256']
    families = [
        ('drift', workspace / 'results/longitudinal_main_flow_20261103/flow-r2/drift-linear',
         31675, 257),
        ('undulator', workspace / 'results/longitudinal_main_flow_20261103/flow-r2/undulator-linear',
         31675, 257),
        ('fel7-w64', workspace / 'results/physical_electron_acceptance_20261106/weight-scan-r2/w64-a',
         178725, 379),
    ]
    prepared = []
    source_inputs = {}
    # Validate all original geometries before creating any case. No padding or
    # changed release/reference clock is allowed to make a requested rank fit.
    for family, original, particles, steps in families:
        card_path = original / 'run.yaml'
        card = yaml.safe_load(card_path.read_text())
        match = re.search(r'Mesh contract: cells=\(\s*(\d+),\s*(\d+),\s*(\d+)\)',
                          (original / 'aprl.log').read_text())
        assert match, str(original)
        active = [int(value) for value in match.groups()]
        assert active[2] >= 8 * ranks, (
            f'{family} has{active[2]} active z cells; parent needs at least8 per rank. '
            f'Requested{ranks} ranks would require geometry changes, which are forbidden.')
        assert card['mesh']['particle_coupling'] == 'conservative-source'
        assert card['energy_ledger']['enabled'] and card['energy_ledger']['discrete_operator_audit']
        source = Path(card['beam']['input']['file'])
        with h5py.File(source, 'r') as handle:
            assert len(handle['particles/records']) == particles
        if str(source) not in source_inputs:
            source_inputs[str(source)] = sha(source)
        headers = detector_headers(original, card, particles)
        prepared.append(dict(family=family, original=original, card=card,
            original_card_sha256=sha(card_path), active=active,
            expected_particles=particles, expected_steps=steps, detector_headers=headers))
    output.mkdir(parents=True, exist_ok=False)
    cases = []
    frozen_cards = {}
    group_records = []
    for family in prepared:
        normalized = []
        names = []
        for mode, suffix, program in [('old', 'old', parent_binary),
                                       ('off', 'newoff', binary), ('on', 'newon', binary)]:
            name = family['family'] + '-' + suffix
            directory = output / name
            directory.mkdir()
            card = rewrite(copy.deepcopy(family['card']), family['original'], directory)
            card['mesh']['longitudinal_shape'] = 'linear-binomial3'
            card['mesh'].pop('mpi_particle_balance', None)
            if mode == 'on':
                card['mesh']['mpi_particle_balance'] = copy.deepcopy(BALANCE)
            card['detectors']['substep_particle_planes'] = True
            for folder in ['detectors', 'ledger']:
                (directory / folder).mkdir()
            if card.get('trajectory', {}).get('enabled'):
                (directory / 'trajectory').mkdir()
            assert not (directory / 'reference').exists()
            card_path = directory / 'run.yaml'
            card_path.write_text(yaml.safe_dump(card, sort_keys=False))
            digest = sha(card_path)
            frozen_cards[str(card_path)] = digest
            canonical = rewrite(copy.deepcopy(card), directory, '@CASE')
            canonical['mesh'].pop('mpi_particle_balance', None)
            canonical_sha = hashlib.sha256(json.dumps(canonical, sort_keys=True,
                                                       allow_nan=False).encode()).hexdigest()
            normalized.append(canonical_sha)
            planes = [dict(name=plane['name'], z=float(plane['z']),
                           expected_records=family['expected_particles'])
                      for plane in card['detectors']['particle_planes']]
            expected = []
            for header in family['detector_headers']:
                row = dict(header)
                row.pop('closed_baseline_file')
                row.pop('baseline_file_bytes')
                expected.append(row)
                if header['kind'] == 'particle_planes':
                    expected.append(dict(kind='substep_particle_planes', name=header['name'],
                        z_m=header['z_m'], filename=header['name'] + '-substeps.h5',
                        count=family['expected_particles'], count_field='records'))
            source = card['beam']['input']['file']
            case = dict(name=name, family=family['family'], mode=mode,
                directory=str(directory), binary=str(program), nodes=1, ranks=ranks, threads=1,
                expected_particles=family['expected_particles'], expected_steps=family['expected_steps'],
                limit_seconds=600, case_limit_seconds=600, card_sha256=digest,
                normalized_physical_card_sha256=canonical_sha, input=source,
                input_sha256=source_inputs[source],
                represented_electrons=float(card['beam']['input']['electrons']),
                shape='linear-binomial3', active_grid=family['active'],
                reference_grid=card['mesh']['cells'], particle_planes=planes,
                expected_detector_files=expected,
                require_all_particles_through_planes=True,
                strict_ledger=True, discrete_operator_audit=True,
                original_card=str(family['original'] / 'run.yaml'),
                original_card_sha256=family['original_card_sha256'],
                balance=copy.deepcopy(BALANCE) if mode == 'on' else None,
                build_record=str(build_record) if mode != 'old' else str(parent_record_path.parent))
            save(directory / 'run-plan.json', case)
            cases.append(case)
            names.append(name)
        assert len(set(normalized)) == 1, 'A physical setting differs across the three versions'
        group_records.append(dict(family=family['family'], cases=names,
            normalized_physical_card_sha256=normalized[0],
            baseline_closed_headers=family['detector_headers'],
            expected_particles=family['expected_particles'], expected_steps=family['expected_steps'],
            active_grid=family['active']))
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(build_record), accepted_for_execution=False, cases=cases,
        families=group_records, parent_binary=str(parent_binary), parent_binary_sha256=parent_sha,
        source_inputs=source_inputs, frozen_cards=frozen_cards,
        nodes=1, ranks=ranks, threads=1, recommended_slurm_minutes=60,
        balance_configuration=BALANCE, submitted=False, physical_acceptance=False,
        common_baseline_changes=['All three versions explicitly use current linear-binomial3 shape',
                                 'All three versions enable native substep particle diagnostics'],
        unchanged=['Original source events and represented charge', 'Original reference and active geometry',
                   'Original dt/spacing/boost', 'Original release and stop', 'Original magnetic elements',
                   'Original CPML/fast reference/inlet', 'Original strict energy-ledger cadence'],
        comparison='Parent versus candidate-off requires unchanged science; candidate-off versus on permits only measured floating-point regrouping from exact integer ownership changes.',
        next_step='Root must validate candidate build, install/freeze runner and analysis, then explicitly accept this plan before submitting.')
    save(output / 'plan.json', plan)
    save(output / 'FROZEN_CARDS.json', frozen_cards)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
        accepted_for_execution=False, cases=len(cases), nodes=1, ranks=ranks, threads=1,
        no_geometry_padding=True, reference_directories_absent=True,
        source_inputs=source_inputs, parent_binary_sha256=parent_sha,
        normalized_physical_cards_identical_within_each_family=True))
    print(json.dumps(dict(complete=True, output=str(output), cases=len(cases), ranks=ranks,
                          threads=1, accepted_for_execution=False), indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--out', '--output', dest='output', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/smoke-r1')
    parser.add_argument('--binary', type=Path,
                        default=WORKSPACE / 'build/mpi-particle-balance-20261106/aprl')
    parser.add_argument('--build-record', '--buildrecord', dest='build_record', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/build-r1')
    parser.add_argument('--ranks', type=int, default=32)
    args = parser.parse_args()
    main(args.workspace, args.output, args.binary, args.build_record, args.ranks)
