"""Prepare a closed150-step prefix of the unchanged60M-electron FEL input.

No job is launched. Only frozen full-run cards/frame metadata and the closed
accepted particle source are read. Running full-run HDF5 outputs are not read.
The prefix is deliberately not a complete-beam or complete-FEL verification.
"""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import h5py
import yaml


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
SOURCE_SHA = '16be0bb3aff3d0d31e139a8d732db32fe237063fdd1d154fce2b17d9de49f125'
PREFIX_STEPS = 150


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def stat(path):
    value = Path(path).stat()
    return dict(bytes=value.st_size, mtime_ns=value.st_mtime_ns, inode=value.st_ino)


def rewrite(value, original, destination):
    if isinstance(value, dict):
        return {key: rewrite(item, original, destination) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, original, destination) for item in value]
    return value.replace(str(original), str(destination)) if isinstance(value, str) else value


def prefix_timing(card, frame):
    gb = float(card['mesh']['boost_gamma'])
    gamma_beta = math.sqrt(gb * gb - 1.)
    bb = gamma_beta / gb
    dz = float(card['mesh']['cell_size'][2])
    dt = dz / 299792458.
    center = float(card['beam']['reference']['initial_center_z'])
    origin = float(frame['predicted_frame_time_origin_lab_s'])
    initial = gb * gb * (center + bb * 299792458. * origin)
    advance = gamma_beta * dz
    physical_stop = float(card['stop']['z'])
    duration = (PREFIX_STEPS - .5) * dt
    accumulated = [0.]
    for unused in range(PREFIX_STEPS):
        accumulated.append(accumulated[-1] + dt)
    assert accumulated[-2] < duration < accumulated[-1]
    assert initial + PREFIX_STEPS * advance < physical_stop
    reference_nz = int(card['mesh']['cells'][2])
    active_nz = int(card['mesh']['moving_window']['cells'])
    reference_origin = float(card['mesh']['center'][2]) - .5 * reference_nz * dz
    release = float(frame['predicted_box_translation_m'])
    start = math.floor((release - reference_origin) / dz) - int(card['mesh']['moving_window']['lower_guard_cells'])
    duration_steps = math.ceil(duration / dt) + 1
    beta_low = float(frame['held_beta_min']) - 1e-12
    beta_high = float(frame['held_beta_max']) + 1e-12
    minimum_query = min(start - 1., start - (1. + beta_high) * duration_steps - 1.)
    maximum_query = max(start + active_nz, start - 1. - beta_low * duration_steps)
    assert start >= 1 and minimum_query >= 0 and maximum_query <= reference_nz - 1
    return dict(dt_box_s=dt, duration_box_s=duration, expected_steps=PREFIX_STEPS,
                configured_duration_step_upper=math.ceil(duration / dt),
                initial_reference_lab_z_m=initial, reference_advance_per_step_m=advance,
                physical_stop_unchanged=True, physical_stop_reference_z_m=physical_stop,
                reference_z_after_prefix_m=initial + PREFIX_STEPS * advance,
                duration_formula='(150-.5)*dt; keep original legal physical stop unchanged',
                repeated_addition_step149_time_s=accumulated[-2],
                repeated_addition_step150_time_s=accumulated[-1],
                expected_exit_code=1, expected_physical_complete=False,
                termination_protocol='Post-close duration exhaustion with all-rank readiness barrier; physical HDF5 complete0, no PIC_COMPLETE.',
                reference_query=dict(reference_nz=reference_nz, active_nz=active_nz,
                    reference_origin_z_m=reference_origin, start_cell=start,
                    duration_steps_for_reference=duration_steps,
                    conservative_minimum_query=minimum_query,
                    conservative_maximum_query=maximum_query,
                    available_minimum=0, available_maximum=reference_nz - 1,
                    passed=True))


def main(workspace, output, binary, build_record, ranks=128, field_cell_weight=0.):
    workspace, output, binary, build_record = [Path(value).resolve()
        for value in [workspace, output, binary, build_record]]
    original = workspace / 'results/binomial3_fel_20261106/full-r1/one4one'
    base = yaml.safe_load((original / 'run.yaml').read_text())
    original_plan = json.loads((original / 'run-plan.json').read_text())
    frame = json.loads((original / 'PREDICTED_FRAME.json').read_text())
    assert original_plan['expected_particles'] == 60360290 == frame['records']
    assert base['mesh']['cells'] == [80, 80, 4352]
    assert base['mesh']['moving_window']['cells'] == 942
    assert base['mesh']['longitudinal_shape'] == 'linear-binomial3'
    assert base['mesh']['particle_coupling'] == 'conservative-source'
    assert base['sources']['magnetic_elements'][0]['periods'] == 37
    assert ranks >= 1 and 942 // 3 >= ranks
    assert math.isfinite(field_cell_weight) and field_cell_weight >= 0.
    accepted = json.loads((build_record / 'COMPLETE.json').read_text())
    assert accepted.get('passed') is True
    binary_sha = sha(binary)
    assert binary_sha == accepted['binary_sha256']
    source = Path(base['beam']['input']['file'])
    source_before = stat(source)
    closed_source = json.loads((source.parent / 'COMPLETE.json').read_text())
    assert closed_source['complete'] and closed_source['electrons'] == 60360290
    assert closed_source['files']['particles.h5']['sha256'] == SOURCE_SHA
    assert sha(source) == SOURCE_SHA == original_plan['input_sha256']
    with h5py.File(source, 'r') as handle:
        assert int(handle['generation/complete'][()]) == 1
        records = handle['particles/records']
        assert len(records) == 60360290
        assert int(records[0]['source_id']) == 1 and int(records[-1]['source_id']) == 60360290
    timing = prefix_timing(base, frame)
    common = copy.deepcopy(base)
    common['mesh']['duration'] = timing['duration_box_s']
    common['stop'] = copy.deepcopy(base['stop'])
    common['detectors']['particle_planes'] = [plane for plane in common['detectors']['particle_planes']
        if plane['name'] in ['entry-particles', 'und-entry-particles']]
    assert len(common['detectors']['particle_planes']) == 2
    common['detectors']['substep_particle_planes'] = True
    static_threshold = 1. - .5 / ranks
    assert static_threshold > 1. - 1. / ranks
    field_cost = dict(requested_field_cell_weight=float(field_cell_weight),
        effective_field_cell_weight=float(field_cell_weight) if field_cell_weight > 0 else
            float(base['mesh']['cells'][0] * base['mesh']['cells'][1]) / 64.,
        zero_means='Automatic scheduling cost Nx*Ny/64 per longitudinal cell',
        scope='Scheduling model only; this value does not change fields, particles or physical equations.')
    output.mkdir(parents=True, exist_ok=False)
    cases, frozen_cards, canonical_hashes = [], {}, []
    for mode, suffix in [('off', 'static'), ('on', 'balanced')]:
        name = 'full-prefix-' + suffix
        directory = output / name
        directory.mkdir()
        card = rewrite(copy.deepcopy(common), original, directory)
        card['beam'] = copy.deepcopy(base['beam'])
        balance = dict(enabled=True, interval_steps=5, minimum_cells=3,
                       minimum_improvement=.1 if mode == 'on' else static_threshold,
                       frozen_weight=.2, field_cell_weight=float(field_cell_weight))
        card['mesh']['mpi_particle_balance'] = balance
        for folder in ['detectors', 'ledger']:
            (directory / folder).mkdir()
        assert not (directory / 'reference').exists()
        path = directory / 'run.yaml'
        path.write_text(yaml.safe_dump(card, sort_keys=False))
        digest = sha(path)
        frozen_cards[str(path)] = digest
        canonical = rewrite(copy.deepcopy(card), directory, '@CASE')
        canonical['mesh']['mpi_particle_balance']['minimum_improvement'] = '@CONTROL'
        canonical_sha = hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()
        canonical_hashes.append(canonical_sha)
        planes = [dict(name=plane['name'], z=float(plane['z']), expected_records=None,
                       partial_records=True) for plane in card['detectors']['particle_planes']]
        expected = []
        for plane in planes:
            for kind, suffix_name in [('particle_planes', '.h5'),
                                      ('substep_particle_planes', '-substeps.h5')]:
                expected.append(dict(kind=kind, name=plane['name'], z_m=plane['z'],
                    filename=plane['name'] + suffix_name, count=None, count_field='records',
                    minimum_count=0, maximum_count=60360290, partial_records=True))
        for key in ['field_planes', 'energy_planes']:
            for plane in card['detectors'].get(key, []):
                expected.append(dict(kind=key, name=plane['name'], z_m=float(plane['z']),
                    filename=plane['name'] + ('-energy.h5' if key == 'energy_planes' else '.h5'),
                    count=None, count_field='time_s' if key == 'field_planes' else 'records',
                    prefix_coverage=True))
        case = dict(name=name, family='full-prefix', mode=mode, directory=str(directory),
            binary=str(binary), binary_sha256=binary_sha, build_record=str(build_record),
            nodes=1, ranks=ranks, threads=1, expected_particles=60360290,
            expected_steps=PREFIX_STEPS, configured_duration_step_upper=timing['configured_duration_step_upper'],
            limit_seconds=1800, case_limit_seconds=1800, card_sha256=digest,
            normalized_physical_card_sha256=canonical_sha, input=str(source), input_sha256=SOURCE_SHA,
            represented_electrons=float(card['beam']['input']['electrons']),
            shape='linear-binomial3', active_grid=[80, 80, 942], reference_grid=[80, 80, 4352],
            particle_planes=planes, expected_detector_files=expected,
            require_all_particles_through_planes=False,
            bounded_duration_prefix=True, expected_run_exit_code=1,
            expected_physical_complete=False,
            verification_scope='closed150-step prefix, not full beam/FEL completion',
            expected_source_id_first=1, expected_source_id_last=60360290,
            input_source_ids_contiguous=True,
            partial_particle_checks=['Unique source_id per plane', 'IDs belong to the unchanged full input source',
                'Compare the exact per-plane source-ID sets between static and balanced controls'],
            strict_ledger=True, discrete_operator_audit=True, profile_interval=0,
            balance=balance, field_scheduling_cost=field_cost,
            expected_balance_repartitions=0 if mode == 'off' else 'at least one',
            original_card=str(original / 'run.yaml'), original_card_sha256=sha(original / 'run.yaml'),
            original_frame_metadata=str(original / 'PREDICTED_FRAME.json'),
            original_frame_metadata_sha256=sha(original / 'PREDICTED_FRAME.json'),
            prefix_timing=timing)
        save(directory / 'run-plan.json', case)
        cases.append(case)
    assert len(set(canonical_hashes)) == 1
    assert stat(source) == source_before
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        binary_sha256=binary_sha, build_record=str(build_record), accepted_for_execution=False,
        cases=cases, nodes=1, ranks=ranks, threads=1, recommended_slurm_minutes=65,
        profile_interval=0, source_inputs={str(source):SOURCE_SHA}, frozen_cards=frozen_cards,
        field_scheduling_cost=field_cost,
        require_all_particles_through_planes=False,
        bounded_duration_prefix=True, expected_run_exit_code=1,
        expected_physical_complete=False,
        verification_scope='closed150-step prefix, not full beam/FEL completion',
        partial_source_id_comparison_required=True, prefix_timing=timing,
        static_control_proof=dict(static_hysteresis=static_threshold,
            universal_nonnegative_gain_upper_bound=1. - 1. / ranks,
            statement='The static threshold exceeds any possible fractional maximum-load improvement; initial valid cuts are retained.'),
        common_prefix_changes=['Duration149.5dt with the original legal physical stop unchanged',
            'Particle output only entry and und-entry, with actual substep records'],
        unchanged=['All60360290 physical electron input events/IDs/weights',
            'Original center/frame origins', 'Original80x80/reference4352/active942 geometry',
            'Original spacing/dt/boost/37-period magnet', 'All original field and energy planes',
            'Conservative-source/three-point coupling, CPML and reference/inlet methods',
            'Strict ledger cadence'],
        source_provenance=dict(path=str(source), sha256=SOURCE_SHA,
            before=source_before, after=stat(source), closed_source_record=str(source.parent / 'COMPLETE.json'),
            verification='Accepted full-file SHA repeated; closed generator marker/count and ID endpoints checked. The accepted source verification establishes contiguous IDs1..N.'),
        physical_acceptance=False, submitted=False,
        next_step='Root must use a partial-prefix-aware reader, freeze the future runner and accept before submission. No complete-cohort detector claim is permitted.')
    save(output / 'plan.json', plan)
    save(output / 'FROZEN_CARDS.json', frozen_cards)
    save(output / 'PREFIX_TIMING.json', timing)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
        accepted_for_execution=False, cases=2, source_unchanged=True,
        expected_input_particles=60360290, expected_steps=PREFIX_STEPS,
        require_all_particles_through_planes=False,
        verification_scope='closed150-step prefix, not full beam/FEL completion',
        bounded_duration_prefix=True, expected_run_exit_code=1,
        expected_physical_complete=False,
        field_scheduling_cost=field_cost,
        reference_directories_absent=True))
    (output / 'README.md').write_text(
        '# Full physical-electron input: closed150-step prefix\n\n'
        'Prepared only. The full60,360,290-electron master, original37-period magnetic '
        'structure and80x80/reference4352/active942 geometry are unchanged. Both controls '
        'execute150 complete field steps using duration149.5dt while retaining the original '
        'legal physical stop. The expected exit is the known post-close duration exception '
        'and MPI_Abort1, after the all-rank output-close readiness barrier.\n\n'
        'This is a closed150-step prefix, not full beam/FEL completion. Only a subset of '
        'the full source has crossed entry/und-entry by this time. The future reader must '
        'check unique valid source IDs and equal per-plane ID sets between controls; it '
        'must not require or claim that all60M particles crossed either plane. Original '
        'field and energy planes remain configured, with prefix coverage only.\n\n'
        'All physical HDF5 completion flags and discrete-field CSV physical flags remain0. '
        'Coupling CSV markers are1 150 because that flag denotes closed records, not physical '
        'completion. Only the dedicated PREFIX_COMPLETE marker may report successful prefix '
        'validation, always with physical_complete=false; PIC_COMPLETE must not be created.\n\n'
        'Static and balanced both enable minimum3-cell support with interval5, frozen '
        f'weight0.2 and field weight{field_cell_weight:g} (effective {field_cost["effective_field_cell_weight"]:g}). '
        'Value0 means automatic Nx*Ny/64; positive values are scheduling costs only '
        'and do not alter physical equations. Static hysteresis1-.5/P guarantees no repartition; '
        'balanced uses0.1. One node128MPIx1, profile0,1800s per case and65min total are '
        'proposed. No running full-run HDF5 was read or modified.\n')
    print(json.dumps(dict(complete=True, output=str(output), cases=2, ranks=ranks,
                          expected_steps=PREFIX_STEPS, prefix_timing=timing), indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--out', '--output', dest='output', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/full-prefix-r2')
    parser.add_argument('--binary', type=Path,
                        default=WORKSPACE / 'build/mpi-particle-balance-r5-20261106/aprl')
    parser.add_argument('--build-record', '--buildrecord', dest='build_record', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/build-r6')
    parser.add_argument('--ranks', type=int, default=128)
    parser.add_argument('--field-cell-weight', type=float, default=0.)
    args = parser.parse_args()
    main(args.workspace, args.output, args.binary, args.build_record, args.ranks,
         args.field_cell_weight)
