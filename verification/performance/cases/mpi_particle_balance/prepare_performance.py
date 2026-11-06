"""Prepare paired MPI-balancing stress tests, without running or submitting.

The source is the unchanged seven-period w8 input. Both controls use the same
942-cell active window and4096-cell reference, explicitly as a short-bunch,
long-window performance stress test. This is not the37-period physical FEL
acceptance, and changing window size is not counted as a balancing speedup.
"""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np
import yaml


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
C = 299792458.
BLOCK = 131072


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


def held_beta_statistics(source, boost_gamma):
    """Read only a closed source; bound the C++ unweighted translation beta."""
    total = np.longdouble(0)
    minimum, maximum = math.inf, -math.inf
    gamma_beta = math.sqrt(boost_gamma - 1.) * math.sqrt(boost_gamma + 1.)
    with h5py.File(source, 'r') as handle:
        selection = handle['selection']
        assert int(selection.attrs['complete']) == 1
        data = handle['particles/records']
        count = len(data)
        for start in range(0, count, BLOCK):
            records = data[start:start + BLOCK]
            assert np.all(np.isfinite(records['proper_velocity']))
            u = records['proper_velocity'][:, 2]
            # The original loader takes longitudinal-only lab momentum,
            # boosts it, and forms u'_z/hypot(1,u'_z). This light-front form
            # is the accepted preparation expression for the same operation.
            light = (np.hypot(1., u) + u) / (boost_gamma + gamma_beta)
            held_u = .5 * (light - 1. / light)
            beta = held_u / np.hypot(1., held_u)
            assert np.all((beta >= 0) & (beta < 1))
            total += np.sum(beta.astype(np.longdouble))
            minimum = min(minimum, float(beta.min()))
            maximum = max(maximum, float(beta.max()))
    return dict(records=count, unweighted_mean_beta=float(total / count),
                minimum_beta=minimum, maximum_beta=maximum,
                numerical_bound_margin_beta=1e-12,
                note='Bounds with1e-12 beta margin cover roundoff in the C++ boost roundtrip and MPI sums; source events are unchanged.')


def reference_query_check(card, frame, beta, reference_nz=4096, active_nz=942):
    """Use initializeMovingWindow's exact index formulas and safe beta bounds."""
    mesh = card['mesh']
    dz = float(mesh['cell_size'][2])
    dt = dz / C
    gb = float(mesh['boost_gamma'])
    activation = float(card['beam']['reference']['activation_plane_z'])
    # Frozen injection sets the Lorentz lab-z origin to the activation plane.
    release = float(frame['predicted_box_translation_m']) + (activation - activation) / gb
    origin = float(mesh['center'][2]) - .5 * reference_nz * dz
    relative = (release - origin) / dz
    start = math.floor(relative) - int(mesh['moving_window']['lower_guard_cells'])
    fractional_distance = min(relative - math.floor(relative), math.ceil(relative) - relative)
    assert fractional_distance > 1e-6, 'Frame prediction too close to an integer cut for a safe static range check'
    duration_steps = math.ceil(float(mesh['duration']) / dt) + 1
    mean = beta['unweighted_mean_beta']
    last_query = start - (1. + mean) * duration_steps - 1.
    max_query = max(start + active_nz, start - 1. - mean * duration_steps)
    margin = beta['numerical_bound_margin_beta']
    low = beta['minimum_beta'] - margin
    high = beta['maximum_beta'] + margin
    conservative_min = min(start - 1., start - (1. + high) * duration_steps - 1.)
    conservative_max = max(start + active_nz, start - 1. - low * duration_steps)
    last_lab = max([activation] + [float(plane['z']) for key in
                   ['particle_planes', 'field_planes', 'energy_planes']
                   for plane in card['detectors'].get(key, [])])
    for magnet in card['sources']['magnetic_elements']:
        assert magnet['type'] == 'planar-undulator' and not magnet.get('gaussian_fringe', False)
        last_lab = max(last_lab, float(magnet['entrance_z']) +
                       float(magnet['period']) * int(magnet['periods']))
    minimum_cells = (math.ceil((last_lab - activation) / (gb * dz)) +
                     int(mesh['moving_window']['lower_guard_cells']) +
                     int(mesh['moving_window']['upper_guard_cells']) + 2)
    assert start >= 1 and conservative_min >= 0
    assert conservative_max <= reference_nz - 1
    assert active_nz >= minimum_cells
    return dict(passed=True, reference_nz=reference_nz, active_nz=active_nz,
                reference_origin_z_m=origin, activation_box_z_at_t0_m=release,
                release_relative_cell=relative, release_fractional_cell_clearance=fractional_distance,
                start_cell=start, duration_steps_for_reference_coverage=duration_steps,
                predicted_unweighted_translation_beta=mean,
                predicted_last_query=last_query, predicted_maximum_query=max_query,
                conservative_minimum_query=conservative_min,
                conservative_maximum_query=conservative_max,
                available_minimum=0, available_maximum=reference_nz - 1,
                minimum_active_cells_for_original_interactions=minimum_cells,
                formulas={'durationSteps': 'ceil(mesh.duration/dt)+1',
                          'startCell': 'floor((release-referenceOrigin)/dz)-lowerGuard',
                          'lastQuery': 'startCell-(1+translationBeta)*durationSteps-1',
                          'maximumQuery': 'max(startCell+activeCells,startCell-1-translationBeta*durationSteps)'},
                source='src/simulation_moving_window.cpp initializeMovingWindow',
                frame_prediction_scope='Accepted original source frame; unchanged center/source/extrema.1e-6-cell clearance guards integer rounding.')


def main(workspace, output, binary, build_record, ranks=128, field_cell_weight=0., repeats=2):
    workspace, output, binary, build_record = [Path(value).resolve()
        for value in [workspace, output, binary, build_record]]
    original = workspace / 'results/physical_electron_acceptance_20261106/weight-scan-r2/w8-a'
    original_plan = json.loads((original / 'run-plan.json').read_text())
    frame = json.loads((original / 'PREDICTED_FRAME.json').read_text())
    original_card = yaml.safe_load((original / 'run.yaml').read_text())
    source = Path(original_card['beam']['input']['file'])
    assert ranks >= 1 and 942 // 3 >= ranks
    assert math.isfinite(field_cell_weight) and field_cell_weight >= 0.
    assert repeats in [1, 2]
    assert original_plan['expected_particles'] == 1429799
    assert original_plan['expected_steps'] == 379
    source_sha = sha(source)
    assert source_sha == original_plan['input_sha256']
    beta = held_beta_statistics(source, float(original_card['mesh']['boost_gamma']))
    assert beta['records'] == original_plan['expected_particles']
    # Shared changes establish one identical stress-test geometry for both
    # versions. They are never credited as a gain from balancing.
    common = copy.deepcopy(original_card)
    common['mesh']['cells'][2] = 4096
    common['mesh']['moving_window']['cells'] = 942
    common['mesh']['longitudinal_shape'] = 'linear-binomial3'
    common['detectors']['particle_planes'] = [plane for plane in common['detectors']['particle_planes']
        if plane['name'] in ['entry-particles', 'und-entry-particles', 'p07-particles']]
    assert len(common['detectors']['particle_planes']) == 3
    query = reference_query_check(common, frame, beta)
    static_threshold = 1. - .5 / ranks
    max_possible_improvement = 1. - 1. / ranks
    assert static_threshold > max_possible_improvement
    base_balance = dict(enabled=True, interval_steps=5, minimum_cells=3,
                        minimum_improvement=.1, frozen_weight=.2,
                        field_cell_weight=float(field_cell_weight))
    field_cost = dict(requested_field_cell_weight=float(field_cell_weight),
        effective_field_cell_weight=float(field_cell_weight) if field_cell_weight > 0 else
            float(common['mesh']['cells'][0] * common['mesh']['cells'][1]) / 64.,
        zero_means='Automatic scheduling cost Nx*Ny/64 per longitudinal cell',
        scope='Scheduling model only; this value does not change fields, particles or physical equations.')
    output.mkdir(parents=True, exist_ok=False)
    cases = []
    frozen_cards = {}
    normalized = []
    order = [('perf-a', 'off', 'static'), ('perf-a', 'on', 'balanced')]
    if repeats == 2:
        order += [('perf-b', 'on', 'balanced'), ('perf-b', 'off', 'static')]
    for family, mode, suffix in order:
        name = family + '-' + suffix
        directory = output / name
        directory.mkdir()
        card = rewrite(copy.deepcopy(common), original, directory)
        balance = copy.deepcopy(base_balance)
        if mode == 'off':
            balance['minimum_improvement'] = static_threshold
        card['mesh']['mpi_particle_balance'] = balance
        for folder in ['detectors', 'ledger']:
            (directory / folder).mkdir()
        assert not (directory / 'reference').exists()
        card_path = directory / 'run.yaml'
        card_path.write_text(yaml.safe_dump(card, sort_keys=False))
        digest = sha(card_path)
        frozen_cards[str(card_path)] = digest
        canonical = rewrite(copy.deepcopy(card), directory, '@CASE')
        canonical['mesh']['mpi_particle_balance']['minimum_improvement'] = '@CONTROL'
        normalized_sha = hashlib.sha256(json.dumps(canonical, sort_keys=True,
                                                   allow_nan=False).encode()).hexdigest()
        normalized.append(normalized_sha)
        particles = [dict(name=plane['name'], z=float(plane['z']), expected_records=1429799)
                     for plane in card['detectors']['particle_planes']]
        expected = []
        for plane in particles:
            for kind, suffix_name in [('particle_planes', '.h5'),
                                      ('substep_particle_planes', '-substeps.h5')]:
                expected.append(dict(kind=kind, name=plane['name'], z_m=plane['z'],
                    filename=plane['name'] + suffix_name, count=1429799, count_field='records'))
        for key in ['field_planes', 'energy_planes']:
            for plane in card['detectors'][key]:
                expected.append(dict(kind=key, name=plane['name'], z_m=float(plane['z']),
                    filename=plane['name'] + ('-energy.h5' if key == 'energy_planes' else '.h5'),
                    count=380, count_field='time_s' if key == 'field_planes' else 'records'))
        case = dict(name=name, family=family, mode=mode, directory=str(directory),
            binary=str(binary), build_record=str(build_record), nodes=1, ranks=ranks, threads=1,
            expected_particles=1429799, expected_steps=379, limit_seconds=1200,
            case_limit_seconds=1200, card_sha256=digest,
            normalized_physical_card_sha256=normalized_sha, input=str(source), input_sha256=source_sha,
            represented_electrons=float(card['beam']['input']['electrons']),
            shape='linear-binomial3', active_grid=[80, 80, 942], reference_grid=[80, 80, 4096],
            particle_planes=particles, expected_detector_files=expected,
            require_all_particles_through_planes=True, strict_ledger=True,
            discrete_operator_audit=True, profile_interval=0,
            balance=balance, field_scheduling_cost=field_cost,
            expected_balance_repartitions=0 if mode == 'off' else 'at least one',
            original_card=str(original / 'run.yaml'), original_card_sha256=sha(original / 'run.yaml'),
            scenario='Seven-period short physical input in the same long field window for both controls; performance stress test only.')
        save(directory / 'run-plan.json', case)
        cases.append(case)
    assert len(set(normalized)) == 1, 'Static/balanced stress physics differs'
    static_proof = dict(ranks=ranks, static_hysteresis=static_threshold,
        universal_nonnegative_load_gain_upper_bound=max_possible_improvement,
        proof='For total W, old maximum M<=W and optimal new maximum Mnew>=W/P. Therefore(M-Mnew)/M<=1-1/P. Static threshold1-.5/P is strictly larger, so the planner always retains valid initial cuts. Zero workload also retains cuts.',
        requires_enabled_flag=True,
        explanation='Both controls need enabled/minimum_cells3 for128 ranks in942 cells; the parent8-layer-only disabled path cannot represent these thin slabs.')
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(build_record), accepted_for_execution=False, cases=cases,
        source_inputs={str(source):source_sha}, frozen_cards=frozen_cards,
        nodes=1, ranks=ranks, threads=1, recommended_slurm_minutes=90 if repeats == 2 else 45,
        repeats=repeats, field_scheduling_cost=field_cost,
        execution_order=[case['name'] for case in cases], profile_interval=0,
        physical_acceptance=False, submitted=False, reference_query_validation=query,
        static_control_proof=static_proof,
        common_baseline_changes=['active_nz342 to942 identically in every control',
            'reference_nz1024 to4096 identically in every control',
            'current matched three-point shape identically in every control',
            'particle output reduced to entry, und-entry and p07 identically in every control'],
        unchanged=['Source events/IDs/weights/charge', 'Reference center and automatic source clock',
            'Spacing/dt/boost', 'Seven-period magnet and original stop379 steps',
            'Transverse CPML and fast reference/inlet', 'Field and energy planes',
            'Conservative particle coupling and strict energy-ledger cadence'],
        scope='Paired short-source long-window performance stress, not37-period physical FEL acceptance; no gain is credited to changing the window.',
        next_step='Root must review reference bounds, candidate integrity, install/freeze a separate runner with profile0, and accept before submitting.')
    save(output / 'plan.json', plan)
    save(output / 'FROZEN_CARDS.json', frozen_cards)
    save(output / 'REFERENCE_QUERY_CHECK.json', query)
    save(output / 'HELD_BETA_STATISTICS.json', beta)
    save(output / 'STATIC_CONTROL_PROOF.json', static_proof)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
        accepted_for_execution=False, cases=len(cases), nodes=1, ranks=ranks, threads=1,
        repeats=repeats, field_scheduling_cost=field_cost,
        reference_query_check_passed=True, static_control_proven=True,
        source_sha256=source_sha, source_unchanged=True,
        normalized_physical_cards_identical=True))
    order_note = ('Four independent output directories run in static/balanced/balanced/static order. '
                  if repeats == 2 else 'Two independent output directories run in static/balanced order. ')
    (output / 'README.md').write_text(
        '# MPI balancing performance stress test\n\n'
        'Prepared only; no simulation has been run or submitted by this script.\n\n'
        'The unchanged1,429,799-particle w8 source traverses its original seven-period magnet. '
        'Both controls use the same942-cell active window and4096-cell reference, with the '
        'original input center, spacing, dt and379-step stop. This short-source, long-window '
        'configuration tests concentrated particle load. It is not the37-period FEL physical '
        'acceptance, and changing window size is not counted as a balancing speedup.\n\n'
        f'{order_note}'
        'Both controls enable thin-slab support; static hysteresis is above the universal '
        'possible gain bound, so its initial uniform cuts never change. The two controls '
        'differ only in that hysteresis. See STATIC_CONTROL_PROOF.json.\n\n'
        'The minimum slab width is3 cells. Field-cell scheduling weight is '
        f'{field_cell_weight:g} (effective {field_cost["effective_field_cell_weight"]:g}). '
        'Value0 means automatic Nx*Ny/64; positive values are scheduling costs only '
        'and never modify physical equations.\n\n'
        'C++ reference-query bounds are checked through ceil(duration/dt)+1 steps, '
        'including conservative held-beta limits. See REFERENCE_QUERY_CHECK.json. '
        'All cases retain strict ledgers, entry/exit field and energy planes, and three '
        'particle planes. Profile collection must remain disabled(profile_interval0). '
        f'Each case has a1200-second limit; the proposed whole allocation is{90 if repeats == 2 else 45}minutes.\n')
    print(json.dumps(dict(complete=True, output=str(output), cases=len(cases), repeats=repeats,
                          ranks=ranks, profile_interval=0, reference_query=query), indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--out', '--output', dest='output', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/performance-r1')
    parser.add_argument('--binary', type=Path,
                        default=WORKSPACE / 'build/mpi-particle-balance-20261106/aprl')
    parser.add_argument('--build-record', '--buildrecord', dest='build_record', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/build-r1')
    parser.add_argument('--ranks', type=int, default=128)
    parser.add_argument('--field-cell-weight', type=float, default=0.)
    parser.add_argument('--repeats', type=int, choices=[1, 2], default=2)
    args = parser.parse_args()
    main(args.workspace, args.output, args.binary, args.build_record, args.ranks,
         args.field_cell_weight, args.repeats)
