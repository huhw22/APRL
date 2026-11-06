"""Prepare the complete current 37-period, physical-electron FEL case.

This utility only reads a closed physical source and writes a new card/plan.
It never launches a simulation, copies/reweights particles, or creates the
reference directory (the blocked-reference writer must create it itself).
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np
import yaml
from scipy.constants import c, e


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
SOURCE = WORKSPACE / 'results/physical_electron_source_20261106/source-37periods-r1'
BASE = WORKSPACE / ('results/energy_critical_control_20261105/bf_origin_review/'
                    'loading_correlation_test_20261105/campaign-r2/original/run.yaml')
EXPECTED_SOURCE_SHA = '16be0bb3aff3d0d31e139a8d732db32fe237063fdd1d154fce2b17d9de49f125'
BLOCK = 131072


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(val) for key, val in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(val) for val in value]
    return value


def save(path, value):
    Path(path).write_text(json.dumps(native(value), indent=2, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def stat(path):
    value = Path(path).stat()
    return dict(bytes=value.st_size, mtime_ns=value.st_mtime_ns, inode=value.st_ino)


def chunks(data):
    for start in range(0, len(data), BLOCK):
        yield start, data[start:min(start + BLOCK, len(data))]


def velocities(records, boost_gamma):
    proper = records['proper_velocity']
    gamma = np.sqrt(1. + np.sum(proper * proper, axis=1))
    beta = proper[:, 2] / gamma
    light = ((np.hypot(1., proper[:, 2]) + proper[:, 2]) /
             (boost_gamma + np.sqrt(boost_gamma - 1) * np.sqrt(boost_gamma + 1)))
    held_u = .5 * (light - 1. / light)
    held_beta = held_u / np.hypot(1., held_u)
    return gamma, beta, held_beta


def source_statistics(path, boost_gamma, dz):
    """Predict the existing native loader's automatic frame in bounded memory.

    The three passes reproduce the accepted prepare expressions: source moments,
    earliest held arrival, then boosted initial extrema. Input events are never
    shifted or regenerated; only the common coordinate/time convention is set.
    """
    total_beta = np.longdouble(0)
    total_beta_time = np.longdouble(0)
    total_gamma = np.longdouble(0)
    tmin, tmax = math.inf, -math.inf
    held_min, held_max = math.inf, -math.inf
    boost_beta = math.sqrt(1. - 1. / boost_gamma**2)
    with h5py.File(path, 'r') as handle:
        assert int(handle['generation/complete'][()]) == 1
        data = handle['particles/records']
        counts = handle['generation/bin_counts'][:].astype(np.int64)
        count = len(data)
        assert int(counts.sum()) == count
        for start, records in chunks(data):
            assert np.all(records['macro_weight'] == 1.)
            assert np.array_equal(records['source_id'],
                                  np.arange(start + 1, start + len(records) + 1, dtype=np.uint64))
            gamma, beta, held = velocities(records, boost_gamma)
            time = records['arrival_time_offset_s'].astype(np.longdouble)
            total_beta += np.sum(beta.astype(np.longdouble))
            total_beta_time += np.sum(beta.astype(np.longdouble) * time)
            total_gamma += np.sum(gamma.astype(np.longdouble))
            tmin = min(tmin, float(time.min()))
            tmax = max(tmax, float(time.max()))
            held_min = min(held_min, float(held.min()))
            held_max = max(held_max, float(held.max()))
        mean_beta = total_beta / count
        mean_beta_time = total_beta_time / count
        center = float(c * (mean_beta * np.longdouble(tmax) - mean_beta_time) + 2 * dz / boost_gamma)
        reference_time = float((np.longdouble(center) / c + mean_beta_time) / mean_beta)
        origin = math.inf
        for _, records in chunks(data):
            _, beta, held = velocities(records, boost_gamma)
            snapshot_z = beta * c * (reference_time - records['arrival_time_offset_s'])
            input_times = -snapshot_z / (c * beta)
            origin = min(origin, float(np.min(input_times - 2 * dz /
                                              (boost_gamma * c * (boost_beta + held)))))
        zmin, zmax = math.inf, -math.inf
        for _, records in chunks(data):
            _, beta, held = velocities(records, boost_gamma)
            snapshot_z = beta * c * (reference_time - records['arrival_time_offset_s'])
            input_times = -snapshot_z / (c * beta)
            z = -boost_gamma * c * (boost_beta + held) * (input_times - origin)
            zmin = min(zmin, float(z.min()))
            zmax = max(zmax, float(z.max()))
    assert reference_time >= tmax
    assert origin < tmin - reference_time
    return dict(records=count, represented_electrons=count,
                mean_beta=float(mean_beta), mean_beta_time_s=float(mean_beta_time),
                mean_lab_gamma=float(total_gamma / count),
                center_lab_z_m=center, reference_time_offset_lab_s=reference_time,
                predicted_frame_time_origin_lab_s=origin,
                predicted_box_translation_m=-.5 * (zmin + zmax),
                held_beta_min=held_min, held_beta_max=held_max,
                arrival_time_min_s=tmin, arrival_time_max_s=tmax), counts


def prepare(output, binary, build_record, acceptance_record=None, nodes=2, case_limit_seconds=259200):
    output = Path(output).resolve()
    assert nodes >= 1
    assert not (output / 'plan.json').exists(), 'Refusing to replace an existing plan'
    assert not (output / 'one4one').exists(), 'Refusing to replace an existing case'
    source = SOURCE / 'particles.h5'
    before = stat(source)
    complete = json.loads((SOURCE / 'COMPLETE.json').read_text())
    source_cfg = json.loads((SOURCE / 'config.json').read_text())
    assert complete['complete'] and before['bytes'] == complete['files']['particles.h5']['bytes']
    assert complete['files']['particles.h5']['sha256'] == EXPECTED_SOURCE_SHA
    assert sha(source) == EXPECTED_SOURCE_SHA, 'Physical source differs from accepted source'
    base = yaml.safe_load(BASE.read_text())
    gb = float(base['mesh']['boost_gamma'])
    bb = math.sqrt(1. - 1. / gb**2)
    dz = float(base['mesh']['cell_size'][2])
    dt = dz / c
    info, counts = source_statistics(source, gb, dz)
    assert info['records'] == 60360290 == source_cfg['electrons']
    assert len(counts) == 111
    wl = float(source_cfg['wavelength_m'])
    periods = 37
    magnet_entrance = .01
    period = .028
    exit_z = magnet_entrance + period * periods
    span = (info['arrival_time_max_s'] - info['arrival_time_min_s'] + periods * wl / c +
            (1. - bb) * exit_z / c + 4e-15)
    steps = int(math.ceil(span * gb / dt))
    active = int(math.ceil(exit_z / (gb * dz))) + 64 + 128 + 2
    translation = info['predicted_box_translation_m']
    extent = max((1. + info['held_beta_max']) * (steps + 2) + 66 - translation / dz,
                 translation / dz - 64 + active + 2)
    reference_nz = int(math.ceil(2 * extent / 256)) * 256
    duration = (steps + .5) * dt
    stop = (gb * gb * (info['center_lab_z_m'] + bb * c *
                       info['predicted_frame_time_origin_lab_s']) +
            math.sqrt(gb * gb - 1) * dz * steps)
    plane_periods = list(range(2, periods, 2)) + [periods]
    particle_planes = [dict(name='entry', z=0.), dict(name='und-entry', z=magnet_entrance)]
    particle_planes += [dict(name=f'p{p:02d}', z=magnet_entrance + p * period)
                        for p in plane_periods]
    assert len(particle_planes) == 21
    prefix = np.r_[0, np.cumsum(counts)]
    first_original = int(round(source_cfg['profile']['window_start_s'] /
                               source_cfg['profile']['bin_width_s']))
    assert first_original == 637
    core_first, core_stop = 37, 74
    core_id_first, core_id_last = int(prefix[core_first]) + 1, int(prefix[core_stop])
    ncore = core_id_last - core_id_first + 1
    geometry = dict(nx=80, ny=80, reference_nz=reference_nz, active_nz=active,
                    dx_m=float(base['mesh']['cell_size'][0]),
                    dy_m=float(base['mesh']['cell_size'][1]), dz_m=dz,
                    dt_box_s=dt, boost_gamma=gb, boost_beta=bb, steps=steps,
                    duration_box_s=duration, first_magnet_z_m=magnet_entrance,
                    last_magnet_z_m=exit_z, periods=periods, source_slices=111,
                    source_nominal_length_m=111 * wl,
                    source_nominal_duration_s=111 * wl / c,
                    nominal_resonant_slippage_m=periods * wl,
                    source_length_over_resonant_slippage=3.,
                    nominal_core_length_m=37 * wl,
                    nominal_guard_wavelengths_each_side=37,
                    scope='Complete current 37-period benchmark, not original six by 77-period lattice.',
                    finite_random_support_note='Random source extrema lie inside nominal support; threefold ratio refers to declared support.',
                    **{key: info[key] for key in ['center_lab_z_m', 'reference_time_offset_lab_s',
                       'predicted_frame_time_origin_lab_s', 'predicted_box_translation_m',
                       'held_beta_min', 'held_beta_max']})
    cohort = dict(definition='Central37 of the111 original master generation bins; unchanged source IDs at every plane.',
                  parent_source=str(source), parent_crop_original_row_start=0,
                  parent_crop_original_row_stop=info['records'],
                  original_optical_start=first_original,
                  original_optical_stop=first_original + 111,
                  core_optical_start=first_original + core_first,
                  core_optical_stop=first_original + core_stop,
                  core_original_source_id_first=core_id_first,
                  core_original_source_id_last=core_id_last,
                  source_plane_z_m=0., parameter_plane_z_m=.01,
                  physical_electrons=info['records'], core_physical_electrons=ncore,
                  front_guard_optical_periods=37, rear_guard_optical_periods=37,
                  temporal_crop_preserves_original_event_values=True,
                  restriction='Margins cover nominal resonant slippage, not arbitrary broadband/quasistatic causal interaction.',
                  particle_planes=particle_planes)
    output.mkdir(parents=True, exist_ok=True)
    case = output / 'one4one'
    case.mkdir()
    for folder in ['detectors', 'ledger']:
        (case / folder).mkdir()
    assert not (case / 'reference').exists()
    cfg = copy.deepcopy(base)
    cfg['output']['manifest'] = str(case / 'run-manifest.yaml')
    cfg['mesh'].update(cells=[80, 80, reference_nz], duration=duration,
                       particle_coupling='conservative-source', longitudinal_shape='linear-binomial3')
    cfg['mesh']['moving_window'].update(cells=active, reference_directory=str(case / 'reference'))
    cfg['beam']['reference']['initial_center_z'] = info['center_lab_z_m']
    cfg['beam']['input'].update(file=str(source), electrons=float(info['records']))
    cfg['sources']['magnetic_elements'][0]['periods'] = periods
    cfg['detectors'].update(directory=str(case / 'detectors'), substep_particle_planes=True,
        particle_planes=[dict(name=plane['name'] + '-particles', z=plane['z'],
                              buffer_records=32768, compression=0) for plane in particle_planes],
        field_planes=[dict(name=name, z=z, rhythm=dt / gb, buffer_samples=8,
                           compression=1, particle_background={'enabled': False})
                      for name, z in [('entry', 0.), ('p37', exit_z)]],
        energy_planes=[dict(name=name + '-transport', z=z, rhythm=dt / gb,
                            buffer_samples=64, compression=0)
                       for name, z in [('entry', 0.), ('p37', exit_z)]])
    cfg['energy_ledger'].update(directory=str(case / 'ledger'), enabled=True,
                                discrete_operator_audit=True)
    cfg['stop'] = dict(mode='reference-center-z', z=float(stop))
    # Probe families used by old quiet-start tests do not apply to independent IDs.
    if 'force_probe' in cfg:
        cfg['force_probe']['enabled'] = False
    (case / 'run.yaml').write_text(yaml.safe_dump(native(cfg), sort_keys=False))
    source_meta = dict(path=str(source), records=info['records'], actual_macro_weight=1.,
                       physical_electrons=info['records'], charge_C=info['records'] * e,
                       file_sha256=EXPECTED_SOURCE_SHA, file_stat=before,
                       original_generation_bin_counts=counts.tolist(),
                       core_records=ncore, core_represented_electrons=float(ncore),
                       expected_full_B_rms=1. / math.sqrt(info['records']),
                       expected_core_B_rms=1. / math.sqrt(ncore),
                       input_selection='All original physical master records unchanged; no thinning/copying/reweighting.')
    plan_case = dict(name='one4one', mode='on', stage='full', directory=str(case),
        binary=str(Path(binary).resolve()), build_record=str(Path(build_record).resolve()),
        acceptance_record=str(Path(acceptance_record).resolve()) if acceptance_record else None,
        nodes=nodes, ranks=16 * nodes, threads=8, input=str(source),
        input_sha256=EXPECTED_SOURCE_SHA, expected_particles=info['records'],
        physical_electrons=info['records'], target_macro_weight=1, actual_macro_weight=1.,
        target_weight=1, actual_weight=1., selection_seed=source_cfg['seed'],
        source_selection_seed=source_cfg['seed'], expected_steps=steps,
        configured_duration_step_upper=math.ceil(duration / dt),
        predicted_stop_step=float(steps), stop_reference_center_z_m=float(stop),
        case_limit_seconds=case_limit_seconds, require_all_particles_through_planes=True,
        wavelength_m=wl, first_original_slice=first_original,
        end_original_slice_exclusive=first_original + 111,
        core_first_slice_inclusive=first_original + core_first,
        core_end_slice_exclusive=first_original + core_stop,
        cohort=str(output / 'COHORT.json'), core_original_source_id_first=core_id_first,
        core_original_source_id_last=core_id_last, expected_core_records=ncore,
        expected_core_represented_electrons=float(ncore), force_probe={'enabled': False},
        source_metadata=source_meta, predicted_frame=str(case / 'PREDICTED_FRAME.json'),
        particle_planes=particle_planes, planes=particle_planes,
        card_sha256=sha(case / 'run.yaml'),
        recommended_slurm_minutes=math.ceil((case_limit_seconds + 3600) / 60))
    output_cost = dict(legacy_particle_bytes_per_record=100, substep_particle_bytes_per_record=164,
                       particle_planes=len(particle_planes), particles_per_plane=info['records'],
                       particle_payload_bytes=info['records'] * 264 * len(particle_planes),
                       particle_payload_GiB=info['records'] * 264 * len(particle_planes) / 2**30,
                       field_payload_upper_GiB=2 * (steps + 1) * 80 * 80 * 48 / 2**30,
                       scope='Uncompressed dataset payload only; HDF5 metadata, reference files, and ledger are additional.')
    save(output / 'GEOMETRY.json', geometry)
    save(output / 'COHORT.json', cohort)
    save(output / 'SOURCES.json', {'one4one': source_meta})
    save(output / 'OUTPUT_COST.json', output_cost)
    save(case / 'PREDICTED_FRAME.json', info)
    save(case / 'run-plan.json', plan_case)
    assert stat(source) == before, 'Source stat changed during read-only preparation'
    save(output / 'SOURCE_PROVENANCE.json', dict(parent=str(source), parent_stat_before=before,
         parent_stat_after=stat(source), parent_sha256=EXPECTED_SOURCE_SHA,
         verification_scope='Full file SHA repeated; closed completion/count/weight1/ID checks and three bounded geometry passes.',
         source_copied_or_changed=False, parameter_plane_events_changed=False))
    plan = dict(workspace=str(WORKSPACE), directory=str(output), binary=plan_case['binary'],
                build_record=plan_case['build_record'], acceptance_record=plan_case['acceptance_record'],
                cases=[plan_case], geometry=geometry,
                cohort=str(output / 'COHORT.json'), accepted_for_execution=False,
                pending='Root must verify accepted binary, user scope choice, frozen scripts/cards and resource layout before submission.',
                purpose='Current complete37-period FEL with physical independent electrons and matched three-point energy-conservative coupling.',
                unchanged_methods=['Cowan-z', 'conservative-source', 'fast reference',
                                   'inlet matching', 'transverse CPML', 'MPI shutdown'],
                controlled_changes=['Matched linear-binomial3 longitudinal shape',
                                    'All111 source bins instead of7-period21-bin crop',
                                    '37-period propagation and two-period BF sampling'],
                force_probe_enabled=False, submitted=False)
    save(output / 'plan.json', plan)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
         accepted_for_execution=False, source_master_unchanged=True, geometry=geometry,
         output_cost=output_cost, card_sha256=plan_case['card_sha256']))
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=WORKSPACE / 'results/binomial3_fel_20261106/full-r1')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-record', type=Path, required=True)
    parser.add_argument('--acceptance-record', type=Path)
    parser.add_argument('--nodes', type=int, default=2)
    parser.add_argument('--case-limit-seconds', type=int, default=259200)
    args = parser.parse_args()
    plan = prepare(args.output, args.binary, args.build_record, args.acceptance_record,
                   args.nodes, args.case_limit_seconds)
    print(json.dumps(plan, indent=2))


if __name__ == '__main__':
    main()
