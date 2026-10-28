#!/usr/bin/env python3
"""Derive a compact moving field window and its finite scalar source domain.

Metadata planning is light. --verify-input streams the original particle input
only inside a Slurm allocation. It never runs APRL or modifies an existing card.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np
import yaml

C = 299792458.0
CHUNK = 262144


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def derive(source_run, ranks=16):
    source_run = Path(source_run).resolve()
    card_path = source_run / 'run.yaml'
    old_path = source_run / 'geometry-preflight.json'
    card = yaml.safe_load(card_path.read_text())
    old = json.loads(old_path.read_text())
    require(old['pass'], 'Historical full-input geometry must have passed')
    mesh, ref = card['mesh'], card['beam']['reference']
    require(card['units'] == {'length': 'm', 'time': 's'}, 'Expected SI units')
    require(ref['frozen_lab_plane_injection'] and ref['hold_transverse_until_activation'], 'Expected held release-event definition')
    require(ref['activation_plane_z'] == ref['input_plane_z'], 'Activation/input plane mismatch')
    require(not mesh['ballistic_outside_transverse'], 'Unexpected outside-transverse mode')
    require(card['stop']['mode'] == 'after-last-element', 'Preserve dynamic detector stop')
    require(len(card['detectors']['field_planes']) == 1, 'Expected one raw field plane')
    require(not card['detectors']['particle_planes'], 'Expected no particle-plane output')
    require(not card.get('particle_retirement', {}).get('enabled', False), 'No retirement layer requested')
    require(ranks > 0, 'Invalid MPI size')
    gamma = float(mesh['boost_gamma'])
    gb = math.sqrt(gamma - 1.) * math.sqrt(gamma + 1.)
    beta = gb / gamma
    dz = float(mesh['cell_size'][2]); dt = dz / C
    lower = int(mesh['longitudinal']['lower_guard_cells'])
    upper = int(mesh['longitudinal']['upper_guard_cells'])
    alignment = int(mesh['longitudinal']['alignment_cells'])
    require((lower, upper) == (64, 128), 'Preserve original 64/128 guards')
    require(alignment > 0 and alignment % 2 == 0, 'Require even source-grid alignment')
    safety = float(mesh['longitudinal']['ballistic_time_safety_factor'])
    duration = float(old['planned_ballistic_time_s'])
    require(math.isclose(duration, old['estimated_ballistic_stop_box_s'] * safety, rel_tol=1e-13), 'Historical duration basis differs')
    query_steps = math.ceil(duration / dt) + 1
    release = float(ref['activation_plane_z'])
    detector = float(card['detectors']['field_planes'][0]['z'])
    for magnet in card['sources']['magnetic_elements']:
        require(magnet['type'] == 'planar-undulator', 'Unsupported prescribed element')
        # The configured detector is well outside this declared Gaussian fringe.
        fringe = magnet['period'] / (2 * math.pi) * math.sqrt(-2 * math.log(magnet['fringe_relative_cutoff']))
        require(detector > magnet['entrance_z'] + magnet['period'] * magnet['periods'] + fringe,
                'Detector is not the final interaction boundary')
    active = math.ceil((detector - release) / (gamma * dz)) + lower + upper + 2
    require(active >= 8 * ranks, 'MPI size would require widening the active field window')
    zmin, zmax = old['initial_bounds_untranslated_m']
    center = float(mesh['center'][2])
    translation = center - .5 * (zmin + zmax)
    # For an explicit mesh Simulation centres the complete source support.
    # All reference-index expressions depend on translation minus mesh centre.
    translated_release_cells = (translation - center) / dz
    release_cell_floor = math.floor(translated_release_cells)
    gamma_max = float(old['input_gamma_range'][1])
    axial_u_upper = math.sqrt((gamma_max - 1.) * (gamma_max + 1.))
    light_front = (gamma_max + axial_u_upper) / (gamma + gb)
    held_beta_upper = (light_front * light_front - 1.) / (light_front * light_front + 1.)
    # An outward rounding allowance, not an extra time/window margin.
    held_beta_upper += 1e-12
    require(0 < held_beta_upper < 1, 'Invalid source held-velocity bound')
    required_half_lower = (1 + held_beta_upper) * query_steps + 1 - release_cell_floor + lower
    required_half_upper = release_cell_floor - lower + active + 1
    required_half_source = max(abs(zmin + translation - center), abs(zmax + translation - center)) / dz + 1
    nz = math.ceil(2 * max(required_half_lower, required_half_upper, required_half_source) / alignment) * alignment
    origin = center - .5 * nz * dz
    offset = math.floor((translation - origin) / dz) - lower
    minimum_query = offset - (1 + held_beta_upper) * query_steps - 1
    maximum_query = offset + active
    require(minimum_query >= 0 and maximum_query <= nz - 1, 'Finite reference does not cover all native queries')
    margin = 4 + math.ceil((1 - beta) * query_steps)
    window_origin = origin + offset * dz
    safe_exit = release + gamma * (window_origin + (active - margin) * dz - translation)
    require(safe_exit > detector + 2 * gamma * dz, 'Right causal collar overlaps detector')
    return {
        'format': 'aprl-moving-window-production-geometry-v1',
        'source_run': str(source_run), 'source_card_sha256': digest(card_path),
        'source_geometry_sha256': digest(old_path),
        'source_input_file': str(Path(old['input_file']).resolve()),
        'source_input_sha256': old['input_sha256'], 'source_records': old['records'],
        'reference_cells_z': nz, 'active_cells_z': active, 'mpi_ranks': ranks,
        'minimum_local_field_cells': active // ranks, 'initial_offset': offset,
        'duration_box_s': duration, 'reference_query_steps': query_steps,
        'minimum_query_bound_cells': minimum_query, 'maximum_query_bound_cells': maximum_query,
        'translation_z_m': translation, 'reference_origin_z_m': origin,
        'active_origin_z_m': window_origin, 'safe_exit_lab_z_m': safe_exit,
        'right_collar_cells': margin, 'bound_held_beta': held_beta_upper,
        'gamma_boost': gamma, 'beta_boost': beta, 'dt_box_s': dt, 'dz_m': dz,
        'source_initial_untranslated_bounds_m': [zmin, zmax],
        'source_frame_origin_lab_time_s': old['frame_origin_lab_time_s'],
        'source_common_time_s': old['recomputed_common_time_s'],
        'source_gamma_range': old['input_gamma_range'],
        'source_ballistic_stop_box_s': old['estimated_ballistic_stop_box_s'],
        'release_lab_z_m': release, 'detector_lab_z_m': detector,
        'source_alignment_cells': alignment, 'active_guard_cells': [lower, upper],
        'inherited_ballistic_time_safety_factor': safety,
        'notes': [
            'Source reference domain stores scalar blocks; it is not the resident E/B field grid.',
            'Duration reuses the old full-input ballistic stop times its original safety factor; dynamic after-last-element stop is unchanged.',
            'No added fs tail or duration padding. The hard bound is not proof that changed self-consistent trajectories always finish before it.',
            'Explicit source centring is a coordinate translation; input/release position and momentum remain the existing event definition.',
            'Held flags encode pending particles and release events; reference generation is the fast blocked Cowan solve, not a long held field pre-evolution.',
            'Upper held beta follows the original complete-input gamma maximum; startup must recheck the actual source and runtime query range.',
        ],
    }


def kinematics(rows, gamma, gb):
    u = rows['proper_velocity']
    gamma_lab = np.sqrt(1 + np.sum(u * u, axis=1))
    mass2 = 1 + u[:, 0] ** 2 + u[:, 1] ** 2
    plus = (gamma_lab + u[:, 2]) / (gamma + gb)
    minus = mass2 / plus
    beta_box = (plus - minus) / (plus + minus)
    axial_plus = (np.hypot(1., u[:, 2]) + u[:, 2]) / (gamma + gb)
    held = (axial_plus ** 2 - 1) / (axial_plus ** 2 + 1)
    return gamma_lab, u[:, 2] / gamma_lab, beta_box, held


def owners(z, origin, dz, cells, ranks, pending=False):
    coordinate = (z - origin) / dz
    require(np.all(coordinate <= cells + 1e-6), 'Source exceeds upper slab boundary')
    if not pending:
        require(np.all(coordinate >= -1e-6), 'Source exceeds lower slab boundary')
    cell = np.clip(np.floor(coordinate), 0, cells - 1).astype(np.int64)
    base, remainder = divmod(cells, ranks)
    long_cells = (base + 1) * remainder
    return np.where(cell < long_cells, cell // (base + 1),
                    remainder + (cell - long_cells) // base)


def check_cic_support(position, origin, spacing, cells, lower, upper, axis):
    """Mirror the blocked generator's nonzero one-dimensional CIC support."""
    coordinate = np.clip((position - origin) / spacing, 0., float(cells))
    cell = np.minimum(np.floor(coordinate), cells - 1).astype(np.int64)
    fraction = np.clip(coordinate - cell, 0., 1.)
    require(np.all((fraction == 1.) | ((cell > lower) & (cell < cells - upper))),
            'Reference CIC lower node touches Dirichlet surface: ' + axis)
    require(np.all((fraction == 0.) | ((cell + 1 > lower) & (cell + 1 < cells - upper))),
            'Reference CIC upper node touches Dirichlet surface: ' + axis)


def verify_full_input(plan):
    require(os.environ.get('SLURM_JOB_ID'), 'Full input scan requires Slurm')
    path = Path(plan['source_input_file']); before = path.stat()
    require(digest(path) == plan['source_input_sha256'], 'Particle input SHA mismatch')
    old_card = yaml.safe_load((Path(plan['source_run']) / 'run.yaml').read_text())
    gamma = plan['gamma_boost']; gb = math.sqrt(gamma - 1) * math.sqrt(gamma + 1)
    beta = gb / gamma; dz = plan['dz_m']; common = plan['source_common_time_s']
    mesh = old_card['mesh']; cpml = old_card['boundary']['cells']
    half = np.array([(mesh['cells'][i] * .5 - cpml[i]) * mesh['cell_size'][i] for i in (0, 1)])
    center = np.array(mesh['center'][:2]); t0 = math.inf
    min_h, max_h = math.inf, -math.inf; beta_sum = np.longdouble(0); count = 0
    with h5py.File(path, 'r') as f:
        group = f['particles']; data = group['records']
        require(int(group.attrs['format_version']) == 4 and len(data) == plan['source_records'], 'Input schema/count mismatch')
        print('Full input pass 1/2: coordinates, source velocity and release frame', flush=True)
        for start in range(0, len(data), CHUNK):
            rows = data[start:start + CHUNK]
            for key in ('plane_position_m', 'arrival_time_offset_s', 'proper_velocity', 'macro_weight'):
                require(np.isfinite(rows[key]).all(), 'Nonfinite particle value')
            require(np.all(rows['macro_weight'] == 1), 'Expected original unit weights')
            require(np.all(np.abs(rows['plane_position_m'] - center) <= half), 'Held source outside physical transverse aperture')
            for axis in (0, 1):
                spacing, cells = mesh['cell_size'][axis], mesh['cells'][axis]
                check_cic_support(rows['plane_position_m'][:, axis], center[axis] - .5 * cells * spacing,
                                  spacing, cells, cpml[axis], cpml[axis], 'xy'[axis])
            _, _, _, held = kinematics(rows, gamma, gb)
            require(np.all(held >= 0) and np.all(held < plan['bound_held_beta']), 'Held velocity violates reference bound')
            t0 = min(t0, float(np.min(rows['arrival_time_offset_s'] - common - 2 * dz / (gamma * C * (beta + held)))))
            min_h = min(min_h, float(np.min(held))); max_h = max(max_h, float(np.max(held)))
            beta_sum += np.sum(held.astype(np.longdouble)); count += len(rows)
        zmin, zmax, ballistic = math.inf, -math.inf, 0.
        reference_hist = np.zeros(plan['mpi_ranks'], dtype=np.int64)
        compact_hist = np.zeros(plan['mpi_ranks'], dtype=np.int64)
        pending_count = 0
        print('Full input pass 2/2: support, field passage bound and finite reference', flush=True)
        for start in range(0, len(data), CHUNK):
            rows = data[start:start + CHUNK]
            _, _, beta_box, held = kinematics(rows, gamma, gb)
            activation = gamma * (rows['arrival_time_offset_s'] - common - t0)
            initial = -C * (beta + held) * activation
            zmin = min(zmin, float(np.min(initial))); zmax = max(zmax, float(np.max(initial)))
            position = initial + plan['translation_z_m']
            check_cic_support(position, plan['reference_origin_z_m'], dz,
                              plan['reference_cells_z'], 0, 0, 'z')
            reference_hist += np.bincount(owners(position, plan['reference_origin_z_m'], dz,
                plan['reference_cells_z'], plan['mpi_ranks']), minlength=plan['mpi_ranks'])
            compact_hist += np.bincount(owners(position, plan['active_origin_z_m'], dz,
                plan['active_cells_z'], plan['mpi_ranks'], pending=True), minlength=plan['mpi_ranks'])
            pending_count += int(np.count_nonzero(position < plan['active_origin_z_m']))
            advance = gamma * C * (beta + beta_box)
            require(np.all(advance > 0), 'Unreachable initial ballistic detector event')
            remaining = plan['detector_lab_z_m'] - (plan['release_lab_z_m'] + gamma * initial)
            ballistic = max(ballistic, float(np.max(np.maximum(0., remaining) / advance)))
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'Input changed during scan')
    require(count == plan['source_records'], 'Incomplete input scan')
    require(abs(t0 - plan['source_frame_origin_lab_time_s']) < 1e-20, 'Release-frame origin changed')
    require(np.max(np.abs(np.array([zmin, zmax]) - plan['source_initial_untranslated_bounds_m'])) < dz * 1e-5,
            'Initial source support differs from frozen metadata')
    require(abs(ballistic - plan['source_ballistic_stop_box_s']) < plan['dt_box_s'] * 1e-5, 'Ballistic stop differs from frozen metadata')
    mean_beta = float(beta_sum / count)
    query_min = plan['initial_offset'] - (1 + mean_beta) * plan['reference_query_steps'] - 1
    require(query_min >= 0, 'Actual mean source velocity exceeds finite reference')
    # Initial pre-compaction redistribution still uses int byte displacements.
    # ParticlePacket is 104 bytes in the tested GNU/x86-64 executable. Compact
    # redistribution is already implemented with bounded 64 MiB packets.
    packet_bytes = 104
    initial_max_send_records = math.ceil(count / plan['mpi_ranks'])
    maximum_initial_send_bytes = initial_max_send_records * packet_bytes
    maximum_initial_receive_bytes = int(reference_hist.max()) * packet_bytes
    require(max(maximum_initial_send_bytes, maximum_initial_receive_bytes) <= 2**31 - 1,
            'Initial full-reference MPI_Alltoallv packet/displacement overflows int')
    return {'pass': True, 'records': count, 'input_sha256': plan['source_input_sha256'],
            'input_unchanged': True, 'frame_origin_lab_time_s': t0,
            'blocked_reference_cic_support_inside_dirichlet_faces': True,
            'initial_untranslated_bounds_m': [zmin, zmax],
            'held_beta_min_max': [min_h, max_h], 'mean_held_beta': mean_beta,
            'actual_mean_beta_query_min_cells': query_min,
            'initial_reference_owner_counts': reference_hist.tolist(),
            'initial_compact_owner_counts': compact_hist.tolist(),
            'initial_pending_upstream_records': pending_count,
            'legacy_initial_alltoallv': {'particle_packet_bytes': packet_bytes,
                'maximum_send_bytes': maximum_initial_send_bytes,
                'maximum_receive_bytes': maximum_initial_receive_bytes,
                'int_limit': 2**31 - 1, 'pass': True},
            'ballistic_stop_box_s': ballistic, 'duration_box_s': plan['duration_box_s'],
            'scope': 'Complete input identity and kinematic domain coverage; not nonlinear trajectory completion or radiation convergence.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--source-run', type=Path)
    mode.add_argument('--verify-plan', type=Path)
    p.add_argument('--ranks', type=int, default=16)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify-input', action='store_true')
    args = p.parse_args()
    require(not args.output.exists(), 'Use a new geometry output path')
    if args.verify_plan:
        plan = json.loads(args.verify_plan.read_text())
        require(derive(Path(plan['source_run']), plan['mpi_ranks']) == plan,
                'Frozen geometry plan no longer matches parent card and planning evidence')
        result = verify_full_input(plan)
    else:
        result = derive(args.source_run, args.ranks)
    if args.verify_input and not args.verify_plan:
        result['full_input_verification'] = verify_full_input(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    keys = ('pass', 'records', 'duration_box_s') if args.verify_plan else ('reference_cells_z', 'active_cells_z', 'initial_offset', 'duration_box_s', 'safe_exit_lab_z_m')
    print(json.dumps({k: result[k] for k in keys}))


if __name__ == '__main__':
    main()
