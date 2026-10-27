#!/usr/bin/env python3
"""Stream the full held input inside Slurm and predict the compact field-only mesh.

This mirrors initial kinematics/layout, not evolved magnetic/self-field motion.
The after-last-element stop remains dynamic and independent of particle output.
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
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(16 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def kinematics(rows, gamma_b, beta_b):
    u = rows['proper_velocity']
    gamma = np.sqrt(1 + np.sum(u * u, axis=1))
    transverse_mass2 = 1 + u[:, 0] ** 2 + u[:, 1] ** 2
    rapidity = gamma_b * (1 + beta_b)
    plus = (gamma + u[:, 2]) / rapidity
    minus = transverse_mass2 / plus
    beta_box = (plus - minus) / (plus + minus)
    held_plus = (np.hypot(1., u[:, 2]) + u[:, 2]) / rapidity
    beta_held = (held_plus ** 2 - 1) / (held_plus ** 2 + 1)
    return gamma, u[:, 2] / gamma, beta_box, beta_held


def analyze(root):
    require(os.environ.get('SLURM_JOB_ID'), 'Full input geometry must run in a Slurm allocation')
    root = root.resolve()
    output = root / 'geometry-preflight.json'
    require(not output.exists(), 'Refusing to overwrite a geometry result')
    plan = json.loads((root / 'run-plan.json').read_text())
    history = json.loads(Path(plan['historical_plan']).read_text())
    card = yaml.safe_load((root / 'run.yaml').read_text())
    variant = plan.get('case_variant', 'undulator')
    require(variant in ('undulator', 'free-drift'), 'Unknown field-only case variant')
    if variant == 'free-drift':
        require(card['sources'].get('magnetic_elements') == [] and
                card['sources'].get('incident_waves') == [], 'Free drift must have no prescribed sources')
        require(plan['field_plane_z_m'] == .5, 'Free-drift field plane must be at 0.5 m')
        require(not card.get('checkpoint') and not plan.get('checkpoint_steps', []),
                'Free-drift comparison does not request a checkpoint')
    mesh, ref = card['mesh'], card['beam']['reference']
    require(card['units'] == {'length': 'm', 'time': 's'}, 'Planner expects SI card units')
    require(mesh['field_solver'] == 'cowan-z' and mesh['cells'][2] == 'auto', 'Expected automatic Cowan grid')
    require(mesh['longitudinal']['mode'] == 'observation-bounded', 'Expected observation-bounded layout')
    require(mesh['longitudinal']['self_field_padding_rms'] == 0, 'Planner supports original zero RMS padding only')
    require(ref['frozen_lab_plane_injection'] and ref['hold_transverse_until_activation'], 'Expected held startup')
    require(ref['activation_plane_z'] == ref['input_plane_z'], 'Activation must coincide with input plane')
    require(card['stop']['mode'] == 'after-last-element', 'Expected dynamic field-plane stop')
    require(not card.get('particle_retirement', {}).get('enabled', False), 'Retirement must be disabled')
    require(not card['detectors'].get('particle_planes', []), 'Particle-plane exports must be absent')
    planes = card['detectors']['field_planes']
    require(len(planes) == 1 and not planes[0].get('particle_background', {}).get('enabled', False), 'Expected one raw field-only detector')
    require(not mesh.get('ballistic_outside_transverse', False), 'Held case uses ordinary transverse boundary handling')
    require(card['beam']['input'].get('position_offset', [0., 0., 0.]) == [0., 0., 0.], 'Unexpected beam offset')
    require(card['boundary']['type'] == 'cpml', 'Expected CPML')

    path = Path(plan['input_file']).resolve()
    before = path.stat()
    require((before.st_size, before.st_mtime_ns) == (plan['input_size_bytes'], plan['input_mtime_ns']), 'Input stat differs from frozen plan')
    require(str(Path(card['beam']['input']['file']).resolve()) == str(path), 'Card input differs from plan')
    sha = digest(path)
    require(sha == plan['input_sha256'], 'Full input checksum differs from frozen plan')

    gamma_b = mesh['boost_gamma']
    beta_b = math.sqrt(gamma_b - 1) * math.sqrt(gamma_b + 1) / gamma_b
    dx, dy, dz = mesh['cell_size']
    nx, ny = mesh['cells'][:2]
    cpml_x, cpml_y, cpml_z = card['boundary']['cells']
    half = np.array([nx * dx / 2 - cpml_x * dx, ny * dy / 2 - cpml_y * dy])
    center_xy = np.array(mesh['center'][:2])
    require(np.all(half > 0), 'No physical transverse aperture')
    z_input, z_reference = ref['input_plane_z'], ref['initial_center_z']
    require(z_input == history['input_plane_z_m'] and z_reference == history['initial_center_z_m'], 'Beam reference differs from historical input plan')
    common = history['common_time_after_reference_s']
    field_z = planes[0]['z']
    require(field_z == plan['field_plane_z_m'] and field_z > z_input, 'Invalid downstream field plane')
    for magnet in card['sources']['magnetic_elements']:
        require(magnet['type'] == 'planar-undulator', 'Unexpected magnet in this case')
        # The production parser is authoritative for the exact fringe extent.
        physical_exit = magnet['entrance_z'] + magnet['period'] * magnet['periods']
        require(field_z > physical_exit + 5 * magnet['period'], 'Field plane not safely beyond the declared undulator fringe')

    count_seen, outside, nonunit = 0, 0, 0
    weight = np.longdouble(0)
    weighted_beta = np.longdouble(0)
    weighted_beta_time = np.longdouble(0)
    t0 = math.inf
    gamma_min, gamma_max = math.inf, -math.inf
    with h5py.File(path, 'r') as file:
        group = file['particles']; data = group['records']; count = len(data)
        require(int(group.attrs['format_version']) == 4, 'Expected HDF5 v4 input')
        require(count == plan['records'] == history['records'], 'Particle count differs from frozen plans')
        required = ['plane_position_m', 'arrival_time_offset_s', 'proper_velocity', 'macro_weight']
        require(set(required).issubset(data.dtype.names), 'Unexpected particle schema')
        print(f'geometry pass 1/3: input validation, held origin, {count} records', flush=True)
        for begin in range(0, count, CHUNK):
            rows = data[begin:begin + CHUNK]
            require(all(np.isfinite(rows[key]).all() for key in required), 'Nonfinite input particle value')
            require(np.all(rows['macro_weight'] > 0) and np.all(rows['proper_velocity'][:, 2] > 0), 'Invalid weight or longitudinal velocity')
            gamma, beta_lab, _, beta_held = kinematics(rows, gamma_b, beta_b)
            alpha = beta_b + beta_held
            require(np.all(alpha > 0), 'Unreachable held longitudinal motion')
            arrival = rows['arrival_time_offset_s'] - common
            t0 = min(t0, float(np.min(arrival - 2 * dz / (gamma_b * C * alpha))))
            gamma_min = min(gamma_min, float(np.min(gamma)))
            gamma_max = max(gamma_max, float(np.max(gamma)))
            outside += int(np.count_nonzero(np.any(np.abs(rows['plane_position_m'] - center_xy) > half, axis=1)))
            nonunit += int(np.count_nonzero(rows['macro_weight'] != 1))
            weights = rows['macro_weight'].astype(np.longdouble)
            wide_beta = beta_lab.astype(np.longdouble)
            weight += np.sum(weights, dtype=np.longdouble)
            weighted_beta += np.sum(weights * wide_beta, dtype=np.longdouble)
            weighted_beta_time += np.sum(weights * wide_beta * rows['arrival_time_offset_s'], dtype=np.longdouble)
            count_seen += len(rows)
        recomputed_common = float(((np.longdouble(z_reference) - np.longdouble(z_input)) / C * weight + weighted_beta_time) / weighted_beta)
        require(abs(recomputed_common - common) < 1e-20, 'Historical common-time reconstruction mismatch')
        require(count_seen == count and outside == 0 and nonunit == 0, 'Expected full unit-weight beam inside held physical aperture')
        require(float(weight) == float(plan['electrons']), 'Represented particle weight mismatch')

        initial_min, initial_max = math.inf, -math.inf
        ballistic_stop = 0.
        print('geometry pass 2/3: initial support and initial-ballistic stop estimate', flush=True)
        for begin in range(0, count, CHUNK):
            rows = data[begin:begin + CHUNK]
            _, _, beta_box, beta_held = kinematics(rows, gamma_b, beta_b)
            activation = gamma_b * (rows['arrival_time_offset_s'] - common - t0)
            initial = -C * (beta_b + beta_held) * activation
            initial_min = min(initial_min, float(np.min(initial)))
            initial_max = max(initial_max, float(np.max(initial)))
            lab_advance = gamma_b * C * (beta_b + beta_box)
            require(np.all(lab_advance > 0), 'Initial ballistic stop unreachable')
            remaining = field_z - (z_input + gamma_b * initial)
            ballistic_stop = max(ballistic_stop, float(np.max(np.maximum(0., remaining) / lab_advance)))
        planned = ballistic_stop * mesh['longitudinal']['ballistic_time_safety_factor']
        ballistic_min, ballistic_max = initial_min, initial_max
        worldtube_min, worldtube_max = initial_min, initial_max
        print('geometry pass 3/3: independent compact layout and held worldtube envelope', flush=True)
        for begin in range(0, count, CHUNK):
            rows = data[begin:begin + CHUNK]
            _, _, beta_box, beta_held = kinematics(rows, gamma_b, beta_b)
            activation = gamma_b * (rows['arrival_time_offset_s'] - common - t0)
            initial = -C * (beta_b + beta_held) * activation
            finish = initial + beta_box * C * planned
            ballistic_min = min(ballistic_min, float(np.min(finish)))
            ballistic_max = max(ballistic_max, float(np.max(finish)))
            switch_time = np.minimum(activation, planned)
            switched = initial + beta_held * C * switch_time
            held_finish = switched + beta_box * C * np.maximum(0., planned - activation)
            worldtube_min = min(worldtube_min, float(np.min(switched)), float(np.min(held_finish)))
            worldtube_max = max(worldtube_max, float(np.max(switched)), float(np.max(held_finish)))

    head_lab = z_input + gamma_b * initial_max
    radiation_advance = (field_z - head_lab) / (gamma_b * (1 + beta_b))
    radiation_upper = initial_max + radiation_advance
    lower = min(initial_min, ballistic_min)
    upper = max(initial_max, ballistic_max, radiation_upper)
    guards = mesh['longitudinal']
    required_lower = lower - guards['lower_guard_cells'] * dz
    required_upper = upper + guards['upper_guard_cells'] * dz
    minimum_physical_cells = math.ceil((np.longdouble(required_upper) - np.longdouble(required_lower)) / np.longdouble(dz))
    alignment = guards['alignment_cells']
    nz = math.ceil((minimum_physical_cells + 2 * cpml_z) / alignment) * alignment
    extent = nz * dz
    physical_lower = mesh['center'][2] - .5 * extent + cpml_z * dz
    physical_upper = mesh['center'][2] + .5 * extent - cpml_z * dz
    translation = physical_lower - required_lower
    dt = dz / C
    estimated_steps = math.ceil(ballistic_stop / dt)
    first_field_box_z = (field_z - z_input) / gamma_b + translation
    first_lab_time = t0 + beta_b * (field_z - z_input) / C
    estimated_last_lab_time = first_lab_time + estimated_steps * dt / gamma_b
    predicted_plane_end = first_field_box_z - beta_b * C * planned
    after = path.stat()
    checks = {
        'input_checksum_verified': sha == plan['input_sha256'],
        'input_unchanged': (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
        'all_records_finite_positive_unit_weight': count_seen == count and nonunit == 0,
        'all_initial_particles_in_physical_xy': outside == 0,
        'held_worldtube_in_physical_z_for_planned_time': worldtube_min + translation >= physical_lower and worldtube_max + translation <= physical_upper,
        'detector_in_physical_z_for_planned_time': first_field_box_z <= physical_upper and predicted_plane_end >= physical_lower,
        'hard_duration_above_planned_time': mesh['duration'] > planned,
        'minimum_two_z_cells_per_rank': nz // plan['ranks'] >= 2,
        'checkpoint_step_inside_estimated_run': all(0 < n < estimated_steps for n in plan.get('checkpoint_steps', [])),
    }
    result = {
        'format': 'aprl-held-field-only-geometry-v1', 'case_variant': variant,
        'pass': all(checks.values()), 'checks': checks,
        'input_file': str(path), 'input_sha256': sha, 'records': count, 'input_weight_sum': float(weight),
        'frame_origin_lab_time_s': t0, 'frame_origin_lab_z_m': z_input,
        'recomputed_common_time_s': recomputed_common,
        'estimated_cells_z': nz, 'estimated_ballistic_steps': estimated_steps,
        'estimated_ballistic_stop_box_s': ballistic_stop, 'planned_ballistic_time_s': planned,
        'dt_box_s': dt, 'detector_z_m': field_z, 'gamma_boost': gamma_b, 'beta_boost': beta_b,
        'input_gamma_range': [gamma_min, gamma_max],
        'physical_bounds_z_m': [physical_lower, physical_upper], 'translation_z_m': translation,
        'initial_bounds_untranslated_m': [initial_min, initial_max],
        'ballistic_bounds_untranslated_m': [ballistic_min, ballistic_max],
        'held_worldtube_bounds_untranslated_m': [worldtube_min, worldtube_max],
        'radiation_upper_untranslated_m': radiation_upper, 'minimum_physical_cells': minimum_physical_cells,
        'detector_window': {'first_lab_time_s': first_lab_time,
            'estimated_ballistic_last_lab_time_s': estimated_last_lab_time,
            'first_box_z_m': first_field_box_z,
            'planned_ballistic_last_box_z_m': predicted_plane_end,
            'first_upper_margin_cells': (physical_upper - first_field_box_z) / dz,
            'planned_last_lower_margin_cells': (predicted_plane_end - physical_lower) / dz},
        'assumptions': [
            'Exact same initial kinematic formulas and compact guards as the simulation; Python/C++ floating-point summation may differ slightly.',
            'The final step is selected by surviving evolved particles crossing the field location; no particle-plane file is needed.',
            'Estimated ballistic steps are not a prescribed final step or a proof of quiet field tail.',
            'No added100fs quiet tail,20fs dynamics allowance,retirement distance or arbitrary enlarged guard.',
            'Held origin depends on input and activation plane only; changing field plane/mesh translation does not change boost time.',
        ],
    }
    if plan.get('checkpoint_steps'):
        result['assumptions'].append('Checkpoint leading-particle position is estimated from prior held trajectories and must be checked from its saved state.')
    if variant == 'free-drift':
        result['assumptions'].append('No prescribed magnetic or incident field is present. Held startup and self-consistent particle motion prevent a claim of exactly zero radiation.')
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    require(result['pass'], 'Geometry preflight failed; inspect checks in ' + str(output))
    print(f'PASS geometry: Nz={nz}, initial-ballistic steps={estimated_steps}, final stop remains dynamic', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    analyze(parser.parse_args().root)
