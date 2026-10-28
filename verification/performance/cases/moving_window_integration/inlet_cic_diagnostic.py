#!/usr/bin/env python3
"""Small algebra audit of a reconstructed held particle's moving-inlet entry.

This reads existing review JSON/YAML only. It does not run PIC, solve a field,
or modify a numerical kernel. Normalized charge coefficients are fractions of
one macro-particle's charge; they are NOT relative electromagnetic-field errors.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import yaml


def require(test, message):
    if not test:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cic(z):
    k = math.floor(z)
    f = z - k
    return {k: 1.0 - f, k + 1: f}


def translate_nodes(values, displacement):
    """The linear native-index interpolation used by the frozen provider."""
    m = math.floor(displacement)
    f = displacement - m
    result = {}
    for k, value in values.items():
        result[k + m] = result.get(k + m, 0.0) + (1.0 - f) * value
        result[k + m + 1] = result.get(k + m + 1, 0.0) + f * value
    return result


def subtract(first, second):
    return {k: first.get(k, 0.0) - second.get(k, 0.0)
            for k in sorted(set(first) | set(second))}


def l1(values):
    return math.fsum(abs(value) for value in values.values())


def z_current(start, end):
    """Pure-z depositCellSegment, summed over the four transverse vertices.

    Jz(k) = coefficient(k) * Q/(dx*dy*dt). Every segment uses the original
    whole Maxwell dt; clipping never changes this denominator. Multiplication
    by the same four transverse CIC weights gives the pointwise 3-D identity.
    """
    require(end >= start, 'This audit is restricted to nonnegative held beta')
    require(end - start <= 2.0, 'Algebra audit segment bound exceeded')
    return {k: max(0.0, min(end, k + 1.0) - max(start, float(k)))
            for k in range(math.floor(start), math.ceil(end))}


def clipping_case(start, end):
    full = z_current(start, end)
    restriction = {k: value for k, value in full.items() if k >= 0}
    clipped = z_current(max(0.0, start), max(0.0, end))
    before, after = cic(start), cic(end)
    keys = set(before) | set(after) | set(full) | {k + 1 for k in full}
    continuity = {k: after.get(k, 0.0) - before.get(k, 0.0)
                  + full.get(k, 0.0) - full.get(k - 1, 0.0) for k in keys}
    error = l1(subtract(restriction, clipped))
    residual = max((abs(value) for value in continuity.values()), default=0.0)
    require(error <= 2e-14 and residual <= 2e-14, 'Pure-z clipping identity failed')
    return {'start_cells': start, 'end_cells': end,
            'full_Jz_coefficients': full, 'window_Jz_coefficients': clipped,
            'window_restriction_L1_error': error,
            'full_continuity_max_normalized_residual': residual,
            'boundary_rho_before_Q': before.get(0, 0.0),
            'boundary_rho_after_Q': after.get(0, 0.0),
            'outside_edge_current_coefficient': full.get(-1, 0.0),
            'boundary_residual_if_external_flux_omitted':
                after.get(0, 0.0) - before.get(0, 0.0) + clipped.get(0, 0.0)}


def translation_case(z0, beta_reference, beta_particle, step):
    displacement = beta_reference * step
    translated = translate_nodes(cic(z0), displacement)
    common_beta_cic = cic(z0 + displacement)
    actual_cic = cic(z0 + beta_particle * step)
    subcell = subtract(translated, common_beta_cic)
    velocity = subtract(common_beta_cic, actual_cic)
    combined = subtract(translated, actual_cic)
    u, f = z0 - math.floor(z0), displacement - math.floor(displacement)
    amplitude = min(u * f, (1.0 - u) * (1.0 - f))
    require(abs(l1(subcell) - 4.0 * amplitude) <= 3e-13,
            'CIC commutator does not match the closed-form second difference')
    require(abs(math.fsum(combined.values())) <= 3e-13,
            'Translation must preserve total macro charge')
    return {'step': step, 'fractional_reference_shift': f,
            'beta_only_position_difference_cells': (beta_particle - beta_reference) * step,
            'reference_translated_initial_CIC': translated,
            'redeposited_CIC_common_beta': common_beta_cic,
            'redeposited_CIC_particle_beta': actual_cic,
            'subcell_commutator_L1_over_abs_Q': l1(subcell),
            'common_beta_difference_L1_over_abs_Q': l1(velocity),
            'combined_L1_over_abs_Q': l1(combined),
            'combined_signed_node_coefficients': combined,
            'commutator_second_difference_amplitude': amplitude,
            'common_beta_CIC_and_reference_same_centroid':
                abs(math.fsum(k * value for k, value in subcell.items())) <= 3e-11}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root', type=Path, required=True)
    parser.add_argument('--source-id', type=int, default=2)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.case_root.resolve()
    paths = {name: root / name for name in (
        'SINGLE_MOVING_REVIEW.json', 'SINGLE_PAIR_REVIEW.json', 'plan.json',
        'moving/run.yaml', 'moving/aprl.log')}
    review = json.loads(paths['SINGLE_MOVING_REVIEW.json'].read_text())
    pair = json.loads(paths['SINGLE_PAIR_REVIEW.json'].read_text())
    plan = json.loads(paths['plan.json'].read_text())
    config = yaml.safe_load(paths['moving/run.yaml'].read_text())
    require(config['units'] == {'length': 'm', 'time': 's'}, 'Expected explicit SI card')
    inlet = review['left_inlet_reconstruction']
    crossings = [event for event in inlet['crossings']
                 if event['source_id'] == args.source_id
                 and event['type'] == 'physical_drift_crosses_old_inlet']
    require(len(crossings) == 1, 'Expected one reconstructed pure-z inlet crossing')
    crossing = crossings[0]
    step = crossing['step']
    require(isinstance(step, int) and 1 <= step <= 4096, 'Audit step bound exceeded')
    initial_offset = review['actual_geometry_vs_template_estimate']['actual_initial_offset']
    index = inlet['initial_behind_ids'].index(args.source_id)
    z0 = initial_offset - inlet['initial_cells_below_inlet'][index]
    start = crossing['relative_before_drift_cells']
    end = crossing['relative_after_drift_before_shift_cells']
    beta_particle = end - start
    beta_reference = inlet['logged_reference_beta']
    require(-1.0 < start < 0.0 < end < 1.0 and 0.0 < beta_particle < 1.0,
            'Selected event is not the expected single-face held crossing')
    old_left = initial_offset - (step - 1)
    require(abs(z0 + (step - 1) * beta_particle - old_left - start) < 1e-10,
            'Independent reconstruction fields disagree')
    require(crossing['release_step_fraction'] > step + 1,
            'Pure-z derivation does not apply to a release within the inlet step')

    current_cases = [clipping_case(a, b) for a, b in
                     [(-1.2, -.97), (-.7, -.47), (start, end),
                      (-.05, .18), (0.0, .23), (.91, 1.14)]]
    translation_cases = [translation_case(z0, beta_reference, beta_particle, n)
                         for n in (step - 1, step, step + 1)]
    controls = {
        'integer_initial_phase_common_beta': translation_case(
            float(math.floor(z0)), beta_reference, beta_reference, step),
        'integer_translation_common_beta': translation_case(z0, .25, .25, 8),
        'equal_beta_original_phase': translation_case(z0, beta_reference, beta_reference, step)}
    require(controls['integer_initial_phase_common_beta']['subcell_commutator_L1_over_abs_Q'] < 3e-13,
            'Integer initial phase control failed')
    require(controls['integer_translation_common_beta']['subcell_commutator_L1_over_abs_Q'] < 3e-13,
            'Integer translation control failed')

    # Reconstruct only a +c longitudinal ray. This is a causal timing check,
    # not propagation of the discretized fields or proof of a radiation mode.
    match = re.search(r'field_nz=(\d+) initial_offset=(\d+)',
                      paths['moving/aprl.log'].read_text())
    require(match is not None and int(match[2]) == initial_offset, 'Geometry log mismatch')
    gamma = float(config['mesh']['boost_gamma'])
    dz = float(config['mesh']['cell_size'][2])
    boost_beta = math.sqrt(gamma - 1.0) * math.sqrt(gamma + 1.0) / gamma
    plane_z = review['field_clock_aperture']['attributes']['plane_z_m']
    outer_z = review['actual_geometry_vs_template_estimate']['actual_field_outer_lab_z_at_t0_m']
    detector_q0 = initial_offset + int(match[1]) + (plane_z - outer_z) / (gamma * dz)
    fraction = -start / beta_particle
    events = [('physical_drift_crossing', old_left, step - 1 + fraction),
              ('tangential_E_boundary_replacement', old_left, float(step)),
              ('new_left_plane_supply', old_left - 1, float(step))]
    timing = []
    for name, position, time in events:
        arrival = (detector_q0 - position + time) / (1.0 + boost_beta)
        timing.append({'event': name, 'native_z_index': position,
                       'boost_step_fraction': time,
                       'on_axis_plus_c_arrival_step': arrival,
                       'first_complete_step_at_or_after_arrival': math.ceil(arrival)})

    output = {
        'format': 'aprl-inlet-cic-algebra-v1',
        'scope': 'Read-only one-particle algebra and existing-result association; no PIC, field solve or causal-mechanism acceptance.',
        'provenance': {
            'case_root': str(root), 'script_sha256': sha(Path(__file__).resolve()),
            'inputs_sha256': {str(path): sha(path) for path in paths.values()},
            'frozen_cpp_sha256': {name: sha(root / 'source/src' / name) for name in
                ('eb_deposition.cpp', 'simulation.cpp', 'simulation_moving_window.cpp', 'moving_window_fields.cpp')},
            'coordinate_precision_limit': 'Initial coordinates come from the independent review using a rounded logged t0; approximately 1e-8 cell reconstruction error, not directly logged held trajectories.'},
        'particle': {'source_id': args.source_id, 'initial_native_z': z0,
                     'initial_subcell_phase': z0 - math.floor(z0),
                     'beta_particle': beta_particle, 'beta_reference': beta_reference,
                     'entry_step': step, 'release_step_fraction': crossing['release_step_fraction'],
                     'macro_electrons': plan['source']['macro_electrons'],
                     'represented_electrons': plan['source']['represented_electrons']},
        'clipped_current': {
            'conclusion': 'For strictly held pure-z motion, clipping at the inlet is algebraically identical to restricting full-domain J to window edges. Whole-step dt is retained. This is a derivation of the C++ formula, not an executed C++ regression.',
            'normalization': 'Jz coefficient multiplies Q/(dx*dy*dt); rho coefficient multiplies Q/(dx*dy*dz). Identical fixed x/y CIC factors make the J identity pointwise at all four transverse vertices.',
            'transverse_currents': 'Jx=Jy=0 exactly for this held segment.',
            'cases': current_cases,
            'boundary_note': 'Before centre crossing, CIC can occupy boundary node0 while all Jz is on outside edge(-1). Omitting that external flux produces an open-boundary continuity residual; it is not a missing in-window Jz. The supplied exterior fields must account for it.'},
        'reference_charge_shape': {
            'separation': 'subcell=translated_initial_CIC - CIC(z0+mean_beta*n); beta=CIC(z0+mean_beta*n)-CIC(z0+particle_beta*n); combined is their signed sum, not the sum of their L1 norms.',
            'formula': 'For initial phase u and fractional shift f, translated_CIC - redeposited_CIC = a*[1,-2,1], a=min(u*f,(1-u)*(1-f)). Charge and centroid are conserved; L1/abs(Q)=4*a.',
            'field_relation': 'In the constant-coefficient interior, native linear translation commutes with the discrete divergence. A translated Gauss-consistent initial E therefore represents translated initial rho, not generally the CIC rho of shifted particles. The actual moving field combines supplied and evolved cells, so these numbers are not its measured Gauss or E/B error.',
            'steps': translation_cases, 'controls': controls},
        'space_time_association': {
            'observed_local_field_context_from_pair_review': pair['local_max_error_context'],
            'observed_global_E_relative_L2': pair['fixed_to_moving']['fields']['electric_V_per_m']['full_output']['relative_l2_reference_norm'],
            'observed_global_B_relative_L2': pair['fixed_to_moving']['fields']['magnetic_T']['full_output']['relative_l2_reference_norm'],
            'detector_native_z_at_step0': detector_q0, 'boost_beta': boost_beta,
            'on_axis_ray_events': timing,
            'limitation': 'The source neighbourhood and predicted arrival agree with the observed peak, making inlet/reference mismatch a falsifiable candidate. Finite transverse propagation, numerical dispersion, and E/B half-time staggering are not simulated here; no amplitude or exclusive-cause prediction is made.'},
        'minimal_falsification_checks_not_run': [
            'Use the unchanged C++ depositor on one held particle in a padded full grid and a clipped grid; compare all retained native Jx/Jy/Jz for the listed cases before any field update.',
            'At steps8/9/10, compare discrete divergence of the actual queried reference E with both linearly translated prepared rho and fresh CIC of the unchanged held source. This tests the precise commutator without modifying PIC.',
            'In a small verification-only oracle, record fields immediately before/after tangential-E replacement and new-plane supply; compare against the matched full evolution. Keep the Cowan/Boris/deposition formulas unchanged.',
            'Equal-beta and integer-phase/integer-translation controls separate velocity spread from CIC phase interpolation; do not silently change the production source or kernel.'],
        'dense_followup': {
            'not_a_fix': 'Lower macro charge reduces the absolute contribution of each isolated macro-particle, but leaves the normalized one-particle CIC commutator unchanged. More particles can add coherently; aggregate field or weak-radiation improvement is not guaranteed.',
            'controlled_inputs': 'Hold total charge, distribution, grid/dt, boundary/guards, release states, origin/centering, observation plane and duration fixed. Preserve an explicit inlet witness and report its actual subcell phase; changed extrema or centering can otherwise confound an N comparison.',
            'measurements': ['Absolute local E and cB difference maps around entry witnesses, with the same spatial/time support and a fixed reference-amplitude floor.',
                             'Local Gauss/source mismatch both in absolute charge units and normalized to individual macro charge.',
                             'Matched particle energy-kick differences: RMS, tails and local cohort, plus transverse self-field response; not only whole-bunch means.',
                             'Whole-plane E/B norms and target spectral-band observables in addition to local maxima.'],
            'acceptance_scope': 'Improved accuracy at production macro weight may support that sampled regime only. Diluting a local defect in a whole-plane norm does not repair the algorithm or establish radiation purity.'},
        'unchanged': ['Cowan advance', 'Boris advance', 'charge/current deposition', 'frozen reference provider', 'all simulation inputs and outputs']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(output, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'output': str(args.output.resolve()), 'scope': output['scope'],
                      'entry_step_combined_charge_shape_L1_over_abs_Q': translation_cases[1]['combined_L1_over_abs_Q'],
                      'script_sha256': output['provenance']['script_sha256']}))


if __name__ == '__main__':
    main()
