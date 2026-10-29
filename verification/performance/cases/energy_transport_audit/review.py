#!/usr/bin/env python3
"""Independent lifecycle, unchanged-physics, energy-accounting and timing review."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import json
from pathlib import Path
import re
import statistics
import h5py
import numpy as np
import yaml
from prepare import canonical, require, save, sha

LD = np.longdouble
EPS = LD(np.finfo(np.float64).eps)


def finite(array):
    a = np.asarray(array)
    if a.dtype.names:
        return all(finite(a[name]) for name in a.dtype.names)
    return bool(np.isfinite(a).all()) if a.dtype.kind in 'fc' else True


def energy_review(path, new, export_csv=True):
    with h5py.File(path) as stream:
        group = stream['energy_ledger']; rows = group['records'][:]
        attrs = {name: value.decode() if isinstance(value, bytes) else value.item() if isinstance(value, np.generic) else value
                 for name, value in group.attrs.items()}
    require(int(attrs.get('format_version', -1)) == (5 if new else 4), 'Unexpected ledger schema version')
    moving = True  # Accepted baseline already includes all v4 moving-volume terms.
    require('boosted computational frame' in str(attrs.get('frame')), 'Unexpected ledger frame')
    require(len(rows) > 1 and np.all(np.diff(rows['step']) > 0), 'Invalid ledger step sequence')
    def a(name):
        require(name in rows.dtype.names, 'Missing ledger component: ' + name)
        return np.asarray(rows[name], dtype=LD)
    kinetic = a('particle_kinetic_active_J') + a('particle_kinetic_removed_J')
    delta_k = kinetic-kinetic[0]
    field = a('field_energy_interior_J'); delta_u = field-field[0]
    faces = a('outward_field_energy_face_J')
    outward = np.sum(faces, axis=1)
    prescribed = a('prescribed_work_J'); boundary = a('boundary_source_work_J')
    raw = delta_k+delta_u+outward-prescribed
    expected = raw-boundary
    identities = {}
    source_scale = abs(kinetic)+abs(kinetic[0])+abs(field)+abs(field[0])+np.sum(abs(faces),axis=1)+abs(prescribed)+abs(boundary)
    if moving:
        inserted = a('moving_inserted_field_J')
        cropped = a('moving_cropped_field_J')
        remap = a('moving_remap_field_change_J')
        correction = a('moving_remap_interface_correction_J')
        overwrite = a('inlet_overwrite_work_J')
        longitudinal = a('longitudinal_boundary_outward_J')
        expected -= remap+overwrite
        source_scale += abs(inserted)+abs(cropped)+abs(remap)+abs(correction)+abs(overwrite)
        identities['remap_equals_insert_minus_crop_plus_interface'] = (remap, inserted-cropped+correction)
        identities['longitudinal_groups_propagation_remap_and_overwrite'] = (longitudinal, np.sum(faces[:,4:],axis=1)-remap-overwrite)
        # Added constrained-current supply fields are handled explicitly once
        # present; they are not inferred from the residual to force a closure.
        boundary_corrected = expected.copy()
        constraint_total = np.zeros_like(expected)
        constraint_absolute_components = np.zeros_like(expected)
        for name in ('frozen_constraint_work_J', 'carrier_constraint_work_J'):
            expected -= a(name); source_scale += abs(a(name)); constraint_total += a(name)
            constraint_absolute_components += abs(a(name))
    identities['outward_equals_sum_of_faces'] = (a('outward_field_energy_J'), outward)
    identities['raw_residual_reconstructed'] = (a('balance_residual_raw_J'), raw)
    identities['corrected_residual_reconstructed'] = (a('balance_residual_J'), expected)
    identity_report = {}
    # This bound covers serialization and cancellation in independently
    # reconstructed sums; it is not a permitted physical conservation error.
    tolerance = 128*EPS*np.maximum(source_scale, LD(np.finfo(np.float64).tiny))
    for name,(reported,reconstructed) in identities.items():
        error = abs(reported-reconstructed)
        passed = bool(np.all(error <= tolerance))
        identity_report[name] = {'pass': passed, 'max_abs_difference_J': float(np.max(error)),
            'max_roundoff_bound_J': float(np.max(tolerance)), 'elements_over_bound': int(np.count_nonzero(error>tolerance))}
    if new:
        stable_k, stable_u = a('delta_particle_kinetic_box_J'), a('delta_field_energy_J')
        event_kinetic = a('particle_kinetic_active_J')+a('particle_kinetic_box_removed_event_J')
        grid, device = a('particle_grid_work_J'), a('prescribed_boris_work_J')
        magnetic, push = a('magnetic_rotation_delta_K_J'), a('boris_push_delta_K_J')
        sources = boundary + remap + overwrite + constraint_total
        absolute_work = a('boris_absolute_work_J')
        work_parts = abs(grid)+abs(device)+abs(magnetic)+abs(push)+abs(a('boris_roundoff_residual_J'))
        absolute_tolerance = 128*EPS*np.maximum(absolute_work, LD(np.finfo(float).tiny))
        identity_report['absolute_work_history_is_nonnegative_monotone_and_bounds_signed_parts'] = {
            'pass':bool(np.all(absolute_work>=0) and np.all(np.diff(absolute_work)>=0) and np.all(work_parts<=absolute_work+absolute_tolerance)),
            'final_absolute_work_J':float(absolute_work[-1]), 'maximum_signed_parts_J':float(np.max(work_parts))}
        removed = a('removed_macroparticles')
        partial = a('partial_work_substeps')
        before_removal = np.ones(len(rows),dtype=bool) if not np.any(removed>0) else np.arange(len(rows))<np.flatnonzero(removed>0)[0]
        identity_report['partial_substeps_only_after_removal'] = {
            'pass':bool(np.all(partial[before_removal]==0) and np.all(partial>=0) and np.all(np.diff(partial)>=0) and np.all(partial==np.floor(partial))),
            'nonzero_samples_before_first_removed':int(np.count_nonzero(partial[before_removal])),
            'final_partial_work_substeps':int(partial[-1]),
            'note':'No positive count required after removal; a terminal crossing can coincide with a complete substep.'}
        additions = {
            'stable_delta_particle_matches_inventory': (stable_k, event_kinetic-event_kinetic[0]),
            'stable_delta_removal_interpolation': (a('particle_removal_interpolation_delta_J'), a('particle_kinetic_box_removed_event_J')-a('particle_kinetic_removed_J')),
            'stable_delta_field_matches_inventory': (stable_u, delta_u),
            'boris_kick_work_decomposition': (a('boris_roundoff_residual_J'), push-grid-device-magnetic),
            'particle_work_residual': (a('particle_work_balance_residual_J'), stable_k-grid-device-magnetic),
            'field_exchange_residual': (a('field_exchange_residual_J'), stable_u+outward-sources+grid),
            'boris_total_residual': (a('balance_residual_boris_J'), stable_k+stable_u+outward-sources-device),
            'split_particle_field_residual': (a('balance_residual_boris_J'), a('particle_work_balance_residual_J')+a('field_exchange_residual_J')+magnetic)}
        for name, (reported, reconstructed) in additions.items():
            # Inventory checks need cancellation bounds; work-only checks use
            # the work scale, never the large initial beam kinetic energy.
            bound = tolerance if name.startswith('stable_delta') else 128*EPS*np.maximum(
                abs(stable_k)+abs(stable_u)+np.sum(abs(faces),axis=1)+abs(boundary)+abs(remap)+abs(overwrite)+
                constraint_absolute_components+abs(grid)+abs(device)+abs(magnetic)+abs(push), LD(np.finfo(float).tiny))
            if name == 'boris_kick_work_decomposition':
                # This sum is measured independently at every Boris substep,
                # before cancellations. Never infer it from this discrepancy.
                bound = absolute_tolerance
            err = abs(reported-reconstructed)
            identity_report[name] = {'pass':bool(np.all(err<=bound)), 'max_abs_difference_J':float(np.max(err)),
                'max_roundoff_bound_J':float(np.max(bound)), 'elements_over_bound':int(np.count_nonzero(err>bound))}
        lab = a('particle_kinetic_lab_active_J')+a('particle_kinetic_lab_removed_J')
        cpz = a('particle_cPz_box_active_J')+a('particle_cPz_box_removed_J')
        lorentz = a('delta_particle_kinetic_lab_J')-a('boost_gamma')*(stable_k+a('boost_beta')*a('delta_particle_cPz_box_J'))
        bound = 128*EPS*np.maximum(abs(lab)+abs(lab[0])+a('boost_gamma')*(abs(kinetic)+abs(kinetic[0])+abs(cpz)+abs(cpz[0])), LD(np.finfo(float).tiny))
        err = abs(a('particle_lorentz_delta_residual_J')-lorentz)
        identity_report['particle_lorentz_delta_record_identity'] = {'pass':bool(np.all(err<=bound)),
            'max_abs_difference_J':float(np.max(err)), 'max_roundoff_bound_J':float(np.max(bound)),
            'elements_over_bound':int(np.count_nonzero(err>bound))}
    residual = a('balance_residual_J')
    exchange = abs(delta_k)+abs(delta_u)+np.sum(abs(faces),axis=1)+abs(prescribed)+abs(boundary)
    code_exchange = abs(delta_k)+abs(delta_u)+abs(outward)+abs(prescribed)+abs(boundary)
    if moving:
        exchange += abs(remap)+abs(overwrite)+constraint_absolute_components
        code_exchange += abs(remap)+abs(overwrite)+constraint_absolute_components
    relative = np.divide(residual, exchange, out=np.zeros_like(residual), where=exchange>0)
    code_relative = np.divide(expected, code_exchange, out=np.zeros_like(expected), where=code_exchange>0)
    reported_relative = a('relative_balance_exchange')
    relative_error = abs(reported_relative-code_relative)
    # Propagate the absolute source-scale summation bound through division.
    ratio_tolerance = np.zeros_like(tolerance)
    nonzero = code_exchange>0
    ratio_tolerance[nonzero] = tolerance[nonzero]/code_exchange[nonzero]*(2+abs(code_relative[nonzero]))+128*EPS
    identity_report['exchange_normalized_residual_reconstructed'] = {
        'pass': bool(np.all(relative_error<=ratio_tolerance)),
        'max_abs_difference': float(np.max(relative_error)), 'max_roundoff_bound': float(np.max(ratio_tolerance)),
        'elements_over_bound': int(np.count_nonzero(relative_error>ratio_tolerance))}
    names = [name for name in rows.dtype.names if rows.dtype[name].shape == () and
             (name.endswith('_J') or name in ('step','time_box_s','closure_valid','removed_macroparticles','partial_work_substeps','released_charge_C','boost_gamma','boost_beta'))]
    if export_csv:
        with path.with_suffix('.csv').open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(names)
            writer.writerows([[row[name].item() for name in names] for row in rows])
    first_exit = np.flatnonzero(rows['removed_macroparticles']>0)
    windows = {'full': np.ones(len(rows), dtype=bool)}
    if len(first_exit):
        windows['before_first_exit'] = np.arange(len(rows))<first_exit[0]
        windows['from_first_exit'] = np.arange(len(rows))>=first_exit[0]
    phases = {}
    for name, mask in windows.items():
        indices = np.flatnonzero(mask)
        if not len(indices):
            continue
        last = indices[-1]
        phases[name] = {'first_step': int(rows['step'][indices[0]]), 'last_step': int(rows['step'][last]),
            'maximum_abs_residual_J': float(np.max(abs(residual[mask]))), 'last_residual_J': float(residual[last]),
            'maximum_abs_code_exchange_normalized_residual': float(np.max(abs(code_relative[mask]))),
            'maximum_abs_facewise_exchange_normalized_residual': float(np.max(abs(relative[mask])))}
    decomposition = {'raw_final_J': float(raw[-1]), 'raw_max_abs_J': float(np.max(abs(raw))),
        'corrected_final_J': float(expected[-1]), 'corrected_max_abs_J': float(np.max(abs(expected)))}
    if moving:
        decomposition.update(boundary_only_final_J=float(boundary_corrected[-1]),
            boundary_only_max_abs_J=float(np.max(abs(boundary_corrected))),
            constraint_work_final_J=float(constraint_total[-1]))
    physical = {}
    net_change = stable_k if new else delta_k
    for label, values in [('legacy', residual)] + ([('boris', a('balance_residual_boris_J'))] if new else []):
        peak = int(np.argmax(abs(values)))
        denominator_change = stable_k if label=='boris' else delta_k
        physical[label] = {'final_J':float(values[-1]), 'max_abs_J':float(abs(values[peak])),
            'max_abs_step':int(rows['step'][peak]),
            'final_over_abs_net_particle_change':float(values[-1]/abs(denominator_change[-1])) if denominator_change[-1] else None,
            'peak_over_abs_final_particle_change':float(abs(values[peak])/abs(denominator_change[-1])) if denominator_change[-1] else None}
    return {'physical_imbalance': physical, 'checks_pass': all(x['pass'] for x in identity_report.values()), 'frame': attrs['frame'],
        'metadata': attrs, 'identities': identity_report, 'closure_valid_values': np.unique(rows['closure_valid']).tolist(),
        'final_components': {name: rows[name][-1].item() for name in names},
        'first_removed_step': int(rows['step'][first_exit[0]]) if len(first_exit) else None,
        'phases': phases, 'final_delta_particle_J': float(net_change[-1]), 'legacy_final_delta_particle_J': float(delta_k[-1]), 'final_delta_field_J': float(delta_u[-1]),
        'residual_decomposition': decomposition,
        'initial_kinetic_J': float(kinetic[0]), 'initial_field_J': float(field[0]),
        'conservation_acceptance': 'NOT_ASSIGNED: physical residual is reported, not gated by a newly chosen percentage.',
        'roundoff_identity_bound': 'Inventory identities use128*float64_eps*sum(abs(source energy terms)); Boris decomposition uses128*float64_eps*independently accumulated boris_absolute_work_J across actual substeps, including canceled historical work, never residual-derived.',
        'exchange_normalization': 'Code denominator uses abs(sum of six outward faces), abs(frozen work)+abs(carrier work); facewise denominator replaces abs(sum faces) by sum(abs(faces)).',
        'scope': 'Net signed control-volume energy, boost frame; no FFT/windowing or radiation-only interpretation.'}


def scalar_attrs(group):
    return {name: value.decode() if isinstance(value, bytes) else value.item() if isinstance(value, np.generic) else value
            for name, value in group.attrs.items()}


def plane_review(path, plan, card):
    """Fixed-lab transport, with unavailable field samples explicitly masked."""
    with h5py.File(path) as stream:
        require(list(stream) == ['energy_plane'], 'Unexpected energy-plane group')
        group = stream['energy_plane']; rows = group['records'][:]; attrs = scalar_attrs(group)
        require(int(group['complete'][()]) == 1 and int(group['committed_records'][()]) == len(rows),
                'Incomplete energy-plane file')
    require(int(attrs['format_version']) == 1 and len(rows)>1, 'Energy-plane schema/sample count')
    require(np.all(np.diff(rows['time_box_s'])>0) and np.all(np.diff(rows['time_lab_s'])>0),
            'Nonmonotone energy-plane times')
    available = rows['field_available'] == 1
    require(np.all((rows['field_available']==0)|available), 'Invalid availability flag')
    masked = {'field_power_lab_W','physical_field_power_lab_W',
              'field_energy_per_length_lab_J_per_m','physical_field_energy_per_length_lab_J_per_m'}
    for name in rows.dtype.names:
        require(finite(rows[name][available] if name in masked else rows[name]), 'Nonfinite available energy quantity: '+name)
        if name in masked:
            require(np.isnan(rows[name][~available]).all(), 'Unavailable field must be explicit NaN: '+name)
    for name in ('field_sample_count','outside_field_window_steps','field_gap_count',
                 'particle_downstream_count','particle_upstream_count',
                 'physical_particle_downstream_count','physical_particle_upstream_count'):
        require(np.all(rows[name][1:]>=rows[name][:-1]), 'Nonmonotone cumulative counter: '+name)
    require(int(rows['particle_downstream_count'][-1]) == plan['particles'] and
            int(rows['particle_upstream_count'][-1]) == 0, 'Missing/duplicate/backwards full particle cohort')
    require(np.all(rows['physical_particle_downstream_count']<=rows['particle_downstream_count']) and
            np.all(rows['physical_particle_upstream_count']<=rows['particle_upstream_count']), 'Invalid aperture counts')
    def a(name):
        require(name in rows.dtype.names, 'Missing energy-plane quantity: '+name)
        value = rows[name].astype(LD)
        if name.startswith(('particle_','physical_particle_')) and name.endswith('_J'):
            require(name+'_low' in rows.dtype.names, 'Missing compensated low component: '+name)
            value += rows[name+'_low'].astype(LD)
        return value
    checks = {}
    def identity(name, lhs, rhs, scale):
        error = abs(lhs-rhs); bound = 256*EPS*np.maximum(scale,LD(np.finfo(float).tiny))
        checks[name] = {'pass':bool(np.all(error<=bound)), 'max_abs_difference':float(np.max(error)),
                        'max_roundoff_bound':float(np.max(bound)), 'elements_over_bound':int(np.count_nonzero(error>bound))}
    gamma, beta = LD(attrs['boost_gamma']), LD(attrs['boost_beta'])
    require(float(gamma)==card['mesh']['boost_gamma'], 'Plane boost differs from physical card')
    # Constants used by the existing code are IEEE double, so convert those
    # exact values into extended precision rather than changing the model.
    me, c, qe = LD(float(9.1093837139e-31)), LD(299792458), LD(float(1.602176634e-19))
    particles = {}
    for prefix in ('particle_', 'physical_particle_'):
        kl, el, pl, kb, eb, pb = (a(prefix+name) for name in
            ('kinetic_lab_J','total_lab_J','cpz_lab_J','kinetic_box_J','total_box_J','cpz_box_J'))
        charge = a(prefix+'charge_C'); rest = -charge/qe*me*c*c
        identity(prefix+'lab_total_kinetic_rest',el-kl,rest,abs(el)+abs(kl)+abs(rest))
        identity(prefix+'box_total_kinetic_rest',eb-kb,rest,abs(eb)+abs(kb)+abs(rest))
        identity(prefix+'lorentz_energy',el,gamma*(eb+beta*pb),abs(el)+gamma*(abs(eb)+abs(beta*pb)))
        identity(prefix+'lorentz_cpz',pl,gamma*(pb+beta*eb),abs(pl)+gamma*(abs(pb)+abs(beta*eb)))
        particles[prefix.rstrip('_')] = {name:float(a(prefix+name)[-1]) for name in
            ('kinetic_lab_J','total_lab_J','cpz_lab_J','kinetic_box_J','total_box_J','cpz_box_J','charge_C')}
        particles[prefix.rstrip('_')].update(downstream_count=int(rows[prefix+'downstream_count'][-1]),
                                           upstream_count=int(rows[prefix+'upstream_count'][-1]))
    expected_charge = -LD(card['beam']['input']['electrons'])*qe
    identity('full_cohort_charge',a('particle_charge_C')[-1],expected_charge,abs(expected_charge))
    require(np.any(available), 'Plane never sampled an available field')
    samples = rows['field_sample_count']; adjacent = available[1:] & available[:-1] & (samples[1:]==samples[:-1]+1)
    require(np.any(adjacent), 'No adjacent field steps for independent quadrature check')
    dt = np.diff(a('time_lab_s'))
    for prefix in ('', 'physical_'):
        power, energy = a(prefix+'field_power_lab_W'), a(prefix+'field_energy_lab_J')
        trapezoid = (power[1:][adjacent]+power[:-1][adjacent])*dt[adjacent]/2
        delta = np.diff(energy)[adjacent]
        identity(prefix+'field_trapezoid',delta,trapezoid,
                 abs(energy[1:][adjacent])+abs(energy[:-1][adjacent])+abs(trapezoid))
        unavailable = ~available[1:]
        if np.any(unavailable):
            identity(prefix+'no_flux_integration_while_unavailable',np.diff(energy)[unavailable],np.zeros(np.count_nonzero(unavailable),dtype=LD),
                     abs(energy[1:][unavailable])+abs(energy[:-1][unavailable]))
    require(int(rows['field_gap_count'][-1]) == 0, 'Discontinuous field coverage within moving-plane crossing')
    valid = np.flatnonzero(available)
    coverage = {'field_sample_count':int(samples[-1]), 'written_available_samples':int(available.sum()),
        'outside_field_window_steps':int(rows['outside_field_window_steps'][-1]),
        'field_gap_count':int(rows['field_gap_count'][-1]), 'quadrature_intervals_checked':int(adjacent.sum()),
        'first_field_time_lab_s':float(rows['first_field_time_lab_s'][-1]),
        'last_field_time_lab_s':float(rows['last_field_time_lab_s'][-1]),
        'first_written_available_time_lab_s':float(rows['time_lab_s'][valid[0]]),
        'last_written_available_time_lab_s':float(rows['time_lab_s'][valid[-1]]),
        'complete_means_full_field_history':False}
    return {'checks_pass':all(v['pass'] for v in checks.values()), 'metadata':attrs, 'identities':checks,
        'coverage':coverage, 'particles':particles,
        'signed_field_energy_lab_J':float(rows['field_energy_lab_J'][-1]),
        'signed_physical_field_energy_lab_J':float(rows['physical_field_energy_lab_J'][-1]),
        'peak_field_energy_per_length_lab_J_per_m':float(np.max(rows['field_energy_per_length_lab_J_per_m'][available])),
        'peak_physical_field_energy_per_length_lab_J_per_m':float(np.max(rows['physical_field_energy_per_length_lab_J_per_m'][available]))}


def plane_pair_review(directory, plan, card):
    filenames = plan['energy_plane_files']
    reports = {name:plane_review(directory/'detectors'/filenames[name],plan,card) for name in ('entry','exit')}
    entry, output = reports['entry'], reports['exit']
    # Reopen only these small diagnostic files; do not round compensated totals
    # through JSON/double before taking the small inlet-to-outlet difference.
    terminal = {}
    for name in reports:
        with h5py.File(directory/'detectors'/filenames[name]) as stream:
            terminal[name] = stream['energy_plane/records'][-1]
    differences = {}
    for prefix in ('particle_', 'physical_particle_'):
        def precise(side, field):
            row = terminal[side]; return LD(row[prefix+field])+LD(row[prefix+field+'_low'])
        changes = {field:float(precise('exit',field)-precise('entry',field)) for field in
                   ('kinetic_lab_J','total_lab_J','cpz_lab_J','kinetic_box_J','total_box_J','cpz_box_J')}
        differences[prefix.rstrip('_')] = changes
    for report in reports.values():
        require(report['particles']['particle']['downstream_count']==plan['particles'], 'Different full transport cohort')
    flux = LD(terminal['entry']['field_energy_lab_J'])-LD(terminal['exit']['field_energy_lab_J'])
    physical_flux = LD(terminal['entry']['physical_field_energy_lab_J'])-LD(terminal['exit']['physical_field_energy_lab_J'])
    loss = -(LD(terminal['exit']['particle_kinetic_lab_J'])+LD(terminal['exit']['particle_kinetic_lab_J_low'])-
             LD(terminal['entry']['particle_kinetic_lab_J'])-LD(terminal['entry']['particle_kinetic_lab_J_low']))
    return {'checks_pass':all(v['checks_pass'] for v in reports.values()), 'planes':reports,
        'particle_exit_minus_entry':differences, 'full_particle_kinetic_loss_lab_J':float(loss),
        'signed_field_energy_entry_minus_exit_lab_J':float(flux),
        'signed_physical_field_energy_entry_minus_exit_lab_J':float(physical_flux),
        'full_particle_plus_field_entry_minus_exit_lab_J':float(loss+flux),
        'same_full_particle_count_and_charge':True,
        'legacy_plane_ID_coverage_checked':bool(card['trajectory']['enabled']),
        'new_plane_IDs_are_stored':False,
        'same_physical_aperture_cohort_proven':False,
        'interpretation':'Fixed-lab-plane transport comparison, not a conservation residual: same-volume storage, lateral flux and sources are not all included; total self fields are not radiation-only.',
        'precision':'Particle plane high+low reconstructed in longdouble before taking energy differences; no sigma-difference estimate of particlewise loss spread.'}


def inspect(directory, item, plan, export_csv=True):
    log = (directory/'aprl.log').read_text(errors='replace')
    require('Configured stop reached:' in log and 'phase=mpi-finalize-complete' in log, 'Missing normal physical completion')
    require(not re.search(r'Simulation failure|\[shutdown-watchdog\]|DUE TO TIME LIMIT|Segmentation fault|OpenFabrics device',log), 'Unexpected failure/warning')
    moving = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)',log)
    require(moving and int(moving[1]) == item['steps'] and int(moving[2]) == plan['particles'], 'Incomplete moving lifecycle')
    released = re.search(r'Frozen lab-plane injection: activations=(\d+),.*still frozen=(\d+)',log)
    require(released and int(released[1]) == plan['particles'] and int(released[2]) == 0, 'Incomplete release lifecycle')
    require('[performance]' not in log, 'Profile must be disabled')
    files = sorted(directory.glob('detectors/*.h5'))+sorted(directory.glob('trajectory/*.h5'))+sorted(directory.glob('ledger/*.h5'))
    card = yaml.safe_load((directory/'run.yaml').read_text())
    science = card['trajectory']['enabled']
    expected_files = (20 if science else 1)+(1 if item['ledger'] else 0)+(2 if item['energy_planes'] else 0)
    require(len(files) == expected_files, 'Unexpected output file count')
    terminal_counts = np.zeros(plan['particles'],dtype=np.int64)
    checks = {}
    for path in files:
        if path.name.endswith('-energy.h5'):
            continue
        with h5py.File(path) as stream:
            require(len(stream) == 1, 'Unexpected HDF5 group count')
            group = stream[next(iter(stream))]
            require(int(group['complete'][()]) == 1, 'Incomplete HDF5: '+str(path))
            for committed, values in (('committed_records','records'),('committed_samples','time_s')):
                if committed in group:
                    require(int(group[committed][()]) == len(group[values]), 'Truncated HDF5')
            for name,dataset in group.items():
                if not isinstance(dataset,h5py.Dataset):
                    continue
                chunks = [dataset[()]] if not dataset.shape else (dataset[i:i+4] if name=='fields' else dataset[i:i+65536]
                    for i in range(0,len(dataset),4 if name=='fields' else 65536))
                for chunk in chunks:
                    require(finite(chunk),'Nonfinite output: '+str(path)+'/'+name)
                    if path.parent.name=='trajectory' and name=='records':
                        final = chunk[chunk['event_type']!=0]
                        require(np.all(final['event_type']==2) and np.all(final['boundary_face']==5), 'Unexpected terminal event')
                        require(np.all((final['source_id']>=1)&(final['source_id']<=plan['particles'])), 'Unexpected terminal ID')
                        np.add.at(terminal_counts,final['source_id'].astype(np.int64)-1,1)
            if path.parent.name=='detectors' and path.stem!='exit-field':
                rows = group['records'][:]
                require(np.array_equal(np.sort(rows['source_id']),np.arange(1,plan['particles']+1)), 'Particle plane coverage differs')
                require(np.all(rows['crossing_direction']==1),'Unexpected backward crossing')
            checks[str(path.relative_to(directory))] = {'complete':True,'bytes':path.stat().st_size}
    if science:
        require(np.all(terminal_counts==1),'Missing/duplicate particle terminal events')
    points = {int(step):float(elapsed) for step,elapsed in re.findall(r'\[resource\] phase=(?:progress|final) step=(\d+) elapsed_s=([\d.eE+\-]+)',log)}
    require(40 in points and 240 in points and item['steps'] in points and points[240]>points[40], 'Missing timing interval')
    result = {'checks_pass':True,'steps':int(moving[1]),'particle_exits':int(moving[2]),'cropped_carriers':int(moving[3]),
        'files':checks,'progress_elapsed_s':points,'steps_40_240_seconds_per_step':(points[240]-points[40])/200,
        'all_steps_seconds_per_step':points[item['steps']]/item['steps'],'profile_enabled':False}
    preparation = json.loads((directory/'reference/preparation.json').read_text())
    result['reference_preparation_seconds'] = preparation['preparationSeconds']
    launched = json.loads((directory/'aprl.log.launch.json').read_text())
    result['launcher_seconds_including_preparation_and_close'] = launched['seconds']
    if item['ledger']:
        result['energy'] = energy_review(directory/'ledger/energy-ledger.h5',item['binary']=='new', export_csv=export_csv)
        require(result['energy']['checks_pass'], 'Energy-record algebra does not match independent reconstruction')
        work = result['energy']['final_components']['prescribed_work_J']
        require(work==0. if item['magnet']=='free-drift' else work!=0., 'External-device work coverage differs')
    if item['energy_planes']:
        result['energy_planes'] = plane_pair_review(directory, plan, card)
        require(result['energy_planes']['checks_pass'], 'Energy plane audit failed')
    return result


def compare(first,second):
    require(canonical(yaml.safe_load((first/'run.yaml').read_text())) ==
            canonical(yaml.safe_load((second/'run.yaml').read_text())), 'Comparison cards differ beyond ledger toggle/output paths')
    relative = [p.relative_to(first) for folder in ('detectors','trajectory') for p in sorted((first/folder).glob('*.h5')) if not p.name.endswith('-energy.h5')]
    require(relative == [p.relative_to(second) for folder in ('detectors','trajectory') for p in sorted((second/folder).glob('*.h5')) if not p.name.endswith('-energy.h5')], 'Compared outputs differ')
    result = {}
    for rel in relative:
        equal = True; samples = 0
        with h5py.File(first/rel) as fa,h5py.File(second/rel) as fb:
            a,b = fa[next(iter(fa))],fb[next(iter(fb))]
            for name in ('time_s','fields','records'):
                if name not in a:
                    require(name not in b,'Compared dataset absent')
                    continue
                require(name in b and a[name].shape == b[name].shape and a[name].dtype == b[name].dtype,'Compared schema differs')
                # Detector record append order may vary with rank; sort the
                # cohort by immutable source ID. Trajectory rank ownership is fixed.
                if rel.parts[0]=='detectors' and name=='records':
                    x,y = a[name][:],b[name][:]
                    x=x[np.argsort(x['source_id'])];y=y[np.argsort(y['source_id'])]
                    equal &= x.tobytes()==y.tobytes(); samples+=len(x)
                else:
                    for i in range(0,len(a[name]),4 if name=='fields' else 65536):
                        n=4 if name=='fields' else 65536
                        x,y=a[name][i:i+n],b[name][i:i+n]
                        equal &= x.tobytes()==y.tobytes();samples+=len(x)
        result[str(rel)] = {'bitwise_equal':bool(equal),'records_or_samples_compared':samples}
    return {'checks_pass':all(v['bitwise_equal'] for v in result.values()),'files':result}



def data_directory(item):
    return Path(item.get('reused_from', item['directory']))


def inspect_reused(item, plan):
    require(item['variant'] in ('old-on','new-on','new-off') and not item['energy_planes'], 'Ineligible reuse variant')
    require(item['reused_binary_sha256']==plan['binary_sha256'], 'Reused binary identity differs')
    require(all(sha(path)==value for path,value in item['reused_evidence_sha256'].items()), 'Reused provenance changed')
    original = data_directory(item)
    require(all(sha(original/relative)==value for relative,value in item['reused_output_sha256'].items()), 'Reused output changed')
    require(canonical(yaml.safe_load((original/'run.yaml').read_text())) ==
            canonical(yaml.safe_load((Path(item['directory'])/'run.yaml').read_text())), 'Reused card differs')
    result = inspect(original, item, plan, export_csv=False)
    result['reuse'] = {'reused_from':str(original), 'original_nodes':item['reused_original_nodes'],
                      'original_binding':item['reused_original_binding'], 'new_computation_launched':False}
    return result


def stage_review(root,stage,plan=None):
    plan = plan or json.loads((root/'plan.json').read_text())
    items=plan['cases'][stage]
    cases={item['name']:json.loads((Path(item['directory'])/'CASE_REVIEW.json').read_text()) for item in items}
    comparisons={}
    for magnet in ('undulator','free-drift'):
        for repeat in sorted({item['repeat'] for item in items}):
            group=[item for item in items if item['magnet']==magnet and item['repeat']==repeat]
            first=next(x for x in group if x['variant']==('new-on' if stage=='mpi' else 'old-on'))
            for other in group:
                if other==first:
                    continue
                comparisons[other['name']]=compare(data_directory(first),data_directory(other))
    result={'format':'aprl-energy-transport-review-v1','stage':stage,'binary_sha256':plan['binary_sha256'],
        'checks_pass':all(x['checks_pass'] for x in cases.values()) and all(x['checks_pass'] for x in comparisons.values()),
        'cases':cases,'physical_output_comparisons':comparisons,'energy_conservation_accepted':False,
        'reused_case_provenance':{item['name']:{'reused_from':item['reused_from'], 'original_nodes':item['reused_original_nodes']} for item in items if item.get('reused_from')},
        'note':'checks_pass covers completion, exact unchanged physical outputs and accounting identities; it does not declare a small physical energy residual.'}
    if stage=='performance':
        metrics={}
        for magnet in ('undulator','free-drift'):
            samples={v:[cases[x['name']]['steps_40_240_seconds_per_step'] for x in items if x['magnet']==magnet and x['variant']==v]
                     for v in ('old-on','new-on','new-off','new-audit-on')}
            summary={v:{'samples_s_per_step':a,'median_s_per_step':statistics.median(a),'min_s_per_step':min(a),'max_s_per_step':max(a)} for v,a in samples.items()}
            ratios={}
            for name,num,den in (('increment_on_new_vs_old','new-on','old-on'),('new_ledger_on_vs_off','new-on','new-off'),
                                 ('new_planes_on_vs_no_planes','new-audit-on','new-on'),('new_all_diagnostics_on_vs_off','new-audit-on','new-off')):
                paired=[100*(a/b-1) for a,b in zip(samples[num],samples[den])]
                ratios[name]={'paired_percent_changes':paired,'median_paired_percent_change':statistics.median(paired),
                    'ratio_of_medians_percent':100*(statistics.median(samples[num])/statistics.median(samples[den])-1)}
            metrics[magnet]={'configurations':summary,'comparisons':ratios}
        result['performance']=metrics
        result['performance_scope']='Same-node fixed16x8, no profile, complete257-step small case; 40–240 steps exclude initialization and early warmup. Not production overhead.'
    save(root/(stage.upper()+'_REVIEW.json'),result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--stage',choices=('science','performance','mpi'),required=True)
    args=parser.parse_args()
    result=stage_review(args.root.resolve(),args.stage)
    require(result['checks_pass'],'Stage review failed')
