#!/usr/bin/env python3
"""Audit a completed full moving-window run on its Slurm compute node.

PASS means configured completion and complete, finite diagnostic output. It does
not certify pure radiation, energy closure, or agreement with the old method.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

import h5py
import numpy as np
import yaml

C = 299792458.0


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    sha = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def string(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def finite(value, name):
    value = np.asarray(value)
    if value.dtype.names:
        for key in value.dtype.names:
            finite(value[key], name + '/' + key)
    elif value.dtype.kind in 'fc':
        require(np.isfinite(value).all(), 'Nonfinite data: ' + name)


def check_identity(group, manifest):
    for key in ('run_id', 'source_revision', 'configuration_digest_fnv1a64'):
        require(string(group.attrs.get(key, '')) == str(manifest[key]),
                'Diagnostic identity mismatch: ' + key)


def check_time_axis(times, cadence):
    require(len(times) > 1 and np.isfinite(times).all() and np.all(np.diff(times) > 0),
            'Invalid field time axis')
    tolerance = 8 * abs(np.spacing(float(np.max(np.abs(times))))) + 64 * np.finfo(float).eps * abs(cadence)
    require(np.all(np.abs(np.diff(times) - cadence) <= tolerance),
            'Field cadence differs beyond absolute timestamp rounding')
    # The solver accumulates box time; bound cumulative floating-point addition,
    # as well as the final conversion to an absolute laboratory timestamp.
    span_tolerance = tolerance + 2 * len(times) * np.finfo(float).eps * len(times) * abs(cadence)
    expected = np.arange(len(times), dtype=np.longdouble) * np.longdouble(cadence)
    residual = times.astype(np.longdouble) - np.longdouble(times[0]) - expected
    require(np.max(np.abs(residual)) <= span_tolerance, 'Field time span disagrees with completed steps')
    return tolerance, span_tolerance


def check_binding(path, plan):
    log = path.read_text()
    rows = re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)', log)
    ranks, threads = plan['ranks'], plan['threads_per_rank']
    require('[binding] result=PASS failures=0' in log, 'Binding probe did not pass')
    require(len(rows) == ranks * threads and all(r[-1] == 'PASS' for r in rows), 'Incomplete thread binding')
    require({(int(r[0]), int(r[1])) for r in rows} ==
            {(r, t) for r in range(ranks) for t in range(threads)}, 'Binding rank/thread coverage differs')
    require(len({r[2] for r in rows}) == plan['nodes'], 'Binding host count differs')
    require(len({(r[2], r[4], r[5]) for r in rows}) == ranks * threads, 'Binding shares physical cores')
    return dict(ranks=ranks, threads_per_rank=threads, hosts=sorted({r[2] for r in rows}),
                physical_cores=ranks * threads, log_sha256=digest(path))


def check_window_csv(directory, ranks, final_step, interval, geometry):
    expected_files = [directory / ('window-rank-%d.csv' % rank) for rank in range(ranks)]
    require(set(directory.glob('window-rank-*.csv')) == set(expected_files), 'Missing or unexpected moving-volume rank CSV')
    expected_steps = np.asarray(list(range(interval, final_step + 1, interval)) + [final_step])
    final_rows, maximum_error = [], 0.
    for path in expected_files:
        rows = np.atleast_1d(np.genfromtxt(path, names=True, delimiter=','))
        finite(rows, str(path))
        require(np.array_equal(rows['step'], expected_steps), 'Moving-volume step coverage differs: ' + path.name)
        require(np.all(rows['energy_measured'] == 1), 'Moving-volume cell energy was not measured')
        expected_origin = np.longdouble(geometry['reference_origin_z_m']) + (
            np.longdouble(geometry['initial_offset']) - rows['step'].astype(np.longdouble)) * np.longdouble(geometry['dz_m'])
        error = float(np.max(np.abs(rows['origin_box_z_m'].astype(np.longdouble) - expected_origin)))
        tolerance = 128 * np.finfo(float).eps * max(1., float(np.max(np.abs(expected_origin))))
        require(error <= tolerance, 'Window did not shift by one -z cell per step: ' + path.name)
        maximum_error = max(maximum_error, error)
        for key in ('particle_exits', 'cropped_carriers', 'inserted_field_J', 'cropped_field_J', 'shift_s', 'migration_s'):
            require(np.all(rows[key] >= 0) and np.all(np.diff(rows[key]) >= 0), 'Invalid cumulative window value: ' + key)
        for key in ('particle_exits', 'cropped_carriers'):
            require(np.all(rows[key] == np.floor(rows[key])), 'Noninteger window count: ' + key)
        if len(rows) > 1 and rows[-1]['step'] == rows[-2]['step']:
            require(all(rows[-1][key] == rows[-2][key] for key in rows.dtype.names), 'Final repeated CSV row differs')
        final_rows.append({key: float(rows[-1][key]) for key in rows.dtype.names})
    return dict(ranks=ranks, rows_per_rank=len(expected_steps), final_step=final_step,
                cell_shift_per_step=-1, maximum_origin_error_m=maximum_error,
                total_particle_exits=int(sum(r['particle_exits'] for r in final_rows)),
                total_cropped_carriers=int(sum(r['cropped_carriers'] for r in final_rows)),
                inserted_field_J=sum(r['inserted_field_J'] for r in final_rows),
                cropped_field_J=sum(r['cropped_field_J'] for r in final_rows), final_rank_rows=final_rows)


def check(root):
    report = dict(format='aprl-moving-production-completion-v1', pass_=False,
                  output_transactions_complete=False, pure_radiation_validated=False,
                  moving_volume_energy_closure_validated=False, particle_plane_output_required=False)
    report['pass'] = False
    try:
        require(bool(os.environ.get('SLURM_JOB_ID')), 'Full output scan must run inside a Slurm allocation')
        plan = json.loads((root / 'run-plan.json').read_text())
        card = yaml.safe_load((root / 'run.yaml').read_text())
        geometry = plan['geometry']
        manifest = yaml.safe_load((root / 'run-manifest.yaml').read_text())
        report['case'] = plan['case']
        require(plan['format'] == 'aprl-moving-production-v1', 'Unexpected production plan')
        require(plan['case'] in ('undulator_4p5m', 'free_drift_0p5m'), 'Unexpected production case')
        require((plan['nodes'], plan['ranks'], plan['threads_per_rank']) == (1, 16, 8), 'Frozen layout differs')
        frozen = json.loads((root / 'frozen.json').read_text())
        for name, expected in frozen.items():
            require(digest(name) == expected, 'Frozen file changed: ' + name)
        require(digest(plan['binary']) == plan['binary_sha256'], 'Binary hash differs')
        preflight = json.loads(Path(plan['preflight_pass']).read_text())
        require(preflight['pass'] and preflight['binary_sha256'] == plan['binary_sha256'] and
                preflight['input_sha256'] == plan['input_sha256'] and
                preflight['card_sha256'][str(root / 'run.yaml')] == digest(root / 'run.yaml'),
                'Preflight does not cover this binary, input and card')
        identity = json.loads((root / 'identity-check.json').read_text())
        require(identity['pass'] and identity['unchanged_input_stat'] and
                identity['binary_sha256'] == plan['binary_sha256'] and
                identity['input_sha256_preflight'] == plan['input_sha256'], 'Production startup identity did not pass')
        stat = Path(plan['input_file']).stat()
        require((stat.st_size, stat.st_mtime_ns) == (plan['input_size_bytes'], plan['input_mtime_ns']), 'Input changed since full preflight')
        report['binding'] = check_binding(root / 'binding.log', plan)
        require(manifest['source_revision'] == plan['compiled_revision'] and manifest['mpi_size'] == plan['ranks'], 'Manifest binary/layout differs')
        require(yaml.safe_load(manifest['configuration_yaml']) == card, 'Manifest configuration differs from frozen card')
        require(card['units'] == {'length': 'm', 'time': 's'}, 'Expected SI card')
        watchdog = manifest['resolved_runtime']['shutdown_watchdog']
        require(watchdog['enabled'] and watchdog['timeout_seconds'] == 600 and
                watchdog['finalize_timeout_seconds'] == 60, 'MPI shutdown protection changed')
        held = manifest['resolved_beam_reference']
        require(held['startup_mode'] == 'held' and held['frozen_lab_plane_injection'] and
                held['hold_transverse_until_activation'], 'Pending-particle release rules changed')
        require(manifest['resolved_initial_self_field']['model'] == 'cowan-relativistic-poisson' and
                manifest['resolved_initial_self_field']['enabled'], 'Fast reference model differs')
        require(card['mesh']['moving_window']['enabled'] and card['boundary']['cells'] == [10, 10, 0], 'Moving-window/CPML card differs')
        require(card['mesh']['cells'] == [160, 160, geometry['reference_cells_z']], 'Reference mesh differs')
        require(card['mesh']['duration'] == geometry['duration_box_s'], 'Hard duration differs')
        mesh = manifest['resolved_mesh']; cells = [160, 160, geometry['active_cells_z']]
        require(mesh['cells'] == cells and mesh['cell_size_m'] == card['mesh']['cell_size'], 'Resolved active mesh differs')
        center = geometry['active_origin_z_m'] + cells[2] * geometry['dz_m'] / 2
        require(np.allclose(mesh['center_m'], [0., 0., center], rtol=0, atol=2e-14), 'Resolved initial window position differs')
        require(not mesh['ballistic_outside_transverse'] and mesh['self_field_feedback'], 'Field feedback/transverse policy differs')
        require(card['stop'] == {'mode': 'after-last-element'}, 'Physical stop changed')
        require(not card.get('checkpoint') and not card['particle_retirement']['enabled'] and
                not card['detectors']['particle_planes'] and not card.get('trajectory', {}).get('enabled', False),
                'Unexpected checkpoint, retirement or particle output')
        field_cards = card['detectors']['field_planes']
        require(len(field_cards) == 1 and field_cards[0]['z'] == plan['field_plane_z_m'], 'Field plane differs')
        field_card = field_cards[0]
        require(not field_card.get('particle_background', {}).get('enabled', False), 'Field subtraction enabled')
        if plan['case'] == 'free_drift_0p5m':
            require(card['sources']['magnetic_elements'] == [] and card['sources']['incident_waves'] == [] and
                    field_card['z'] == .5, 'Free drift has prescribed fields or wrong plane')
        log = (root / 'aprl.log').read_text(errors='replace')
        require('Simulation failure' not in log and 'mesh.duration was exhausted' not in log, 'Run failed or exhausted duration')
        stop = re.findall(r'Configured stop reached: all (\d+) valid particles passed the last interaction boundary at z=([^ ]+) m', log)
        require(len(stop) == 1 and int(stop[0][0]) > 0 and
                np.isclose(float(stop[0][1]), field_card['z'], rtol=2e-6, atol=0), 'Missing expected physical stop')
        final = re.findall(r'\[resource\] phase=final step=(\d+)\b', log)
        require(len(final) == 1 and 'phase=solve-complete' in log and 'phase=mpi-finalize-complete' in log,
                'Normal solve/finalization not confirmed')
        final_step = int(final[0]); dt = geometry['dt_box_s']
        require(0 < final_step < geometry['reference_query_steps'], 'Final step outside planned duration')
        moving = re.findall(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)', log)
        require(len(moving) == 1 and int(moving[0][0]) == final_step, 'Missing moving-window completion')
        reference = Path(card['mesh']['moving_window']['reference_directory'])
        preparation = json.loads((reference / 'preparation.json').read_text())
        require(preparation['complete'] and preparation['particles'] == plan['records'] and
                preparation['interiorNodes'] == [139, 139, geometry['reference_cells_z'] - 1] and
                preparation['blockLayers'] == card['mesh']['moving_window']['block_layers'], 'Reference preparation incomplete or different')
        require('status complete' in (reference / 'metadata.txt').read_text(), 'Reference completion metadata missing')
        for key, filename in (('sourceFileBytes', 'rho.bin'), ('checkpointFileBytes', 'checkpoints.bin'), ('potentialFileBytes', 'phi.bin')):
            require((reference / filename).stat().st_size == preparation[key], 'Reference scalar file length differs: ' + filename)
        interval = card['energy_ledger']['sample_interval_steps']
        report['moving_window'] = check_window_csv(reference, plan['ranks'], final_step, interval, geometry)
        require(report['moving_window']['total_particle_exits'] == int(moving[0][1]) and
                report['moving_window']['total_cropped_carriers'] == int(moving[0][2]), 'CSV/log particle accounting differs')
        field_path = root / 'detectors' / (field_card['name'] + '.h5')
        require(sorted((root / 'detectors').glob('*.h5')) == [field_path], 'Unexpected/missing detector files')
        with h5py.File(field_path) as file:
            group = file['field_plane']; check_identity(group, manifest)
            require(string(group.attrs['frame']) == 'laboratory' and string(group.attrs['external_background_subtracted']) == 'none', 'Expected raw lab fields')
            require(string(group.attrs['particle_background_status']) == 'disabled; raw total field retained', 'Particle background subtraction differs')
            require('adjacent MPI half-cell planes included' in string(group.attrs['spatial_sampling']), 'Corrected MPI interpolation not identified')
            require('stored leapfrog half time' in string(group.attrs['magnetic_time_stagger']), 'Native magnetic staggering differs')
            require(float(group.attrs['plane_z_m']) == field_card['z'] and
                    int(group.attrs['nx']) == 160 and int(group.attrs['ny']) == 160 and
                    float(group.attrs['dx_m']) == card['mesh']['cell_size'][0] and
                    float(group.attrs['dy_m']) == card['mesh']['cell_size'][1], 'Field plane geometry differs')
            count = int(group['committed_samples'][()]); data = group['fields']; times = group['time_s'][:]
            require(int(group['complete'][()]) == 1 and count == final_step + 1 and
                    data.shape == (count, 160, 160) and times.shape == (count,), 'Incomplete field transaction/extent')
            require(data.dtype.names == ('electric_V_per_m', 'magnetic_T') and
                    all(data.dtype.fields[k][0].shape == (3,) for k in data.dtype.names), 'Unexpected E/B dataset schema')
            cadence = dt / geometry['gamma_boost']
            require(np.isclose(field_card['rhythm'], cadence, rtol=1e-12, atol=0), 'Field rhythm is not every full step')
            tolerance, span_tolerance = check_time_axis(times, cadence)
            first_time = geometry['source_frame_origin_lab_time_s'] + geometry['beta_boost'] * (field_card['z'] - geometry['release_lab_z_m']) / C
            require(abs(float(times[0]) - first_time) <= tolerance, 'Initial lab time differs from release frame')
            first_box = (field_card['z'] - geometry['release_lab_z_m']) / geometry['gamma_boost'] + geometry['translation_z_m']
            relative = first_box - geometry['active_origin_z_m'] + (1 - geometry['beta_boost']) * np.array([0., final_step * dt]) * C
            require(np.min(relative) >= 0 and np.max(relative) <= cells[2] * geometry['dz_m'], 'Detector left moving physical mesh')
            chunk = max(1, (32 * 1024 * 1024) // (160 * 160 * data.dtype.itemsize))
            for begin in range(0, count, chunk):
                finite(data[begin:min(count, begin + chunk)], 'field frames %d...' % begin)
            report['field'] = dict(file=str(field_path), complete=True, samples=count,
                first_lab_time_s=float(times[0]), last_lab_time_s=float(times[-1]), cadence_s=cadence,
                cadence_tolerance_s=tolerance, cumulative_time_tolerance_s=span_tolerance,
                finite_check_scope='every E/B scalar in every frame', magnetic_time_stagger=string(group.attrs['magnetic_time_stagger']))
        ledger_path = root / 'ledger' / card['energy_ledger']['filename']
        require(sorted((root / 'ledger').glob('*.h5')) == [ledger_path], 'Unexpected/missing ledger files')
        with h5py.File(ledger_path) as file:
            group = file['energy_ledger']; check_identity(group, manifest)
            records = group['records'][:]; count = int(group['committed_records'][()])
            require(int(group['complete'][()]) == 1 and count == len(records) and count > 0 and
                    int(group.attrs['mpi_size']) == plan['ranks'], 'Incomplete ledger transaction/layout')
            finite(records, 'ledger')
            expected_steps = list(range(0, final_step + 1, interval))
            if expected_steps[-1] != final_step: expected_steps.append(final_step)
            require(np.array_equal(records['step'], expected_steps), 'Ledger step coverage differs')
            require(np.allclose(records['time_box_s'], records['step'] * dt, rtol=2e-11, atol=0), 'Ledger box time differs')
            require(np.all(records['active_macroparticles'] + records['removed_macroparticles'] == plan['records']), 'Ledger particle inventory differs')
            prescribed = float(np.max(np.abs(records['prescribed_work_J'])))
            if plan['case'] == 'free_drift_0p5m': require(prescribed == 0., 'Free drift has nonzero prescribed-field work')
            report['ledger'] = dict(file=str(ledger_path), complete=True, records=count, final_step=final_step,
                prescribed_work_max_abs_J=prescribed, closure_valid_values=np.unique(records['closure_valid']).astype(int).tolist(),
                interpretation='Stationary-volume ledger is diagnostic only; moving inserted/cropped energy is reported separately.')
        if (root / 'EXIT_CODE').exists(): require((root / 'EXIT_CODE').read_text().strip() == '0', 'Recorded job exit was nonzero')
        report.update(pass_=True, output_transactions_complete=True, final_step=final_step,
            physical_stop=dict(valid_particles=int(stop[0][0]), logged_lab_z_m=float(stop[0][1])),
            reference_cells_z=geometry['reference_cells_z'], active_cells=cells,
            binary_sha256=plan['binary_sha256'], input_hash_scope='Full preflight SHA plus unchanged startup/completion stat',
            limitations=['Validates the configured stop, complete output transactions and finite stored fields; not pure radiation or numerical convergence.',
                'The full-input moving-window spectrum still requires comparison with the earlier full-window runs.',
                'After-last-element stops at the detector; particles need not all reach the downstream causal crop.',
                'Native half-step B is preserved; temporal E/B alignment belongs to subsequent spectral analysis.'])
        report['pass'] = True
        if plan['case'] == 'free_drift_0p5m':
            report['limitations'].append('No prescribed magnet/wave and exactly zero prescribed work do not certify zero radiation from startup or self-consistent motion.')
    except (ValueError, KeyError, OSError, TypeError, IndexError) as error:
        report['error'] = str(error) or type(error).__name__
    (root / 'completion-check.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    result = check(args.root.resolve())
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result['pass'] else 1)


if __name__ == '__main__':
    main()
