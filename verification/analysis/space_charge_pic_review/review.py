#!/usr/bin/env python3
"""Read completed sphere PIC runs; independently check frames and momentum bias.

No PIC kernel is called. Stored detector rows and earlier reports are immutable.
Reference particles are sampled both at the actual PIC event time and at their
own reference crossing of the same plane. A momentum recentering is a diagnostic
only; the reported original PIC energies are never replaced by it.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-space-charge-review-mpl')
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import h5py
import numpy as np

SOURCE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SOURCE / 'performance/cases/space_charge_physics'))
from model import Sphere, G_BEAM, c, e, m_e
from scipy.optimize import brentq

LD = np.longdouble
MC2 = LD(m_e) * LD(c)**2 / LD(e)
G = LD(G_BEAM)
GB = np.sqrt((G-1)*(G+1))
BETA = GB/G


def read(path):
    with h5py.File(path, 'r') as h:
        group = h['particle_plane']
        assert group.attrs['frame'] in (b'laboratory', 'laboratory')
        assert int(group['complete'][()]) == 1
        arr = group['records'][:]
        assert len(arr) == int(group['committed_records'][()]) == 32768
    arr = arr[np.argsort(arr['source_id'])]
    assert np.array_equal(arr['source_id'], np.arange(1, len(arr)+1))
    assert np.all(arr['crossing_direction'] == 1)
    for name in arr.dtype.names:
        if arr[name].dtype.kind == 'f':
            assert np.isfinite(arr[name]).all()
    return arr


def gamma(u):
    return np.sqrt(1 + np.sum(u*u, axis=1))


def energy(u):
    u2 = np.sum(u*u, axis=1)
    return MC2*u2/(np.sqrt(1+u2)+1)


def inverse_boost(u):
    u = np.asarray(u, dtype=LD)
    # Independent direct long-double Lorentz transform.
    v = u.copy()
    v[:, 2] = G*u[:, 2] - GB*gamma(u)
    # Compare with the light-cone transform used by the existing analyzer.
    plus = (gamma(u)+u[:, 2])/(G+GB)
    light_z = (plus-(1+np.sum(u[:, :2]**2, axis=1))/plus)/2
    assert np.max(abs(v[:, 2]-light_z)) < LD('2e-12')
    return v, float(np.max(abs(v[:, 2]-light_z)))


def moments(x, w):
    avg = np.sum(w*x)
    std = np.sqrt(np.sum(w*(x-avg)**2))
    return float(avg), float(std)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def run(workspace, out):
    base = workspace/'results/energy_ledger_test_record_20261030'
    phys = base/'physical_benchmarks'
    paths = {k:phys/'space_charge_pic_20261030_r2'/k for k in ('uniform', 'converging')}
    paths['gaussian'] = phys/'space_charge_pic_gaussian_20261030_r4/gaussian'
    sources = {}
    table, summary = [], {}
    transform_max = identity_max = existing_report_max = 0.0
    for kind, d in paths.items():
        config = json.loads((d/'case.json').read_text())
        result = json.loads((d/'analysis/RESULT.json').read_text())
        assert result['execution_integrity_pass']
        params = config['model']
        assert params['beam_gamma_lab'] == G_BEAM
        model = Sphere(kind,params['electrons'],params['sigma_m'],params['nr'],params['nmu'],params['nphi'])
        with np.load(d/'reference_input/expected_events.npz') as h:
            expected = {k:h[k] for k in h.files}
        arr0 = read(d/'detectors/entry.h5')
        w = np.asarray(abs(arr0['charge_C']), dtype=LD)
        w /= np.sum(w)
        assert np.max(abs(w-model.weights)) < 1e-13
        offset = np.median(np.asarray(expected['entry__lab_time_s'], dtype=LD)-arr0['time_s'])
        u0, check0 = inverse_boost(arr0['proper_velocity'])
        k0 = energy(u0)
        uref0 = np.asarray(expected['entry__rest_u'], dtype=LD)
        kref0 = energy(uref0)
        lab0 = np.asarray(arr0['proper_velocity'], dtype=LD)
        labg0 = gamma(lab0)
        refg0 = gamma(np.asarray(expected['entry__lab_u'], dtype=LD))
        kind_rows = []
        for plane, earlier in zip(config['planes'], result['planes']):
            name = plane['name']
            path = d/'detectors'/f'{name}.h5'
            sources[str(path.relative_to(workspace))] = sha(path)
            arr = read(path)
            assert np.array_equal(arr['charge_C'], arr0['charge_C'])
            labu = np.asarray(arr['proper_velocity'], dtype=LD)
            u, check = inverse_boost(labu)
            transform_max = max(transform_max, check, check0)
            t = np.asarray(arr['time_s'], dtype=LD)+offset
            tau = (LD(model.entry_center_time)+G*(t-BETA*arr['position_m'][:, 2]/LD(c)))/LD(model.t0)
            r, p = model.shell_state(model.shell, np.asarray(tau, dtype=float))
            refu = np.asarray((np.sqrt(model.eta)*p)[:, None]*model.directions, dtype=LD)
            nominalu = np.asarray(expected[name+'__rest_u'], dtype=LD)
            gain = energy(u)-k0
            refgain = energy(refu)-kref0
            nominalgain = energy(nominalu)-kref0
            labgain = MC2*(gamma(labu)-labg0)
            reflabgain = G*(refgain+BETA*MC2*(refu[:, 2]-uref0[:, 2]))
            nominallabgain = MC2*(gamma(np.asarray(expected[name+'__lab_u'], dtype=LD))-refg0)
            delta_u = (u-u0)-(refu-uref0)
            avg_du = np.sum(w[:, None]*delta_u, axis=0)
            actual_pos = np.asarray(arr['position_m'], dtype=LD).copy()
            actual_pos[:, 2] = G*(actual_pos[:, 2]-BETA*LD(c)*t)
            reference_pos = np.asarray((r*model.sigma)[:, None]*model.directions, dtype=LD)
            position_error = actual_pos-reference_pos
            centroid_error = np.sum(w[:, None]*position_error, axis=0)
            centered_position_l2 = float(np.sqrt(np.sum(w*np.sum((position_error-centroid_error)**2,axis=1))
                /np.sum(w*np.sum(reference_pos**2,axis=1))))
            rest_error, _ = moments(gain-refgain, w)
            lab_error, _ = moments(labgain-reflabgain, w)
            kinetic_term = float(G)*rest_error
            momentum_term = float(G*BETA*MC2*avg_du[2])
            identity_error = abs(lab_error-kinetic_term-momentum_term)
            identity_max = max(identity_max, identity_error)
            assert identity_error < 0.001  # eV, vs 1e5 eV physical discrepancies
            # Independent gamma change must reproduce the earlier light-cone result.
            mean, std = moments(gain, w)
            labmean, labstd = moments(labgain, w)
            existing_report_max = max(existing_report_max, abs(mean-earlier['PIC_mean_rest_gain_eV']))
            assert abs(mean-earlier['PIC_mean_rest_gain_eV']) < 1e-5
            assert abs(labmean-earlier['PIC_mean_lab_gain_eV']) < 0.001
            # Remove just the measured common proper-momentum error for a diagnostic.
            # This is neither a changed simulation nor a corrected physics result.
            recentered_gain = energy(u-avg_du)-k0
            recentered_mean, recentered_std = moments(recentered_gain, w)
            error_rms = float(np.sqrt(np.sum(w*np.sum(delta_u**2, axis=1))))
            centered_rms = float(np.sqrt(np.sum(w*np.sum((delta_u-avg_du)**2, axis=1))))
            reference_rms = float(np.sqrt(np.sum(w*np.sum(refu**2, axis=1))))
            row = dict(case=kind, **earlier)
            row.update(mean_rest_gain_error_pct=100*(mean-earlier['reference_mean_rest_gain_eV'])/abs(earlier['reference_mean_rest_gain_eV']) if name!='entry' else 0,
                std_rest_gain_error_pct=100*(std/earlier['reference_std_rest_gain_eV']-1) if name!='entry' else 0,
                sigma_x_error_pct=100*(earlier['PIC_sigma_x_m']/earlier['reference_sigma_x_m']-1),
                sigma_z_error_pct=100*(earlier['PIC_sigma_z_rest_m']/earlier['reference_sigma_z_rest_m']-1),
                mean_lab_gain_error_eV=lab_error,
                rest_energy_error_contribution_lab_eV=kinetic_term,
                rest_longitudinal_momentum_error_contribution_lab_eV=momentum_term,
                gain_transform_identity_error_eV=identity_error,
                mean_delta_rest_ux_error=float(avg_du[0]),
                mean_delta_rest_uy_error=float(avg_du[1]),
                mean_delta_rest_uz_error=float(avg_du[2]),
                delta_rest_u_error_rms=error_rms,
                delta_rest_u_error_centered_rms=centered_rms,
                reference_rest_u_rms=reference_rms,
                proper_velocity_error_rms_fixed_v0_pct=100*error_rms/np.sqrt(model.eta),
                centered_proper_velocity_error_fixed_v0_pct=100*centered_rms/np.sqrt(model.eta),
                nominal_crossing_reference_mean_rest_gain_eV=moments(nominalgain,w)[0],
                nominal_crossing_reference_mean_lab_gain_eV=moments(nominallabgain,w)[0],
                mean_rest_event_time_shift_tau=moments(tau-expected[name+'__tau'],w)[0],
                maximum_absolute_rest_event_time_shift_tau=float(np.max(abs(tau-expected[name+'__tau']))),
                centroid_rest_x_error_m=float(centroid_error[0]),
                centroid_rest_y_error_m=float(centroid_error[1]),
                centroid_rest_z_error_m=float(centroid_error[2]),
                centered_position_relative_L2=centered_position_l2,
                diagnostic_recentered_mean_rest_gain_eV=recentered_mean,
                diagnostic_recentered_std_rest_gain_eV=recentered_std)
            table.append(row)
            kind_rows.append(row)
        n = json.loads((base/f'relative_ledger_20261030/physical_{kind}/normalization.json').read_text())
        normalization_path=base/f'relative_ledger_20261030/physical_{kind}/normalization.json'
        sources[str(normalization_path.relative_to(workspace))]=sha(normalization_path)
        summary[kind] = dict(directory=str(d.relative_to(workspace)),
            particles_per_plane=32768, number_of_planes=6, steps=result['steps'],
            model=config['model'], final=kind_rows[-1],
            maximum_absolute_mean_rest_gain_error_pct=max(abs(r['mean_rest_gain_error_pct']) for r in kind_rows[1:]),
            maximum_absolute_std_rest_gain_error_pct=max(abs(r['std_rest_gain_error_pct']) for r in kind_rows[1:]),
            maximum_absolute_sigma_error_pct=max(max(abs(r['sigma_x_error_pct']),abs(r['sigma_z_error_pct'])) for r in kind_rows),
            sampled_minimum_sigma_x_plane=min(kind_rows,key=lambda r:r['PIC_sigma_x_m'])['plane'],
            ledger=n)
        if kind == 'converging':
            radial_weights = np.asarray(model.weights.reshape(model.nr,-1).sum(axis=1))
            def envelope_derivative(tau):
                state = model.solution.sol(tau)
                x,p = state[:model.nr],state[model.nr:]
                return np.sum(radial_weights*x*p/np.sqrt(1+model.eta*p*p))
            turn_tau = brentq(envelope_derivative,0,3)
            xturn = model.solution.sol(turn_tau)[:model.nr]
            summary[kind]['simultaneous_reference_rms_turn_tau'] = turn_tau
            summary[kind]['simultaneous_reference_rms_minimum_m'] = float(model.sigma*np.sqrt(np.sum(radial_weights*xturn*xturn)/3))
        for f in (d/'case.json', d/'analysis/RESULT.json', d/'analysis/comparison.csv', d/'reference_input/expected_events.npz'):
            sources[str(f.relative_to(workspace))] = sha(f)
    out.mkdir(parents=True, exist_ok=False)
    with (out/'comparison_all_planes.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=table[0].keys())
        writer.writeheader(); writer.writerows(table)
    (out/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    checks = dict(status='ANALYSIS_VALIDATED_PHYSICS_HAS_UNRESOLVED_DISCREPANCIES',
        completed_particle_planes=18, records_per_plane=32768,
        inverse_boost_independent_max_u_difference=transform_max,
        lab_gain_transform_max_identity_error_eV=identity_max,
        independent_rest_mean_gain_max_difference_from_previous_eV=existing_report_max,
        raw_file_hashes_unchanged=all(sha(workspace/p)==v for p,v in sources.items()),
        physics_methods_modified=False, simulations_submitted=False)
    (out/'VALIDATION.json').write_text(json.dumps(checks,indent=2)+'\n')
    (out/'INPUT_SHA256.json').write_text(json.dumps(sources,indent=2)+'\n')
    print(json.dumps(checks,indent=2))
    for k,s in summary.items():
        r=s['final']
        print(k, {n:r[n] for n in ('PIC_mean_rest_gain_eV','reference_mean_rest_gain_eV','std_rest_gain_error_pct','sigma_x_error_pct','sigma_z_error_pct','mean_lab_gain_error_eV','mean_delta_rest_uz_error','rest_energy_error_contribution_lab_eV','rest_longitudinal_momentum_error_contribution_lab_eV','nominal_crossing_reference_mean_lab_gain_eV','diagnostic_recentered_std_rest_gain_eV','proper_velocity_error_rms_fixed_v0_pct')})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.workspace.resolve(), args.output.resolve())
