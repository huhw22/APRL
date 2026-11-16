"""Existing-data/code audit only: no particle or Maxwell evolution."""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import h5py
import numpy as np

LD = np.longdouble
C = LD(299792458)


def read(p):
    return json.loads(Path(p).read_text())


def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def identity(p):
    s = Path(p).stat()
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def cache_review(campaign, out):
    previous = campaign / 'momentum-ledger-20261116-r1'
    with np.load(previous / 'boost/cache.npz') as f:
        a = f['records']
    with np.load(previous / 'whole/energy-momentum/whole_cache.npz') as f:
        wj = np.r_[LD(0), np.cumsum(f['fixed_end_WJ_particle_J'], dtype=LD)]
        pj = np.r_[LD(0), np.cumsum(f['fixed_end_PJ_particle_J'], dtype=LD)]
        old_gap = f['particle_vs_deposited_current_Lorentz_gap_J']
    G, beta = LD(a['boost_gamma'][0]), LD(a['boost_beta'][0])
    external = a['prescribed_boris_work_J'].astype(LD)
    # Algebraic static-lab-magnet cancellation, not a newly measured force
    # trace: dKext' + beta*d(cPext') = 0. Do not call the inferred impulse an
    # independently measured electric/magnetic substep force decomposition.
    energy = G * (a['delta_particle_kinetic_box_J'].astype(LD) - external - wj)
    momentum = G * beta * (a['delta_particle_cPz_box_J'].astype(LD) + external / beta + pj)
    assert np.max(abs(energy + momentum - old_gap)) < LD('2e-18')
    source_energy = G * (a['particle_grid_work_J'].astype(LD) - wj)
    assert np.max(abs(source_energy - energy)) < LD('1e-15')
    states = np.genfromtxt(campaign / 'six-aprl/energy-momentum/worldtube/states.csv',
                           delimiter=',', names=True, dtype=LD)
    full = abs(a['released_charge_C'][:-1]) >= abs(LD(a['released_charge_C'][-1])) * (1 - LD('1e-12'))
    no_frozen = states['frozen_count'] == 0
    no_exit = a['removed_macroparticles'][1:] == 0
    clean = full & no_frozen & no_exit
    first_exit = int(np.flatnonzero(~no_exit)[0] + 1)
    full_step = int(np.flatnonzero(full & no_frozen)[0] + 1)
    phases = {}
    step = a['step'][1:]
    for name, mask in [('before_first_exit', step < first_exit),
                       ('injection_and_exit_overlap', (step >= first_exit) & (step < full_step)),
                       ('after_full_release', step >= full_step)]:
        phases[name] = dict(first_step=int(step[mask][0]), last_step=int(step[mask][-1]),
            delta_lab_gap_J=float(np.diff(old_gap)[mask].sum(dtype=LD)),
            delta_energy_pairing_gap_J=float(np.diff(energy)[mask].sum(dtype=LD)),
            delta_longitudinal_pairing_gap_J=float(np.diff(momentum)[mask].sum(dtype=LD)))
    with (out / 'energy_momentum_pairing.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['step', 'time_box_s', 'energy_pairing_after_static_external_cancel_J',
                         'longitudinal_pairing_after_static_external_cancel_J', 'lab_gap_J'])
        writer.writerows(zip(a['step'], a['time_box_s'], energy, momentum, old_gap))
    answer = dict(
        external_boost_work_projected_Gamma_J=float(G * external[-1]),
        energy_pairing_gap_after_static_external_cancel_J=float(energy[-1]),
        longitudinal_pairing_gap_after_static_external_cancel_J=float(momentum[-1]),
        total_Lorentz_gap_J=float(old_gap[-1]),
        independent_grid_work_vs_WJ_gap_J=float(source_energy[-1]),
        pairing_split_identity_max_J=float(np.max(abs(energy + momentum - old_gap))),
        full_release_and_zero_exit_common_steps=int(clean.sum()), phases=phases,
        static_external_impulse_inferred_algebraically_not_measured=True,
        missing_actual_particle_Ez_and_v_cross_B_force_trace=True)
    save(out / 'PAIRING.json', answer)
    return answer


class IndependentRays:
    """Direct inverse-Lorentz rays and Simpson quadrature, no polygon moments.

    Recorded currents are bilinear in t_box and moving-grid index. Split each
    ray at every time and spatial cell crossing. Its restriction is quadratic,
    so endpoint/midpoint Simpson is exact on each resulting segment.
    """
    def __init__(self, directory):
        m = json.loads((directory / 'METADATA.json').read_text(), parse_float=LD)
        self.m = m
        states = np.genfromtxt(directory / 'states.csv', delimiter=',', names=True, dtype=LD)
        self.t, self.o = states['time_box_s'], states['origin_z_box_m']
        self.g = LD(m['gamma'])
        self.b = np.sqrt((self.g - 1) * (self.g + 1)) / self.g
        self.eta = 1 / (self.g * (1 + self.b))
        self.q = np.memmap(directory / 'field_lines.bin', dtype='<f8', mode='r',
                           shape=(len(states), int(m['nz']), 6))
        self.dz = LD(m['dz_box_m'])
        self.bo = LD(m['frame_box_z_origin_m'])
        self.zo = LD(m['frame_lab_z_origin_m'])

    def integrate(self, kind, position, low, high):
        if kind == 'port':
            tlo = self.g * low + self.eta * (position - self.zo) / C
            thi = self.g * high + self.eta * (position - self.zo) / C
            intercept, slope = self.bo + (position - self.zo) / self.g, -self.b * C
            channel, jacobian = 1, 1 / self.g
        else:
            tau = position
            tlo = ((low - self.zo) / (self.g * C) + tau / self.eta) / (1 + self.b)
            thi = ((high - self.zo) / (self.g * C) + tau / self.eta) / (1 + self.b)
            intercept, slope = self.bo - C * tau / self.eta, C
            channel, jacobian = 0, self.g * C * (1 + self.b)
        if thi <= tlo:
            return LD(0), LD(0)
        assert tlo >= self.t[0] - LD('1e-25') and thi <= self.t[-1] + LD('1e-25')
        first = max(0, int(np.searchsorted(self.t, tlo, side='right')) - 1)
        last = min(len(self.t) - 2, int(np.searchsorted(self.t, thi)))
        total, coverage = LD(0), LD(0)
        for n in range(first, last + 1):
            t0, t1 = self.t[n:n+2]
            left, right = max(t0, tlo), min(t1, thi)
            if right <= left:
                continue
            odot = (self.o[n+1] - self.o[n]) / (t1 - t0)
            k0 = (intercept + slope * t0 - self.o[n]) / self.dz - LD('.5')
            kdot = (slope - odot) / self.dz
            kl, kr = k0 + kdot * (left - t0), k0 + kdot * (right - t0)
            cuts = [left, right]
            if kdot:
                for k in range(int(np.floor(min(kl, kr))) + 1, int(np.ceil(max(kl, kr)))):
                    time = t0 + (k - k0) / kdot
                    if left < time < right:
                        cuts.append(time)
            cuts.sort()
            for l, h in zip(cuts[:-1], cuts[1:]):
                vals = []
                for time in [l, (l + h) / 2, h]:
                    kt = k0 + kdot * (time - t0)
                    k = min(int(self.m['nz']) - 2, max(0, int(np.floor(kt))))
                    assert -LD('1e-10') <= kt <= LD(self.m['nz']) - 1 + LD('1e-10')
                    x, y = (time - t0) / (t1 - t0), kt - k
                    q00, q10, q01, q11 = [LD(v) for v in (
                        self.q[n,k,channel], self.q[n+1,k,channel],
                        self.q[n,k+1,channel], self.q[n+1,k+1,channel])]
                    vals.append((1-x) * ((1-y)*q00+y*q01) + x*((1-y)*q10+y*q11))
                total += (h-l) * (vals[0]+4*vals[1]+vals[2]) / 6 * jacobian
                coverage += (h-l) * jacobian
        expected = high-low if kind == 'port' else high-low
        assert abs(coverage / expected - 1) < LD('3e-12')
        return total, coverage


def plane_bounds(case, name, tau_ref):
    with h5py.File(case / 'detectors' / (name + '-particles-substeps.h5')) as f:
        g = f['particle_substep_plane']
        assert int(g['complete'][()]) == 1
        events = g['records'].fields(['time_s', 'crossing_direction', 'source_id'])[:]
        assert len(events) == 1733612 and np.all(events['crossing_direction'] == 1)
        assert len(np.unique(events['source_id'])) == len(events)
        tau = events['time_s'].astype(LD) - LD(g.attrs['plane_z_m']) / C - tau_ref
    return dict(minimum_tau_s=float(tau.min()), maximum_tau_s=float(tau.max()))


def quadrature_review(campaign, out):
    previous = campaign / 'momentum-ledger-20261116-r1'
    case = campaign / 'six-aprl/energy-momentum'
    physical = read(previous / 'lab/energy-momentum/LAB_RESULT.json')
    low, high = map(LD, physical['common_relative_tau_s'])
    # Config retains all long-double digits used by the original integrator.
    original = (previous / 'lab/energy-momentum/tube-config.txt').read_text().split()
    low, high = LD(original[7]), LD(original[8])
    zs = [LD(x) for x in original[10:]]
    surfaces = np.genfromtxt(previous / 'lab/energy-momentum/field_surfaces.csv',
                              delimiter=',', names=True, dtype=LD)
    rays = IndependentRays(case / 'worldtube')
    checks = []
    for i, z in enumerate(zs):
        port, coverage = rays.integrate('port', z, low, high)
        early = late = LD(0)
        if i:
            early, _ = rays.integrate('cap', low, zs[0], z)
            late, _ = rays.integrate('cap', high, zs[0], z)
        pairs = [('field_exit_signed_J', port), ('time_cap_early_signed_J', -early),
                 ('time_cap_late_signed_J', late)]
        for name, value in pairs:
            error = value - LD(surfaces[name][i])
            assert abs(error) < LD('1e-15'), (i, name, error)
            checks.append(dict(plane=i, channel=name, independently_integrated_J=float(value),
                               previous_J=float(surfaces[name][i]), difference_J=float(error)))
    save(out / 'INDEPENDENT_PORT_CAP_CHECK.json', dict(
        passed=True, method='Independent inverse Lorentz rays, cell crossings and exact Simpson quadrature',
        maximum_absolute_difference_J=max(abs(c['difference_J']) for c in checks), checks=checks,
        scope='Quadrature of already-recorded currents; does not validate native-to-continuum field reconstruction.'))
    tau_ref = LD(rays.m['frame_lab_t_origin_s']) - rays.zo / C
    bounds = {name: plane_bounds(case, name, tau_ref) for name in ['audit-entry', 'p1', 'p6']}
    save(out / 'PARTICLE_TIME_BOUNDS.json', bounds)
    integrator = previous / 'lab/worldtube_integrator'
    tests = [('empty_early', zs[1:], low, LD(bounds['p1']['minimum_tau_s']) - LD('2e-18')),
             ('empty_late', zs, LD(bounds['p6']['maximum_tau_s']) + LD('2e-18'), high)]
    vacuum = {}
    for label, planes, lo, hi in tests:
        if hi <= lo:
            vacuum[label] = dict(available=False)
            continue
        d = out / label
        d.mkdir()
        cfg = original[:]
        cfg[7:10] = [str(lo), str(hi), str(len(planes))]
        cfg = cfg[:10] + [str(z) for z in planes]
        (d / 'config.txt').write_text(' '.join(cfg) + '\n')
        subprocess.run([str(integrator), str(case / 'worldtube'), str(d/'config.txt'), str(d/'field_surfaces.csv')], check=True)
        data = np.genfromtxt(d / 'field_surfaces.csv', delimiter=',', names=True, dtype=LD)
        terms = [k for k in data.dtype.names if k.endswith('_J')]
        total = sum(data[k] for k in terms)
        vacuum[label] = dict(available=True, common_tau_s=[float(lo),float(hi)],
            no_particles_by_subluminal_monotone_tau_and_full_crossing_coverage=True,
            final_field_only_sum_J=float(total[-1]), maximum_abs_field_only_sum_J=float(max(abs(total))),
            final_field_components_J={k:float(data[k][-1]) for k in terms},
            caution='An empty-region nonzero sum tests recorded field discretization/reconstruction independently of local particle pushing; it is not a decomposition of the full-gate sum.')
    save(out / 'PARTICLE_FREE_FIELD_CHECK.json', vacuum)
    return vacuum


def main(campaign, out):
    campaign, out = Path(campaign).resolve(), Path(out).resolve()
    guard = read(out / 'PROTECTED_OUTPUTS.json')
    def unchanged():
        for path, before in guard.items():
            assert identity(path) == before, path
    unchanged()
    pairing = cache_review(campaign, out)
    vacuum = quadrature_review(campaign, out)
    unchanged()
    save(out / 'RESULT.json', dict(complete=True, existing_data_only=True,
        no_particle_or_Maxwell_evolution=True, pairing=pairing, particle_free_fields=vacuum,
        original_output_identities_unchanged=len(guard)))
    print('AUDIT_COMPLETE', flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])
