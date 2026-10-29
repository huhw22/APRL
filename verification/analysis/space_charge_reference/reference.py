#!/usr/bin/env python3
"""Cold spherical Coulomb expansion references; no APRL simulation is run.

Nonrelativistic continuum model in the bunch mean rest frame. Gaussian
trajectories are valid only before shell crossing. The uniform sphere uses
the exact self-similar solution of that model. All output is labeled theory.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.constants import c, epsilon_0, e, m_e
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erf

os.environ.setdefault('MPLCONFIGDIR', '/tmp/aprl-space-charge-reference-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KCOULOMB = 1 / (4 * np.pi * epsilon_0)
SQRT2PI = np.sqrt(2 / np.pi)


def gaussian_cdf(x):
    x = np.asarray(x, dtype=float)
    result = erf(x / np.sqrt(2)) - SQRT2PI * x * np.exp(-x*x/2)
    # Avoid cancellation near the origin.
    small = np.abs(x) < .05
    series = np.zeros_like(x)
    import math
    for k in range(9):
        series += (-1)**k * x**(2*k+3) / (2**k * math.factorial(k) * (2*k+3))
    return np.where(small, SQRT2PI*series, result)


def gaussian_pdf(x):
    return SQRT2PI*x*x*np.exp(-x*x/2)


def scale_solution(b):
    """Invert y*sqrt(1+y^2)+asinh(y)=b, with a=1+y^2."""
    b = np.asarray(b, dtype=float)
    y = np.where(b < 1, b/2, np.sqrt(b))
    for _ in range(24):
        delta = (y*np.sqrt(1+y*y)+np.arcsinh(y)-b)/(2*np.sqrt(1+y*y))
        y = np.maximum(0, y-delta)
        if np.max(np.abs(delta)) < 1e-14:
            break
    assert np.max(np.abs(y*np.sqrt(1+y*y)+np.arcsinh(y)-b)) < 1e-11
    return 1+y*y, y


def gaussian_shells(x, tau):
    """r/sigma, v/(sigma/t0), dr/dr0; tau=t/t0."""
    f = gaussian_cdf(x)
    b = tau*np.sqrt(2*f/x**3)
    a, y = scale_solution(b)
    speed = np.sqrt(2*f/x)*y/np.sqrt(a)
    jac = a + x*y/np.sqrt(a)*b/2*(gaussian_pdf(x)/f-3/x)
    return a*x, speed, jac


def minimum_jacobian(tau):
    x = np.linspace(.01, 20, 2000)
    jac = gaussian_shells(x, tau)[2]
    i = int(np.argmin(jac))
    lo, hi = x[max(0, i-1)], x[min(len(x)-1, i+1)]
    result = minimize_scalar(lambda z: float(gaussian_shells(z, tau)[2]),
                             bounds=(lo, hi), method='bounded')
    return float(result.fun), float(result.x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--electrons', type=float, default=1e5)
    parser.add_argument('--sigma-m', type=float, default=1e-4)
    parser.add_argument('--refresh', action='store_true', help='Refresh this theory-only output; never an APRL result directory')
    args = parser.parse_args()
    out = args.output.resolve()
    if args.refresh:
        assert json.loads((out/'reference.json').read_text())['status']=='THEORY_REFERENCE_ONLY_NOT_APRL_PHYSICS_ACCEPTANCE'
    else:
        out.mkdir(parents=True, exist_ok=False)
    (out/'figures').mkdir(exist_ok=args.refresh)
    n, sigma = args.electrons, args.sigma_m
    assert n > 0 and sigma > 0
    charge = n*e
    energy = KCOULOMB*n*e*e/sigma  # Per physical electron, not per macro-particle.
    t0 = np.sqrt(m_e*sigma**3/(KCOULOMB*n*e*e))
    nr_strength = energy/(m_e*c*c)
    uniform_radius = np.sqrt(5.)  # Equal initial rms sigma_x in both cases.
    ug0, uu0 = 1/(2*np.sqrt(np.pi)), 3/(5*uniform_radius)
    nodes, weights = np.polynomial.legendre.leggauss(600)
    x, wx = 5*(nodes+1), 5*weights
    f, df = gaussian_cdf(x), gaussian_pdf(x)
    mass_weights = wx*df
    assert abs(np.sum(mass_weights)-1) < 1e-12
    assert abs(np.sum(mass_weights*x*x)/3-1) < 1e-12
    rows = []
    for tau in np.linspace(0, 3, 121):
        radius, speed, jac = gaussian_shells(x, tau)
        kinetic = .5*speed*speed  # In units of per-electron energy.
        mean = np.sum(mass_weights*kinetic)
        std = np.sqrt(np.sum(mass_weights*(kinetic-mean)**2))
        u = np.sum(mass_weights*f/radius)
        # Independent field-energy integral, including analytic exterior tail.
        uradial = .5*np.sum(wx*f*f/radius**2*jac)+.5*float(gaussian_cdf(10))**2/float(gaussian_shells(10., tau)[0])
        a = float(scale_solution(tau*np.sqrt(2/uniform_radius**3))[0])
        uniform_kmax = (1-1/a)/uniform_radius
        assert abs(mean+u-ug0) < 2e-13
        assert abs(uradial-u) < 2e-13
        assert np.min(jac) > 0
        rows.append({
            'tau_t_over_t0': float(tau), 'rest_time_s': float(tau*t0),
            'gaussian_sigma_x_over_initial': float(np.sqrt(np.sum(mass_weights*radius*radius)/3)),
            'gaussian_mean_gain_eV': float(mean*energy/e),
            'gaussian_std_gain_eV': float(std*energy/e),
            'gaussian_kinetic_J': float(n*energy*mean),
            'gaussian_field_J': float(n*energy*u),
            'gaussian_field_energy_fraction': float(u/ug0),
            'gaussian_field_integral_check_J': float(n*energy*uradial),
            'uniform_scale_a': a,
            'uniform_mean_gain_eV': float(.6*uniform_kmax*energy/e),
            'uniform_std_gain_eV': float(np.sqrt(12/175)*uniform_kmax*energy/e),
            'uniform_kinetic_J': float(n*energy*uu0*(1-1/a)),
            'uniform_field_J': float(n*energy*uu0/a),
            'uniform_field_energy_fraction': 1/a,
        })

    # A separate time-domain ODE checks the implicit shell solution.
    probe_x = np.linspace(.15, 6, 64)
    rhs_f = gaussian_cdf(probe_x)
    def rhs(t, state):
        return np.r_[state[64:], rhs_f/state[:64]**2]
    sol = solve_ivp(rhs, [0, 3], np.r_[probe_x, np.zeros(64)],
                    rtol=2e-12, atol=2e-13, method='DOP853', t_eval=[.5, 1, 2, 3])
    assert sol.success
    ode_difference = 0.
    for k, tau in enumerate(sol.t):
        analytic = gaussian_shells(probe_x, tau)[0]
        ode_difference = max(ode_difference, float(np.max(np.abs(sol.y[:64, k]/analytic-1))))
    assert ode_difference < 1e-10
    poisson_u0 = quad(lambda z: .5*float(gaussian_cdf(z))**2/z**2, 0, np.inf,
                     epsabs=1e-12, epsrel=1e-12)[0]
    assert abs(poisson_u0-ug0) < 1e-12
    short_r = gaussian_shells(x, 1e-3)[0]
    virial = (np.sum(mass_weights*short_r**2)/3-1)/1e-6
    assert abs(virial-ug0/3) < 1e-7
    crossing = brentq(lambda t: minimum_jacobian(t)[0], 3, 30)
    min_j, min_x = minimum_jacobian(3)

    # A uniform sphere also stays self-similar with a coherent initial chirp.
    # Its independently known turning radius and energy invariant check the ODE.
    chirp = -.2  # h0*t0, an illustrative reference parameter, not production.
    def turn_event(t, state):
        return state[1]
    turn_event.direction = 1
    chirp_sol = solve_ivp(lambda t, s: [s[1], 1/(uniform_radius**3*s[0]**2)],
                          [0, 6], [1, chirp], rtol=2e-12, atol=2e-13,
                          method='DOP853', dense_output=True, events=turn_event)
    assert chirp_sol.success and len(chirp_sol.t_events[0])==1
    turning_tau = float(chirp_sol.t_events[0][0])
    turning_scale = float(chirp_sol.sol(turning_tau)[0])
    turning_exact = 1/(1+.5*chirp**2*uniform_radius**3)
    assert abs(turning_scale-turning_exact) < 1e-10
    chirp_rows=[]
    initial_k = .3*uniform_radius**2*chirp**2
    max_chirp_energy_error=0.
    for tau in np.linspace(0, 6, 241):
        a, adot = chirp_sol.sol(tau)
        kinetic, field = .3*uniform_radius**2*adot**2, uu0/a
        max_chirp_energy_error=max(max_chirp_energy_error,abs(kinetic+field-initial_k-uu0))
        chirp_rows.append({'tau_t_over_t0':float(tau),'rest_time_s':float(tau*t0),
           'scale_a':float(a),'kinetic_J':float(n*energy*kinetic),
           'field_J':float(n*energy*field),'delta_kinetic_J':float(n*energy*(kinetic-initial_k)),
           'delta_field_J':float(n*energy*(field-uu0)),
           'mean_particle_energy_change_eV':float(energy/e*(kinetic-initial_k))})
    assert max_chirp_energy_error < 1e-11
    with (out/'converging_uniform.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=chirp_rows[0].keys())
        writer.writeheader();writer.writerows(chirp_rows)

    with (out/'reference.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    selected = [rows[i] for i in (0, 20, 40, 80, 120)]
    summary = {
        'status': 'THEORY_REFERENCE_ONLY_NOT_APRL_PHYSICS_ACCEPTANCE',
        'model': 'Cold isolated spherical continuum; primary Gaussian and uniform cases initially at rest; separate uniform inward-chirp case. Nonrelativistic mean field before shell crossing.',
        'frame': 'bunch mean rest frame, common time',
        'example_parameters_not_a_production_card': {'electrons': n, 'charge_C': charge, 'sigma_x_y_z_m': sigma, 'uniform_radius_m': uniform_radius*sigma},
        't0_s': t0, 'energy_unit_per_electron_J': energy,
        'nonrelativistic_strength': nr_strength,
        'gaussian_initial_full_space_field_J': n*energy*ug0,
        'uniform_initial_full_space_field_J': n*energy*uu0,
        'gaussian_max_total_gain_per_electron_eV_if_fully_dispersed': energy*ug0/e,
        'uniform_energy_gain_std_over_mean_nonrelativistic': np.sqrt(4/21),
        'gaussian_shell_crossing_estimate_tau': crossing,
        'gaussian_shell_crossing_estimate_s': crossing*t0,
        'sampled_time_interval_tau': [0, 3],
        'minimum_gaussian_mapping_jacobian_at_final_time': min_j,
        'jacobian_minimum_initial_r_over_sigma': min_x,
        'gaussian_field_fraction_outside_8sigma_initial': float(quad(lambda z: .5*float(gaussian_cdf(z))**2/z**2, 8, np.inf, epsabs=1e-13)[0]/ug0),
        'independent_reference_checks': {'maximum_implicit_vs_ODE_relative_radius_difference': ode_difference,
             'gaussian_initial_energy_Poisson_integral_dimensionless': poisson_u0,
             'short_time_rms_squared_coefficient_measured': virial,
             'short_time_rms_squared_coefficient_analytic': ug0/3},
        'converging_uniform_sphere': {
            'initial_chirp_h_times_t0':chirp,'turning_tau':turning_tau,
            'turning_rest_time_s':turning_tau*t0,'minimum_scale_a_numerical':turning_scale,
            'minimum_scale_a_analytic':turning_exact,'initial_kinetic_J':n*energy*initial_k,
            'minimum_mean_energy_change_eV':-energy/e*initial_k,
            'max_energy_invariant_error_in_N_times_energy_units':max_chirp_energy_error,
            'physical_scope':'Self-similar coherent contraction then expansion in the nonrelativistic continuum model.'},
        'selected_predictions': selected,
        'new_APRL_jobs': 0,
    }
    (out/'reference.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), layout='constrained')
    tau = [r['tau_t_over_t0'] for r in rows]
    for prefix, color, line, label in [('gaussian', '#4b78a4', '-', 'Cold spherical Gaussian'),
                                      ('uniform', '#d38b39', '--', 'Cold uniform sphere')]:
        size = 'gaussian_sigma_x_over_initial' if prefix=='gaussian' else 'uniform_scale_a'
        axes[0].plot(tau, [r[size] for r in rows], line, color=color, label=label)
        axes[1].plot(tau, [r[prefix+'_mean_gain_eV'] for r in rows], line, color=color, label=label)
        axes[2].plot(tau, [r[prefix+'_field_energy_fraction'] for r in rows], line, color=color, label=label)
    titles = ['Rms expansion', 'Mean kinetic-energy gain', 'Remaining electrostatic energy']
    units = [r'$\sigma_x(t)/\sigma_x(0)$', 'eV per electron', r'$U(t)/U(0)$']
    for ax, title, unit in zip(axes, titles, units):
        ax.set_title(title); ax.set_ylabel(unit); ax.set_xlabel(r'$t/t_0$'); ax.grid(alpha=.2)
        ax.set_xlim(0, 3); ax.set_xticks([0, 1, 2, 3])
    axes[0].legend(fontsize=8)
    fig.suptitle(f'Theory only; no APRL data. N={n:g}, initial rms={sigma*1e6:g} um, t0={t0*1e12:.1f} ps.', fontsize=11)
    fig.savefig(out/'figures/expansion_energy_reference.png', dpi=180)
    fig.savefig(out/'figures/expansion_energy_reference.pdf'); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
    chirp_tau=[r['tau_t_over_t0'] for r in chirp_rows]
    axes[0].plot(chirp_tau,[r['scale_a'] for r in chirp_rows],color='#4b78a4')
    axes[0].axhline(turning_exact,color='#777777',linestyle=':',label=f'Analytic minimum {turning_exact:.4f}')
    axes[0].set_ylabel(r'$R(t)/R(0)$');axes[0].set_title('Coherent contraction and re-expansion');axes[0].legend(fontsize=8)
    axes[1].plot(chirp_tau,[r['delta_kinetic_J']/1e-15 for r in chirp_rows],color='#4b78a4',label='Particle kinetic-energy change')
    axes[1].plot(chirp_tau,[r['delta_field_J']/1e-15 for r in chirp_rows],color='#d38b39',linestyle='--',label='Electrostatic-energy change')
    axes[1].set_ylabel('Energy change (fJ)');axes[1].set_title('Energy transferred between particles and field');axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xlabel(r'$t/t_0$');ax.set_xlim(0,6);ax.set_xticks([0,2,4,6]);ax.grid(alpha=.2)
        ax.axvline(turning_tau,color='#888888',linestyle=':',linewidth=1)
    fig.suptitle(f'Theory only; no APRL data. N={n:g}, initial rms={sigma*1e6:g} um, h0*t0={chirp:g}.',fontsize=10)
    fig.savefig(out/'figures/converging_uniform_reference.png',dpi=180)
    fig.savefig(out/'figures/converging_uniform_reference.pdf');plt.close(fig)
    (out/'SOURCE_SHA256.json').write_text(json.dumps({'reference.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2)+'\n')
    print(json.dumps({k: summary[k] for k in ['status', 't0_s', 'gaussian_initial_full_space_field_J',
        'gaussian_shell_crossing_estimate_tau', 'independent_reference_checks']}, indent=2))


if __name__ == '__main__':
    main()
