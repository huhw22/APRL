#!/usr/bin/env python3
"""Check the written derivation independently; no PIC, H5 or result rewrite.

Integrate analytic moments and Gaussian self-energy; then solve a separate
dimensional radial momentum ODE at several beam speeds and weak charges.
The latter distinguishes the fixed-plane coefficient 2/5 from the earlier
ultrarelativistic auxiliary expression (3 + beta**2)/10.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import c, e, epsilon_0, m_e
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import erf


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixed_plane_coefficient(strength, beta):
    """Independent units c=R0=m=1, omega^2=strength, cold uniform sphere."""
    nr, na = 20, 20
    radial_node, wr = np.polynomial.legendre.leggauss(nr)
    radius0 = (radial_node + 1) / 2
    # Integrate in r0 with density 3*r0^2. This integrates the leading radial
    # moments exactly, unlike a small quadrature in F with fractional powers.
    wr = wr / 2 * 3 * radius0**2
    q = radius0**3
    mu, wa = np.polynomial.legendre.leggauss(na)
    wa = wa / 2
    # Exact relativistic momentum equations; do not use the series being checked.
    def rhs(t, state):
        radius, momentum = state[:nr], state[nr:]
        return np.r_[momentum / np.sqrt(1 + momentum**2), strength * q / radius**2]

    sol = solve_ivp(rhs, [0, 8], np.r_[radius0, np.zeros(nr)],
                    method='DOP853', rtol=5e-13, atol=1e-15, dense_output=True)
    assert sol.success
    # Plane center times T0=3, T1=5; both planes are downstream of every cold start.
    means = []
    for center_time in (3., 5.):
        values = []
        for i in range(nr):
            for j, cosine in enumerate(mu):
                event = brentq(lambda t: cosine * sol.sol(t)[i] + beta * (t-center_time),
                               0, 8, xtol=2e-14)
                u = np.longdouble(sol.sol(event)[nr+i])
                gamma = np.sqrt(1 + u*u)
                # Stable gamma-1, and all angular cancellations in long double.
                energy_over_Gamma = u*u/(gamma+1) + np.longdouble(beta)*cosine*u
                values.append(np.longdouble(wr[i])*wa[j]*energy_over_Gamma)
        means.append(sum(values))
    return float((means[1]-means[0]) / (np.longdouble(strength)**2 * (25-9)))


def check(bundle):
    physical = bundle / 'physical_benchmarks/space_charge_scaling_20261030'
    weak_path = physical / 'WEAK_LIMIT_CHECK.json'
    case_path = physical / 'q0p001/case.json'
    before = {str(p.relative_to(bundle)): sha(p) for p in (weak_path, case_path)}
    moments = [quad(lambda r: 3*r**(p+2), 0, 1, epsabs=1e-13)[0] for p in (2, 4)]
    assert abs(moments[0]-3/5) < 1e-13 and abs(moments[1]-3/7) < 1e-13
    ratio = np.sqrt(moments[1]-moments[0]**2)/moments[0]
    assert abs(ratio-np.sqrt(4/21)) < 1e-13

    def enclosed(x):
        return erf(x/np.sqrt(2))-np.sqrt(2/np.pi)*x*np.exp(-x*x/2)

    assembly_integral = quad(lambda x: enclosed(x)*np.sqrt(2/np.pi)*x*np.exp(-x*x/2),
                             0, np.inf, epsabs=1e-13)[0]
    field_integral = .5*quad(lambda x: (enclosed(x)/x)**2,
                             0, np.inf, epsabs=1e-13)[0]
    assert abs(assembly_integral-1/(2*np.sqrt(np.pi))) < 1e-12
    assert abs(field_integral-assembly_integral) < 1e-12

    params = json.loads(case_path.read_text())['model']
    Gamma = params['beam_gamma_lab']; beta = np.sqrt(1-1/Gamma**2)
    scans = []
    for speed in (.6, .9, float(beta)):
        coefficients = [fixed_plane_coefficient(lam, speed) for lam in (1e-4, 1e-5, 1e-6)]
        extrapolated = (10*coefficients[2]-coefficients[1])/9
        assert abs(extrapolated-.4) < 5e-7
        assert abs(coefficients[2]-.4) < abs(coefficients[1]-.4) < abs(coefficients[0]-.4)
        scans.append({'beta': speed, 'omega_squared': [1e-4, 1e-5, 1e-6],
                      'measured_coefficients': coefficients, 'extrapolated_coefficient': extrapolated,
                      'derived_coefficient': .4, 'previous_auxiliary_coefficient': (3+speed**2)/10})

    weak = json.loads(weak_path.read_text())
    radius = np.sqrt(5)*params['sigma_m']
    omega2 = params['electrons']*e*e/(4*np.pi*epsilon_0*m_e*radius**3)
    corrected = .4*Gamma*m_e*omega2**2*radius**2*(weak['T1_s']**2-weak['T0_s']**2)/e
    old = next(x['cold_weak_limit_gain_eV'] for x in weak['rows'] if x['charge_factor']==.001)
    relative = corrected/old-1
    assert abs(relative-(4/(3+beta**2)-1)) < 1e-13
    assert all(sha(bundle/rel)==value for rel,value in before.items())
    return {'checks_pass': True, 'scope': 'Checks of the written physical derivation; not a PIC accuracy acceptance.',
            'uniform_radial_moments': moments, 'uniform_gain_std_over_mean': float(ratio),
            'gaussian_self_energy_assembly_integral': assembly_integral,
            'gaussian_self_energy_field_integral': field_integral,
            'fixed_plane_weak_limit_checks': scans,
            'current_case_auxiliary_correction': {'beam_gamma': Gamma,
                'previous_auxiliary_weak_gain_eV': old, 'corrected_weak_gain_eV': float(corrected),
                'absolute_difference_eV': float(corrected-old), 'relative_difference': float(relative)},
            'input_sha256_unchanged': before, 'new_pic_simulations': 0, 'h5_opened': False}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    result = check(a.bundle.resolve())
    with a.output.open('x') as f:
        f.write(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
