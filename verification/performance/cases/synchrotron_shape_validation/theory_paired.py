"""Independent paired finite-bend minus straight-worldline LW reference."""
import numpy as np
from scipy.constants import c, e, hbar
from scipy.ndimage import gaussian_filter1d
from scipy.special import roots_hermitenorm

from theory import (arc_spectrum, energy_density_spectrum,
                    gaussian_point_reference, uniform_velocity_fields)


def straight_velocity_spectrum(points, energy_eV, p, time_first=-.9e-15,
                               time_span=20e-15, ntime=8192,
                               sigma_time=0.):
    """Fourier transform of the incoming straight worldline on the PIC record."""
    dt = time_span / ntime
    u = time_first + np.arange(ntime) * dt
    margin = int(np.ceil(8 * sigma_time / dt)) if sigma_time else 0
    ue = time_first + np.arange(-margin, ntime + margin) * dt
    omega = np.asarray(energy_eV) * e / hbar
    phase = np.exp(1j * omega[:, None] * u[None, :]) * dt
    eout = np.zeros((len(omega), len(points), 3), complex)
    bout = np.zeros_like(eout)
    for start in range(0, len(points), 64):
        ef, bf = uniform_velocity_fields(points[start:start + 64], ue, p, False)
        if margin:
            scale = sigma_time / dt
            ef = gaussian_filter1d(ef, scale, axis=0, truncate=8,
                                   mode='nearest')[margin:margin + ntime]
            bf = gaussian_filter1d(bf, scale, axis=0, truncate=8,
                                   mode='nearest')[margin:margin + ntime]
        stop = start + ef.shape[1]
        eout[:, start:stop] = np.einsum('ft,tpd->fpd', phase, ef,
                                        optimize=True)
        bout[:, start:stop] = np.einsum('ft,tpd->fpd', phase, bf,
                                        optimize=True)
    return eout, bout


def gaussian_straight_reference(x, y, energies, p, sigma_xy=20e-6,
                                sigma_z=20e-9, order=9, **window):
    """Transverse Gaussian average of the matched straight control."""
    nodes, weights = roots_hermitenorm(order)
    weights /= np.sqrt(2 * np.pi)
    xx, yy = np.meshgrid(nodes * sigma_xy, nodes * sigma_xy, indexing='xy')
    quadrature = np.outer(weights, weights).ravel()
    points = np.c_[x - xx.ravel(), y - yy.ravel(),
                   np.full(xx.size, p['screen'])]
    ef, bf = straight_velocity_spectrum(
        points, energies, p, sigma_time=sigma_z / (p['beta'] * c), **window)
    return (np.einsum('fpd,p->fd', ef, quadrature),
            np.einsum('fpd,p->fd', bf, quadrature))


def gaussian_paired_reference(x, y, energies, p, sigma_xy=20e-6,
                              sigma_z=20e-9, order=9, **window):
    """Bend minus straight Fourier fields for an identical Gaussian bunch."""
    acceleration, acceleration_b, velocity, velocity_b = gaussian_point_reference(
        x, y, energies, p, sigma_xy=sigma_xy, sigma_z=sigma_z,
        order=order, **window)
    straight, straight_b = gaussian_straight_reference(
        x, y, energies, p, sigma_xy=sigma_xy, sigma_z=sigma_z,
        order=order, **window)
    return (acceleration + velocity - straight,
            acceleration_b + velocity_b - straight_b,
            acceleration, acceleration_b, velocity, velocity_b,
            straight, straight_b)


def self_check():
    """Small algebra/quadrature check independent of APRL output."""
    from theory import constants
    p = constants(length=.0024, screen=.3024)
    energy = np.array([8., 10., 12.])
    window = dict(time_first=-1e-15, time_span=5e-15, ntime=512)
    values = gaussian_paired_reference(0., 0., energy, p, order=4, **window)
    paired, paired_b = values[:2]
    reconstructed = values[2] + values[4] - values[6]
    reconstructed_b = values[3] + values[5] - values[7]
    scale = max(np.linalg.norm(paired), 1e-300)
    error = np.linalg.norm(paired - reconstructed) / scale
    scale_b = max(np.linalg.norm(paired_b), 1e-300)
    error_b = np.linalg.norm(paired_b - reconstructed_b) / scale_b
    if not (np.isfinite(paired).all() and np.isfinite(paired_b).all() and
            error < 1e-14 and error_b < 1e-14):
        raise RuntimeError('Paired theory self-check failed')
    density = energy_density_spectrum(paired, paired_b)
    return dict(passed=True, algebra_relative_error=float(error),
                magnetic_algebra_relative_error=float(error_b),
                finite_spectral_density=bool(np.isfinite(density).all()))


if __name__ == '__main__':
    import json
    print(json.dumps(self_check(), indent=2))
