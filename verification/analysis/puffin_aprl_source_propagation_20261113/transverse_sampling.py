"""Collocation and transverse sampling audit on a closed Puffin field cache.

For four complex corner values z_i, the exact identity
    mean(|z_i|^2) - |mean(z_i)|^2 = sum_{i<j}|z_i-z_j|^2/16
splits into amplitude and phase terms via
    |z_i-z_j|^2 = (|z_i|-|z_j|)^2
                  + 2*(|z_i||z_j|-Re(z_i*conj(z_j))).
These are interpolation diagnostics, not corrections to a simulated energy.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(value, message):
    if not value:
        raise RuntimeError(message)


def main(root):
    root = Path(root).resolve()
    workspace = root.parents[1]
    require(read(root / 'physical-route-control-r1/COMPLETE.json')['complete'], 'Parent incomplete')
    out = root / 'sampling-audit-r1'
    out.mkdir(exist_ok=True)
    cache = root / 'physical-route-control-r1/analysis/fields.npz'
    with np.load(cache) as f:
        frequency, x, y, field = f['frequency_Hz'], f['x_m'], f['y_m'], f['Puffin']
    require(np.isfinite(field).all(), 'Nonfinite closed cache')
    dx, dy = x[1] - x[0], y[1] - y[0]
    ix = np.flatnonzero(abs(x) <= .00056 + 1e-12)
    iy = np.flatnonzero(abs(y) <= .00056 + 1e-12)
    require(max(abs(x[ix[0]] + .00056), abs(x[ix[-1]] - .00056),
                abs(y[iy[0]] + .00056), abs(y[iy[-1]] - .00056)) < 1e-12,
            'Aperture must coincide with native cell edges')
    nodes = field[:, iy][:, :, ix]
    corners = [nodes[:, :-1, :-1], nodes[:, :-1, 1:], nodes[:, 1:, :-1], nodes[:, 1:, 1:]]
    reference = sum(abs(z) ** 2 for z in corners) / 4
    average = sum(corners) / 4
    amplitude = np.zeros_like(reference)
    phase = np.zeros_like(reference)
    for i in range(4):
        for j in range(i + 1, 4):
            ai, aj = abs(corners[i]), abs(corners[j])
            amplitude += (ai - aj) ** 2 / 16
            phase += (ai * aj - (corners[i] * corners[j].conj()).real) / 8

    def integrate(values):
        return float(simpson(values.sum(axis=(1, 2, 3)) * dx * dy, x=frequency))

    native = integrate(reference)
    bilinear = integrate(abs(average) ** 2)
    amp_loss, phase_loss = integrate(amplitude), integrate(phase)
    closure = native - bilinear - amp_loss - phase_loss
    require(abs(closure) < 1e-13 * native, 'Complex pairwise variance identity failed')
    stages = read(root / 'physical-route-control-r1/analysis/MEASUREMENT_RESULT.json')['stages']
    energies = {row['stage']: row['band_energy_J'] for row in stages}
    require(abs(native / energies['Puffin native'] - 1) < 1e-12, 'Native quadrature mismatch')
    require(abs(bilinear / energies['Both-axis linear sampling'] - 1) < 1e-12,
            'Four-node and Fourier linear sampling disagree')
    # Pure coordinate registration on the complete periodic grid changes only
    # each transverse Fourier mode's phase, not its norm.
    kx = 2 * np.pi * np.fft.fftfreq(len(x), dx)
    ky = 2 * np.pi * np.fft.fftfreq(len(y), dy)
    shift = np.exp(1j * (ky[:, None] * dy / 2 + kx[None, :] * dx / 2))
    shifted = np.fft.ifft2(np.fft.fft2(field, axes=(1, 2)) * shift[None, :, :, None], axes=(1, 2))
    full0, full1 = integrate(abs(field) ** 2), integrate(abs(shifted) ** 2)
    require(abs(full1 / full0 - 1) < 1e-12, 'Pure periodic Fourier shift changed total norm')

    # These are the actual frozen sources used for the existing Puffin binary.
    puf = workspace / 'results/puffin_same_sample_20261113_r1/source'
    frozen = read(workspace / 'results/puffin_same_sample_20261113_r3/FROZEN.json')
    code_hashes = {}
    for name in ['puffin/lib/interaction/system_interpolation.f90',
                 'puffin/lib/interaction/rhs.f90', 'puffin/lib/io/hdf5_puff_coll.f90']:
        path = puf / name
        code_hashes[str(path)] = sha(path)
        require(code_hashes[str(path)] == frozen[str(path)], 'Frozen Puffin source changed')
    result = dict(complete=True, no_new_simulations=True, no_new_seeds=True,
        field_cache=str(cache), field_cache_sha256=sha(cache), source_hashes=code_hashes,
        Puffin_transverse_components_collocated=True,
        basis='Same sAr/sAi node indices and lis_GR weights for gather and deposition; nodal output.',
        native_band_energy_J=native, bilinear_band_energy_J=bilinear,
        total_interpolation_loss_J=native - bilinear,
        amplitude_variation_loss_J=amp_loss, phase_variation_loss_J=phase_loss,
        phase_fraction_of_interpolation_loss=phase_loss / (amp_loss + phase_loss),
        total_loss_fraction_of_native=(native - bilinear) / native,
        amplitude_loss_fraction_of_native=amp_loss / native,
        phase_loss_fraction_of_native=phase_loss / native,
        variance_identity_residual_J=closure,
        pure_registration_finite_aperture_change_fraction=energies['Fourier positions only'] / native - 1,
        pure_periodic_registration_total_norm_relative=full1 / full0 - 1,
        observed_stages=[dict(stage=row['stage'], band_energy_J=row['band_energy_J'],
                             change_fraction_from_native=row['band_energy_J'] / native - 1)
                         for row in stages if row['stage'] != 'PIC stored detector'],
        limitations=[
            'Phase/amplitude decomposition is an exact four-corner interpolation identity on this Puffin cache.',
            'It is not the missing native PIC energy and not a decomposition of the solvers physical output difference.',
            'Finite-aperture Fourier registration assumes the stored periodic band-limited representation.',
            'Longitudinal/time interpolation and feedback are outside this audit.'])
    (out / 'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    theory = []
    for delta in [.1, .2, .4, .6, .8, 1., 1.2]:
        theory.append(dict(neighbor_phase_step_rad=delta, one_axis_midpoint_power_ratio=float(np.cos(delta / 2) ** 2),
                           relative_power_loss=float(np.sin(delta / 2) ** 2)))
    with (out / 'plane_wave_estimate.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(theory[0]))
        writer.writeheader()
        writer.writerows(theory)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main(sys.argv[1])
