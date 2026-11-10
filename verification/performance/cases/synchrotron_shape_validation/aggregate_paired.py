"""Aggregate four closed bend-minus-straight spectra against paired LW theory."""
import csv
import json
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, h, hbar, mu_0
from scipy.special import roots_hermitenorm

from prepare_paired import require, save, sha
from theory import arc_spectrum, energy_density_spectrum, velocity_spectrum
from theory_paired import (gaussian_paired_reference,
                           straight_velocity_spectrum, self_check)


def load_field(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    with h5py.File(path, 'r') as handle:
        group = handle['field_plane']
        require(int(group['complete'][()]) == 1, 'Incomplete field plane')
        frame = group.attrs['frame']
        frame = frame.decode() if isinstance(frame, bytes) else frame
        require(frame == 'laboratory', 'Wrong field frame')
        n = int(group['committed_samples'][()])
        require(n == len(group['time_s']) == len(group['fields']) and n > 32,
                'Uncommitted field tail')
        times = group['time_s'][:n]
        dt = float((times[-1] - times[0]) / (n - 1))
        require(dt > 0, 'Nonpositive sample spacing')
        jitter = float(np.max(abs(times - times[0] - np.arange(n) * dt)) / dt)
        require(jitter < 1e-4, 'Nonuniform timestamps')
        nx, ny = int(group.attrs['nx']), int(group.attrs['ny'])
        dx, dy = float(group.attrs['dx_m']), float(group.attrs['dy_m'])
        x = float(group.attrs['x_first_m']) + np.arange(nx) * dx
        y = float(group.attrs['y_first_m']) + np.arange(ny) * dy
        ix, iy = np.flatnonzero(abs(x) < 50e-6), np.flatnonzero(abs(y) < 50e-6)
        require(len(ix) and len(iy) and np.all(np.diff(ix) == 1) and
                np.all(np.diff(iy) == 1), 'Invalid analysis aperture')
        cube = np.empty((n, len(iy), len(ix), 6))
        for start in range(0, n, 16):
            stop = min(n, start + 16)
            records = group['fields'][start:stop]
            electric, magnetic = (records['electric_V_per_m'],
                                  records['magnetic_T'])
            require(np.isfinite(electric).all() and np.isfinite(magnetic).all(),
                    'Nonfinite field record')
            cube[start:stop, :, :, :3] = electric[:, iy[0]:iy[-1] + 1,
                                                       ix[0]:ix[-1] + 1]
            cube[start:stop, :, :, 3:] = magnetic[:, iy[0]:iy[-1] + 1,
                                                    ix[0]:ix[-1] + 1]
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            'Closed field file changed')
    return dict(times=times, dt=dt, jitter=jitter, x=x[ix], y=y[iy],
                dx=dx, dy=dy, cube=cube)


def transform(data, first):
    n, dt = len(data['times']), data['dt']
    energy = h * np.fft.rfftfreq(n, dt) / e
    phase = np.exp(1j * energy * e / hbar * first)[:, None, None, None]
    fourier = np.conj(np.fft.rfft(data['cube'], axis=0)) * dt * phase
    density = energy_density_spectrum(fourier[..., :3], fourier[..., 3:])
    density[0] *= .5
    if n % 2 == 0:
        density[-1] *= .5
    return energy, fourier, density


def parseval(cube, fourier, density, dt, de, dx, dy):
    time_energy = float(np.sum(cube[..., 0] * cube[..., 4] -
                               cube[..., 1] * cube[..., 3]) * dt * dx * dy / mu_0)
    fft_energy = float(np.sum(density) * de * dx * dy)
    relative = abs(fft_energy - time_energy) / max(abs(time_energy),
                                                   abs(fft_energy), 1e-300)
    return dict(time_domain_J=time_energy, frequency_domain_J=fft_energy,
                relative_error=float(relative))


def load_particles(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    with h5py.File(path, 'r') as handle:
        group = handle['particle_plane']
        require(int(group['complete'][()]) == 1, 'Incomplete particle plane')
        records = group['records'][:]
        require(int(group['committed_records'][()]) == len(records),
                'Uncommitted particle records')
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            'Closed particle file changed')
    records = records[np.argsort(records['source_id'])]
    require(np.array_equal(records['source_id'], np.arange(1, len(records) + 1)),
            'Particle source IDs differ')
    require(np.all(records['crossing_direction'] == 1), 'Backward crossing')
    weight = abs(records['charge_C']).astype(np.longdouble)
    weight /= weight.sum()
    proper = records['proper_velocity'].astype(np.longdouble)
    angle = np.arctan2(proper[:, 0], proper[:, 2])
    gamma = np.sqrt(1 + np.sum(proper * proper, axis=1))
    mean = lambda value: float(np.sum(weight * value))
    return dict(count=len(records), mean_angle_rad=mean(angle),
                mean_gamma=mean(gamma), mean_x_m=mean(records['position_m'][:, 0]))


def gaussian_grid_reference(x, y, energies, p, sigma_xy, sigma_z, order,
                            time_first, time_span, ntime):
    """Vectorized Gaussian reference for all detector pixels."""
    nodes, weights = roots_hermitenorm(order)
    weights /= np.sqrt(2 * np.pi)
    ox, oy = np.meshgrid(nodes * sigma_xy, nodes * sigma_xy, indexing='xy')
    quadrature = np.outer(weights, weights).ravel()
    gx, gy = np.meshgrid(x, y, indexing='xy')
    points = np.stack((gx.ravel()[:, None] - ox.ravel()[None, :],
                       gy.ravel()[:, None] - oy.ravel()[None, :],
                       np.full((gx.size, ox.size), p['screen'])), axis=-1)
    flat = points.reshape(-1, 3)
    acceleration, acceleration_b = arc_spectrum(flat, energies, p)
    velocity, velocity_b = velocity_spectrum(
        flat, energies, p, time_first=time_first, time_span=time_span,
        ntime=ntime, sigma_time=sigma_z / (p['beta'] * c))
    straight, straight_b = straight_velocity_spectrum(
        flat, energies, p, time_first=time_first, time_span=time_span,
        ntime=ntime, sigma_time=sigma_z / (p['beta'] * c))
    target = (len(energies), gx.size, ox.size, 3)
    form = np.exp(-.5 * (np.asarray(energies) * e / hbar *
                         sigma_z / (p['beta'] * c))**2)
    average = lambda value: np.einsum(
        'fpqd,q->fpd', value.reshape(target), quadrature, optimize=True)
    acceleration = average(acceleration) * form[:, None, None]
    acceleration_b = average(acceleration_b) * form[:, None, None]
    velocity, velocity_b = average(velocity), average(velocity_b)
    straight, straight_b = average(straight), average(straight_b)
    shape = (len(energies), len(y), len(x), 3)
    reshape = lambda value: value.reshape(shape)
    return tuple(map(reshape, (acceleration + velocity - straight,
                               acceleration_b + velocity_b - straight_b,
                               acceleration, acceleration_b, velocity,
                               velocity_b, straight, straight_b)))


def centroid(image, x, y):
    positive = np.maximum(np.asarray(image, dtype=np.longdouble), 0)
    total = positive.sum(dtype=np.longdouble)
    require(total > 0, 'Nonpositive registered-band image')
    return (float(np.sum(positive * x[None, :], dtype=np.longdouble) / total),
            float(np.sum(positive * y[:, None], dtype=np.longdouble) / total))


def main(root):
    root = Path(root).resolve()
    campaign = json.loads((root / 'campaign.json').read_text())
    frozen = json.loads((root / 'FROZEN.json').read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            'Frozen input changed')
    shapes = campaign['shapes']
    require(all((root / f'ARRAY_COMPLETE_{index}.json').is_file()
                for index in range(len(shapes))), 'Array is incomplete')
    common, p = campaign['common'], campaign['common']['parameters']
    first = common['screen_first_relative_u_s']
    records, reference_grid = {}, None
    for shape in shapes:
        pair = {}
        for kind in ['straight', 'bend']:
            directory = root / shape / kind
            field = load_field(directory / 'detectors/screen.h5')
            grid = (field['times'], field['x'], field['y'], field['dx'], field['dy'])
            if reference_grid is None:
                reference_grid = grid
            else:
                require(all(np.array_equal(a, b) if isinstance(a, np.ndarray) else a == b
                            for a, b in zip(reference_grid, grid)),
                        'Bend/straight or shape grids differ')
            energy, fourier, density = transform(field, first)
            pair[kind] = dict(field=field, fourier=fourier, density=density)
        differential_cube = pair['bend']['field']['cube'] - pair['straight']['field']['cube']
        differential_fourier = pair['bend']['fourier'] - pair['straight']['fourier']
        differential_density = energy_density_spectrum(
            differential_fourier[..., :3], differential_fourier[..., 3:])
        differential_density[0] *= .5
        if len(pair['bend']['field']['times']) % 2 == 0:
            differential_density[-1] *= .5
        de = energy[1] - energy[0]
        check = parseval(differential_cube, differential_fourier,
                         differential_density, pair['bend']['field']['dt'], de,
                         pair['bend']['field']['dx'], pair['bend']['field']['dy'])
        require(check['relative_error'] < 2e-11, 'Differential Parseval failed')
        records[shape] = dict(pair=pair, differential_cube=differential_cube,
                              differential_fourier=differential_fourier,
                              differential_density=differential_density,
                              parseval=check)

    x, y, dx, dy = reference_grid[1], reference_grid[2], reference_grid[3], reference_grid[4]
    roi = ((abs(y[:, None]) < 10e-6) & (abs(x[None, :]) < 10e-6))
    require(abs(float(roi.sum()) * dx * dy - 4e-10) < 1e-20,
            'Registered 20 um square aperture changed')
    selected = (energy >= 4.) & (energy <= 18.)
    band = ((energy >= common['registered_band_eV'][0]) &
            (energy <= common['registered_band_eV'][1]))
    require(np.count_nonzero(band) >= 4, 'Insufficient registered bins')
    theory_energy = energy[selected]
    span = len(reference_grid[0]) * records[shapes[0]]['pair']['bend']['field']['dt']
    theory_values = gaussian_grid_reference(
        x, y, theory_energy, p, common['sigma_x_m'], common['sigma_z_m'], 8,
        first, span, 4 * len(reference_grid[0]))
    theory_e, theory_b = theory_values[:2]
    number = common['represented_electrons']
    theory_e *= number
    theory_b *= number
    theory_density = energy_density_spectrum(theory_e, theory_b)
    local_band = band[selected]
    aperture = lambda value: np.sum(value[:, roi], axis=1) * dx * dy
    theory_spectrum = aperture(theory_density)
    de = energy[1] - energy[0]
    theory_band_j = float(np.sum(theory_spectrum[local_band]) * de)
    require(theory_band_j > 0, 'Nonpositive paired theory band')

    center_j, center_i = np.argwhere(roi)[len(np.argwhere(roi)) // 2]
    refined = gaussian_paired_reference(
        x[center_i], y[center_j], theory_energy[local_band], p,
        sigma_xy=common['sigma_x_m'], sigma_z=common['sigma_z_m'], order=12,
        time_first=first, time_span=span, ntime=8 * len(reference_grid[0]))
    base_e = theory_e[local_band, center_j, center_i] / number
    base_b = theory_b[local_band, center_j, center_i] / number
    convergence = float(np.sqrt(
        np.linalg.norm(refined[0] - base_e)**2 +
        (c * np.linalg.norm(refined[1] - base_b))**2) /
        max(np.sqrt(np.linalg.norm(refined[0])**2 +
                    (c * np.linalg.norm(refined[1]))**2), 1e-300))

    rows, spectra, images = [], {}, {}
    theory_image = np.sum(theory_density[local_band], axis=0) * de
    theory_cx, theory_cy = centroid(theory_image, x, y)
    for shape in shapes:
        record = records[shape]
        spectrum = aperture(record['differential_density'][selected])
        spectra[shape] = spectrum
        image = np.sum(record['differential_density'][band], axis=0) * de
        images[shape] = image
        cx, cy = centroid(image, x, y)
        pic_band_j = float(np.sum(spectrum[local_band]) * de)
        raw_bend = float(np.sum(aperture(record['pair']['bend']['density'])[band]) * de)
        raw_straight = float(np.sum(aperture(record['pair']['straight']['density'])[band]) * de)
        l2 = float(np.linalg.norm(spectrum[local_band] - theory_spectrum[local_band]) /
                   np.linalg.norm(theory_spectrum[local_band]))
        straight_particles = load_particles(root / shape / 'straight/detectors/exit.h5')
        bend_particles = load_particles(root / shape / 'bend/detectors/exit.h5')
        rows.append(dict(shape=shape, paired_PIC_band_J=pic_band_j,
                         paired_theory_band_J=theory_band_j,
                         PIC_over_theory=pic_band_j / theory_band_j,
                         registered_spectral_relative_L2=l2,
                         raw_bend_band_J=raw_bend,
                         raw_straight_band_J=raw_straight,
                         differential_Parseval_relative_error=record['parseval']['relative_error'],
                         PIC_centroid_x_angle_rad=cx / p['screen'],
                         PIC_centroid_y_angle_rad=cy / p['screen'],
                         theory_centroid_x_angle_rad=theory_cx / p['screen'],
                         theory_centroid_y_angle_rad=theory_cy / p['screen'],
                         straight_exit_angle_rad=straight_particles['mean_angle_rad'],
                         bend_exit_angle_rad=bend_particles['mean_angle_rad'],
                         bend_exit_angle_error_rad=bend_particles['mean_angle_rad'] - p['theta']))

    with (root / 'summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez(root / 'paired_spectra.npz', photon_eV=theory_energy,
             theory_J_per_eV=theory_spectrum,
             **{shape.replace('-', '_') + '_J_per_eV': spectra[shape]
                for shape in shapes}, x_m=x, y_m=y,
             theory_registered_image_J_per_m2=theory_image,
             **{shape.replace('-', '_') + '_registered_image_J_per_m2': images[shape]
                for shape in shapes})

    omega = theory_energy * e / hbar
    sigma_t = common['sigma_z_m'] / (p['beta'] * c)
    form2 = np.exp(-(omega * sigma_t)**2)
    physical_over_continuum = (number + number * (number - 1) * form2) / (number**2 * form2)
    colors = {'linear': '#4c78a8', 'quintic-compensated': '#f58518',
              'linear-notch': '#54a24b', 'linear-binomial3': '#b279a2'}
    labels = {'linear': 'linear', 'quintic-compensated': 'quintic',
              'linear-notch': '7-point', 'linear-binomial3': '3-point'}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    axes[0, 0].plot(theory_energy, theory_spectrum * 1e18, 'k-', lw=2,
                    label='paired finite-arc LW')
    for shape in shapes:
        axes[0, 0].plot(theory_energy, spectra[shape] * 1e18, marker='o', ms=3,
                        color=colors[shape], label=labels[shape])
        ratio = np.divide(spectra[shape], theory_spectrum,
                          out=np.full_like(theory_spectrum, np.nan),
                          where=theory_spectrum != 0)
        axes[0, 1].plot(theory_energy, ratio, marker='o', ms=3,
                        color=colors[shape], label=labels[shape])
    for axis in axes[0]:
        axis.axvspan(*common['registered_band_eV'], color='grey', alpha=.12)
        axis.set_xlabel('Photon energy (eV)')
        axis.grid(alpha=.2)
    axes[0, 0].set(title='Bend minus straight: central 20 x 20 um aperture',
                   ylabel='aJ / eV')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(title='Paired PIC / paired finite-arc LW', ylabel='ratio')
    axes[0, 1].legend(fontsize=8)
    xloc = np.arange(len(shapes))
    axes[1, 0].bar(xloc, [row['PIC_over_theory'] for row in rows],
                   color=[colors[shape] for shape in shapes])
    axes[1, 0].axhline(1, color='black', lw=1)
    axes[1, 0].set(title='Registered 8-12 eV absolute energy ratio',
                   ylabel='paired PIC / theory')
    axes[1, 1].bar(xloc - .18,
                   [row['PIC_centroid_x_angle_rad'] * 1e6 for row in rows],
                   width=.36, label='PIC x')
    axes[1, 1].bar(xloc + .18,
                   [row['PIC_centroid_y_angle_rad'] * 1e6 for row in rows],
                   width=.36, label='PIC y')
    axes[1, 1].axhline(theory_cx / p['screen'] * 1e6, color='black',
                       ls='--', lw=1, label='theory x')
    axes[1, 1].axhline(0, color='grey', ls=':', lw=1, label='+z axis')
    axes[1, 1].set(title='Registered-band centroid about selected +z tangent',
                   ylabel='angle (urad)')
    axes[1, 1].legend(fontsize=8)
    for axis in axes[1]:
        axis.set_xticks(xloc, [labels[shape] for shape in shapes], rotation=15)
        axis.grid(alpha=.2)
    fig.suptitle('Gaussian finite-dipole synchrotron benchmark: four longitudinal shapes\n'
                 'matched bend-minus-straight complex-field subtraction; no fitted normalization')
    fig.savefig(root / 'synchrotron_four_shapes.png', dpi=200)
    fig.savefig(root / 'synchrotron_four_shapes.pdf')
    plt.close(fig)

    final = dict(complete=True, cases=rows, theory_self_check=self_check(),
                 reference_GH8_time4N_vs_GH12_time8N_relative_field_L2=convergence,
                 radiation_axis=common['radiation_axis'],
                 registered_band_eV=common['registered_band_eV'],
                 included_bin_centers_eV=energy[band].tolist(),
                 subtraction=('Complex bend Fourier fields minus matched straight Fourier fields '
                              'before Poynting spectrum; no constant, fitted background, or scalar normalization.'),
                 continuum_reference=True,
                 physical_shot_noise_correction=dict(
                     maximum_physical_over_continuum_in_registered_band=float(
                         np.max(physical_over_continuum[local_band])),
                     minimum_form_factor_squared_in_registered_band=float(np.min(form2[local_band])),
                     note=('Physical N+N(N-1)|F|^2 divided by simulated N^2|F|^2; '
                           'reported only and not added to PIC.')),
                 theory_reference=('Independent finite hard-edge circular arc plus exact finite-record '
                                   'velocity fields minus the identical incoming straight worldline; '
                                   'finite distance, aperture, and Gaussian form factor; no fitted amplitude.'),
                 physical_accuracy_accepted=False,
                 acceptance_note=('Single mesh and time step. Absolute agreement is measured; a later '
                                  'resolution study is required before physical accuracy acceptance.'))
    save(root / 'FINAL.json', final)
    require(all(sha(path) == digest for path, digest in frozen.items()),
            'Frozen input changed during analysis')
    print(json.dumps(final, indent=2), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
