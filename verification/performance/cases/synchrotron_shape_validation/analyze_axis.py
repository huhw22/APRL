"""Closed-only exact +z-axis spectrum using symmetric complex-field interpolation."""
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

from prepare_paired import require, save, sha
from theory import energy_density_spectrum
from theory_paired import gaussian_paired_reference


SHAPES = ['linear', 'quintic-compensated', 'linear-notch', 'linear-binomial3']


def load_axis(path):
    before = (path.stat().st_size, path.stat().st_mtime_ns)
    with h5py.File(path, 'r') as handle:
        group = handle['field_plane']
        require(int(group['complete'][()]) == 1, 'Incomplete field plane')
        n = int(group['committed_samples'][()])
        require(n == len(group['time_s']) == len(group['fields']), 'Uncommitted field tail')
        time = group['time_s'][:]
        nx, ny = int(group.attrs['nx']), int(group.attrs['ny'])
        dx, dy = float(group.attrs['dx_m']), float(group.attrs['dy_m'])
        x = float(group.attrs['x_first_m']) + np.arange(nx) * dx
        y = float(group.attrs['y_first_m']) + np.arange(ny) * dy
        ix = np.argsort(abs(x))[:2]
        iy = np.argsort(abs(y))[:2]
        require(np.allclose(np.sort(x[ix]), [-dx / 2, dx / 2], rtol=0, atol=1e-18) and
                np.allclose(np.sort(y[iy]), [-dy / 2, dy / 2], rtol=0, atol=1e-18),
                'The four nearest samples do not symmetrically bracket +z')
        ix, iy = np.sort(ix), np.sort(iy)
        records = group['fields'][:]
        electric = records['electric_V_per_m'][:, iy][:, :, ix].mean(axis=(1, 2))
        magnetic = records['magnetic_T'][:, iy][:, :, ix].mean(axis=(1, 2))
        require(np.isfinite(electric).all() and np.isfinite(magnetic).all(),
                'Nonfinite axis field')
    require(before == (path.stat().st_size, path.stat().st_mtime_ns),
            'Closed field file changed')
    dt = float((time[-1] - time[0]) / (len(time) - 1))
    jitter = float(np.max(abs(time - time[0] - np.arange(len(time)) * dt)) / dt)
    require(jitter < 1e-4, 'Nonuniform time samples')
    return dict(time=time, dt=dt, electric=electric, magnetic=magnetic,
                x_samples_m=x[ix].tolist(), y_samples_m=y[iy].tolist(),
                interpolation='equal-weight bilinear value at x=y=0')


def fourier(data, first):
    n, dt = len(data['time']), data['dt']
    energy = h * np.fft.rfftfreq(n, dt) / e
    phase = np.exp(1j * energy * e / hbar * first)[:, None]
    ef = np.conj(np.fft.rfft(data['electric'], axis=0)) * dt * phase
    bf = np.conj(np.fft.rfft(data['magnetic'], axis=0)) * dt * phase
    return energy, ef, bf


def density(ef, bf, n):
    value = energy_density_spectrum(ef, bf)
    value[0] *= .5
    if n % 2 == 0:
        value[-1] *= .5
    return value


def field_difference(first, second):
    return float(np.sqrt(np.linalg.norm(first[0] - second[0])**2 +
                         (c * np.linalg.norm(first[1] - second[1]))**2) /
                 max(np.sqrt(np.linalg.norm(second[0])**2 +
                             (c * np.linalg.norm(second[1]))**2), 1e-300))


def main(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    require(not output.exists(), 'Axis output already exists')
    output.mkdir(parents=True)
    campaign = json.loads((root / 'campaign.json').read_text())
    frozen = json.loads((root / 'FROZEN.json').read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()),
            'Frozen campaign input changed')
    require(all((root / f'ARRAY_COMPLETE_{index}.json').is_file()
                for index in range(4)), 'PIC campaign is incomplete')
    common, p = campaign['common'], campaign['common']['parameters']
    first = common['screen_first_relative_u_s']
    values, time_reference, energy = {}, None, None
    for shape in SHAPES:
        pair = {}
        for kind in ['straight', 'bend']:
            data = load_axis(root / shape / kind / 'detectors/screen.h5')
            if time_reference is None:
                time_reference = data['time']
            else:
                require(np.array_equal(data['time'], time_reference),
                        'Paired timestamps differ')
            photon, ef, bf = fourier(data, first)
            if energy is None:
                energy = photon
            else:
                require(np.array_equal(photon, energy), 'Frequency grids differ')
            pair[kind] = dict(data=data, ef=ef, bf=bf,
                              spectrum=density(ef.copy(), bf, len(data['time'])))
        ef, bf = pair['bend']['ef'] - pair['straight']['ef'], pair['bend']['bf'] - pair['straight']['bf']
        spectrum = density(ef.copy(), bf, len(time_reference))
        e_time = pair['bend']['data']['electric'] - pair['straight']['data']['electric']
        b_time = pair['bend']['data']['magnetic'] - pair['straight']['data']['magnetic']
        time_integral = float(np.sum(e_time[:, 0] * b_time[:, 1] -
                                     e_time[:, 1] * b_time[:, 0]) *
                              pair['bend']['data']['dt'] / mu_0)
        frequency_integral = float(np.sum(spectrum) * (energy[1] - energy[0]))
        parseval = abs(time_integral - frequency_integral) / max(
            abs(time_integral), abs(frequency_integral), 1e-300)
        require(parseval < 2e-11, 'Axis Parseval failed')
        values[shape] = dict(pair=pair, ef=ef, bf=bf, spectrum=spectrum,
                             parseval=parseval)

    # The card's explicit radiation-resolution contract ends at 12 eV.
    # Frequencies above it are neither accepted nor displayed.
    selected = ((energy >= 4.) &
                (energy <= common['registered_band_eV'][1]))
    band = ((energy >= common['registered_band_eV'][0]) &
            (energy <= common['registered_band_eV'][1]))
    local_band = band[selected]
    photon = energy[selected]
    span = len(time_reference) * values[SHAPES[0]]['pair']['bend']['data']['dt']
    references = {}
    for order, multiplier in [(20, 8), (24, 4), (24, 8)]:
        result = gaussian_paired_reference(
            0., 0., photon, p, sigma_xy=common['sigma_x_m'],
            sigma_z=common['sigma_z_m'], order=order, time_first=first,
            time_span=span, ntime=multiplier * len(time_reference))
        references[(order, multiplier)] = (result[0], result[1])
    final_e, final_b = references[(24, 8)]
    transverse_convergence = field_difference(references[(20, 8)],
                                              references[(24, 8)])
    time_convergence = field_difference(references[(24, 4)],
                                        references[(24, 8)])
    number = common['represented_electrons']
    theory = energy_density_spectrum(final_e * number, final_b * number)
    theory_band = float(np.sum(theory[local_band]) * (energy[1] - energy[0]))
    require(theory_band > 0, 'Nonpositive axis theory')

    rows, spectra = [], {}
    for shape in SHAPES:
        spectrum = values[shape]['spectrum'][selected]
        spectra[shape] = spectrum
        pic_band = float(np.sum(spectrum[local_band]) * (energy[1] - energy[0]))
        straight_transverse = np.sqrt(
            np.sum(abs(values[shape]['pair']['straight']['ef'][selected][local_band, :2])**2) +
            c**2 * np.sum(abs(values[shape]['pair']['straight']['bf'][selected][local_band, :2])**2))
        paired_transverse = np.sqrt(
            np.sum(abs(values[shape]['ef'][selected][local_band, :2])**2) +
            c**2 * np.sum(abs(values[shape]['bf'][selected][local_band, :2])**2))
        rows.append(dict(shape=shape, axis_PIC_band_J_per_m2=pic_band,
                         axis_theory_band_J_per_m2=theory_band,
                         PIC_over_theory=pic_band / theory_band,
                         registered_spectral_relative_L2=float(
                             np.linalg.norm(spectrum[local_band] - theory[local_band]) /
                             np.linalg.norm(theory[local_band])),
                         straight_control_transverse_field_over_paired=float(
                             straight_transverse / paired_transverse),
                         axis_Parseval_relative_error=values[shape]['parseval']))

    with (output / 'axis_summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output / 'axis_spectra.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['photon_eV', 'theory_J_per_eV_m2'] + SHAPES)
        writer.writerows(zip(photon, theory, *(spectra[shape] for shape in SHAPES)))
    np.savez(output / 'axis_spectra.npz', photon_eV=photon,
             theory_J_per_eV_m2=theory,
             **{shape.replace('-', '_') + '_J_per_eV_m2': spectra[shape]
                for shape in SHAPES})

    colors = {'linear': '#4c78a8', 'quintic-compensated': '#f58518',
              'linear-notch': '#54a24b', 'linear-binomial3': '#b279a2'}
    labels = {'linear': 'linear', 'quintic-compensated': 'quintic',
              'linear-notch': '7-point', 'linear-binomial3': '3-point'}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    axes[0, 0].plot(photon, theory * 1e12, 'k-', lw=2, label='finite-arc LW')
    for shape in SHAPES:
        axes[0, 0].plot(photon, spectra[shape] * 1e12, marker='o', ms=3,
                        color=colors[shape], label=labels[shape])
        ratio = np.divide(spectra[shape], theory,
                          out=np.full_like(theory, np.nan), where=theory != 0)
        axes[0, 1].plot(photon, ratio, marker='o', ms=3,
                        color=colors[shape], label=labels[shape])
    for axis in axes[0]:
        axis.axvspan(*common['registered_band_eV'], color='grey', alpha=.12)
        axis.set_xlabel('Photon energy (eV)')
        axis.grid(alpha=.2)
    axes[0, 0].set(title='Exact +z-axis spectral fluence',
                   ylabel='pJ / (eV m2)')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(title='Paired PIC / finite-arc LW', ylabel='ratio')
    axes[0, 1].legend(fontsize=8)
    xloc = np.arange(len(SHAPES))
    axes[1, 0].bar(xloc, [row['PIC_over_theory'] for row in rows],
                   color=[colors[shape] for shape in SHAPES])
    axes[1, 0].axhline(1, color='black', lw=1)
    axes[1, 0].set(title='Integrated 8-12 eV on-axis fluence ratio',
                   ylabel='PIC / theory')
    axes[1, 1].bar(xloc,
                   [row['straight_control_transverse_field_over_paired'] for row in rows],
                   color=[colors[shape] for shape in SHAPES])
    axes[1, 1].set(title='Residual straight-control transverse field',
                   ylabel='straight / paired field norm')
    for axis in axes[1]:
        axis.set_xticks(xloc, [labels[shape] for shape in SHAPES], rotation=15)
        axis.grid(alpha=.2)
    fig.suptitle('Gaussian finite-dipole synchrotron spectrum on the automatically selected +z tangent\n'
                 'symmetric complex-field interpolation; no fitted normalization')
    fig.savefig(output / 'axis_four_shapes.png', dpi=200)
    fig.savefig(output / 'axis_four_shapes.pdf')
    plt.close(fig)

    omega = photon * e / hbar
    sigma_t = common['sigma_z_m'] / (p['beta'] * c)
    form2 = np.exp(-(omega * sigma_t)**2)
    physical = (number + number * (number - 1) * form2) / (number**2 * form2)
    final = dict(complete=True, status='EXPLORATORY_SINGLE_GRID_AXIS_COMPARISON',
                 cases=rows, radiation_axis=common['radiation_axis'],
                 axis_interpolation=dict(x_samples_m=values[SHAPES[0]]['pair']['bend']['data']['x_samples_m'],
                                         y_samples_m=values[SHAPES[0]]['pair']['bend']['data']['y_samples_m'],
                                         target_m=[0., 0.], method='bilinear complex E/B before Poynting spectrum'),
                 theory_order20_vs_order24_relative_field_L2=transverse_convergence,
                 theory_time4N_vs_time8N_relative_field_L2=time_convergence,
                 displayed_and_analyzed_photon_range_eV=[4., common['registered_band_eV'][1]],
                 included_bin_centers_eV=energy[band].tolist(),
                 subtraction='bend minus straight complex Fourier E/B; no fitted normalization',
                 physical_shot_noise_maximum_correction=float(np.max(physical[local_band])),
                 physical_accuracy_accepted=False,
                 acceptance_note=('The on-axis theory quadrature is checked, but this is one PIC mesh/time step. '
                                  'A resolution scan is still required for physical accuracy acceptance.'),
                 supersedes_for_axis_comparison=str(root / 'FINAL.json'),
                 superseded_reason=('The earlier 20x20 um area integral includes the straight control beam crossing '
                                    'the aperture; its transverse Gaussian quadrature did not converge.'))
    save(output / 'FINAL.json', final)
    require(all(sha(path) == digest for path, digest in frozen.items()),
            'Frozen campaign changed during axis analysis')
    print(json.dumps(final, indent=2), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
