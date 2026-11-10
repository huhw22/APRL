"""Prepare paired bend/straight frozen-injection cases for four shapes."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil

import h5py
import numpy as np
from scipy.constants import c, e, hbar, m_e
from scipy.special import roots_hermitenorm
import yaml


SHAPES = ['linear', 'quintic-compensated', 'linear-notch', 'linear-binomial3']
GBOOST = 929.1305492858543
DZ = 1.2543262415359034e-6
DT = DZ / c


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value


def save(path, value):
    Path(path).write_text(json.dumps(native(value), indent=2, allow_nan=False) + '\n')


def rewrite(value, old, new):
    if isinstance(value, dict):
        return {key: rewrite(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, old, new) for item in value]
    if isinstance(value, str):
        return value.replace(str(old), str(new))
    return value


def parameters(gamma=1174.171072271572, rho=3., length=.0024, screen=.3024):
    beta = math.sqrt((gamma - 1) * (gamma + 1)) / gamma
    theta = math.asin(length / rho)
    return dict(gamma=gamma, beta=beta, rho=rho, length=length, screen=screen,
                theta=theta, field_T=gamma * beta * m_e * c / (e * rho),
                critical_eV=1.5 * hbar * c * gamma**3 / (rho * e))


def prepare(workspace, output, binary, build_record):
    workspace, output, binary, build_record = [Path(x).resolve() for x in
                                                (workspace, output, binary, build_record)]
    require(not output.exists(), 'Output directory already exists')
    identity = json.loads(build_record.read_text())
    require(identity.get('passed') is True and Path(identity['binary']).resolve() == binary,
            'Unaccepted binary identity')
    require(sha(binary) == identity['binary_sha256'], 'Binary hash differs')
    output.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    prior = here.parent / 'radiation_method_validation'
    scripts = output / 'scripts'
    scripts.mkdir()
    for name in ['README.md', 'prepare_paired.py', 'run_paired.py',
                 'aggregate_paired.py', 'theory_paired.py']:
        shutil.copy2(here / name, scripts / name)
    # The preserved launcher imports the campaign helpers under this historic
    # module name.  Freeze an identical alias rather than editing the launcher.
    shutil.copy2(here / 'prepare_paired.py', scripts / 'prepare.py')
    for name in ['theory.py', 'verify_reference.py', 'launcher.py']:
        shutil.copy2(prior / name, scripts / name)

    p = parameters()
    sx = sy = 20e-6
    sz = 20e-9
    represented = 1e4
    nodes, weights = [], []
    for order in [8, 8, 16]:
        node, weight = roots_hermitenorm(order)
        nodes.append(node)
        weights.append(weight / math.sqrt(2 * math.pi))
    xx, yy, zz = np.meshgrid(*nodes, indexing='ij')
    wx, wy, wz = np.meshgrid(*weights, indexing='ij')
    weight = (wx * wy * wz).ravel()
    weight /= weight.sum()
    count = len(weight)
    dtype = np.dtype([('plane_position_m', '<f8', (2,)),
                      ('arrival_time_offset_s', '<f8'),
                      ('proper_velocity', '<f8', (3,)),
                      ('source_id', '<u8'), ('macro_weight', '<f8')])
    particles = np.zeros(count, dtype=dtype)
    particles['plane_position_m'] = np.c_[xx.ravel() * sx, yy.ravel() * sy]
    particles['arrival_time_offset_s'] = zz.ravel() * sz / (p['beta'] * c)
    particles['proper_velocity'][:, 2] = math.sqrt((p['gamma'] - 1) * (p['gamma'] + 1))
    particles['source_id'] = np.arange(1, count + 1)
    particles['macro_weight'] = weight * count
    inputs = output / 'inputs'
    inputs.mkdir()
    input_file = inputs / 'gaussian_plane_events.h5'
    with h5py.File(input_file, 'x') as handle:
        group = handle.create_group('particles')
        group.attrs['format_version'] = 4
        group.attrs['scope'] = ('Product Gauss-Hermite weighted Gaussian continuum; cold parallel '
                                'fixed-lab-plane crossing events shared by bend and straight controls')
        group.create_dataset('records', data=particles)

    gb = math.sqrt((GBOOST - 1) * (GBOOST + 1))
    bb = gb / GBOOST
    uz_lab = np.longdouble(particles['proper_velocity'][0, 2])
    plus = ((np.sqrt(1 + uz_lab * uz_lab) + uz_lab) /
            (np.longdouble(GBOOST) + np.longdouble(gb)))
    uz_box = (plus - 1 / plus) / 2
    held = float(uz_box / np.sqrt(1 + uz_box * uz_box))
    arrival = particles['arrival_time_offset_s'].astype(np.longdouble)
    mean_arrival = np.sum(weight.astype(np.longdouble) * arrival)
    center = float(max(np.longdouble(1e-7),
                       np.longdouble(p['beta']) * c * (np.max(arrival) - mean_arrival) +
                       256 * np.finfo(float).eps))
    reference_time = float(np.longdouble(center) / (np.longdouble(p['beta']) * c) + mean_arrival)
    frame_time = float(np.min(arrival - reference_time -
                              2 * DZ / (GBOOST * c * (bb + held))))
    source_z = (-np.longdouble(GBOOST) * c * (bb + held) *
                (arrival - reference_time - frame_time))
    translation = float(-.5 * (source_z.min() + source_z.max()))
    entrance = .05
    screen_z = entrance + p['screen']
    steps = math.ceil(5e-15 * GBOOST / DT)
    duration = (steps + .5) * DT
    active = math.ceil(screen_z / (GBOOST * DZ)) + 64 + 128 + 2
    required = max((1 + held) * (steps + 2) + 66 - translation / DZ,
                   translation / DZ - 64 + active + 2)
    full = math.ceil(2 * required / 256) * 256
    origin = -.5 * full * DZ
    offset = math.floor((translation - origin) / DZ) - 64
    require(offset - (1 + held) * (steps + 2) - 1 > 0 and offset + active < full,
            'Moving-window reference extent is insufficient')
    initial_reference = GBOOST * GBOOST * (center + bb * c * frame_time)
    stop_z = initial_reference + gb * DZ * steps
    form_errors = []
    for photon in [8., 10., 12.]:
        for nx in [-2e-4, 0., 2e-4]:
            omega = photon * e / hbar
            k = omega / c
            phase = omega * particles['arrival_time_offset_s']
            phase -= k * nx * particles['plane_position_m'][:, 0]
            finite = np.sum(weight * np.exp(1j * phase))
            exact = np.exp(-.5 * ((omega * sz / (p['beta'] * c))**2 +
                                  (k * nx * sx)**2))
            form_errors.append(float(abs(finite - exact)))
    require(max(form_errors) < 1e-8, 'Gaussian form-factor quadrature failed')

    old_template = workspace / 'results/radiation_method_validation_20261103/synchrotron-r1/template'
    base = yaml.safe_load((old_template / 'run.yaml').read_text())
    common = dict(binary=str(binary), binary_sha256=sha(binary),
        build_record=str(build_record), binding_probe=str(workspace / 'build/cluster-tools/binding_probe'),
        macroparticles=count, represented_electrons=represented, input=str(input_file),
        input_sha256=sha(input_file), parameters=p, sigma_x_m=sx, sigma_y_m=sy,
        sigma_z_m=sz, magnet_entrance_z_m=entrance, screen_lab_z_m=screen_z,
        reference_time_offset_s=reference_time, reference_frame_time_origin_s=frame_time,
        screen_first_relative_u_s=(frame_time + bb * screen_z / c -
                                  (-reference_time + entrance / (p['beta'] * c)) - p['screen'] / c),
        registered_band_eV=[8., 12.], registered_aperture_m=[-10e-6, 10e-6, -10e-6, 10e-6],
        radiation_axis=dict(selected_emission_z_m=entrance, tangent=[0., 0., 1.],
                            screen_intersection_m=[0., 0., screen_z], aperture_center_m=[0., 0.],
                            automatically_computed=True),
        paired_subtraction=('Complex Fourier field of bend minus same-shape straight control before '
                            'forming Poynting spectrum; theory uses the identical full-LW difference.'),
        continuum_reference=('N^2 Gaussian coherent continuum; physical N(1-|F|^2) shot term is not simulated'),
        geometry=dict(reference_nz=full, active_nz=active, nx=160, ny=128,
                      field_steps=steps, mesh_duration_s=duration,
                      stop_reference_lab_z_m=stop_z),
        form_factor_quadrature_max_absolute_error=max(form_errors), timeout_seconds=1800,
        physical_acceptance=False)
    cases = []
    for shape in SHAPES:
        for kind in ['straight', 'bend']:
            directory = output / shape / kind
            directory.mkdir(parents=True)
            cfg = rewrite(copy.deepcopy(base), old_template, directory)
            cfg['output']['manifest'] = str(directory / 'run-manifest.yaml')
            cfg['mesh'].update(cells=[160, 128, full], duration=duration,
                               particle_coupling='conservative-source', longitudinal_shape=shape)
            cfg['mesh']['moving_window'].update(cells=active,
                                                reference_directory=str(directory / 'reference'))
            cfg['runtime']['resource_monitor']['progress_interval_steps'] = 50
            cfg['beam']['reference'].update(input_plane_z=0., initial_center_z=center,
                                             frozen_lab_plane_injection=True,
                                             hold_transverse_until_activation=True,
                                             activation_plane_z=0.)
            cfg['beam']['input'].update(file=str(input_file), electrons=represented)
            zero = dict(type='planar-undulator', entrance_z=.02, period=.009, periods=1,
                        strength_parameter=0., gaussian_fringe=False)
            dipole = dict(type='uniform-dipole', entrance_z=entrance, length=p['length'],
                           field_T=p['field_T'], polarization_angle_rad=math.pi / 2)
            cfg['sources']['magnetic_elements'] = ([dipole, zero] if kind == 'bend' else [zero])
            cfg['detectors'].update(directory=str(directory / 'detectors'),
                field_planes=[dict(name='screen', z=screen_z, rhythm=4 * DT / GBOOST,
                                   buffer_samples=8, compression=0,
                                   particle_background={'enabled': False})],
                particle_planes=[dict(name=name, z=z, buffer_records=8192, compression=0)
                                 for name, z in [('entry', 0.),
                                                 ('bend-exit', entrance + p['length']),
                                                 ('exit', screen_z)]],
                energy_planes=[dict(name=name, z=z, rhythm=4 * DT / GBOOST,
                                    buffer_samples=64, compression=0)
                               for name, z in [('entry-energy', 0.), ('exit-energy', screen_z)]])
            cfg['radiation_resolution']['maximum_photon_energy_eV'] = 12.
            cfg['energy_ledger'].update(directory=str(directory / 'ledger'), enabled=True,
                                        sample_interval_steps=10, discrete_operator_audit=True)
            cfg['trajectory']['enabled'] = False
            cfg['stop'] = {'mode': 'reference-center-z', 'z': stop_z}
            (directory / 'run.yaml').write_text(yaml.safe_dump(native(cfg), sort_keys=False))
            case_plan = dict(common, shape=shape, kind=kind, case_directory=str(directory),
                             expected_steps=steps, expected_particles=count)
            save(directory / 'PLAN.json', case_plan)
            cases.append(dict(shape=shape, kind=kind, directory=str(directory), timeout_seconds=1800))
    campaign = dict(common=common, shapes=SHAPES, cases=cases, binary=str(binary),
                    binding_probe=common['binding_probe'],
                    scope=('Four same-coupling longitudinal shapes; paired bend/straight complex-field '
                           'subtraction; selected radiation tangent is +z'))
    save(output / 'campaign.json', campaign)
    frozen_paths = [output / 'campaign.json', input_file, binary, build_record]
    frozen_paths += [path for path in scripts.iterdir() if path.is_file()]
    for case in cases:
        frozen_paths += [Path(case['directory']) / 'run.yaml', Path(case['directory']) / 'PLAN.json']
    save(output / 'FROZEN.json', {str(path): sha(path) for path in frozen_paths})
    save(output / 'PREPARED.json', dict(complete=True, paired_cases=len(cases), common=common,
                                        production_default_changed=False, source_core_changed=False))
    print(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-record', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.workspace, args.output, args.binary, args.build_record)
