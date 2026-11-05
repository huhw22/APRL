"""Freeze a central three-slippage SPARC first-segment control pair."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil

import h5py
import numpy as np
from scipy.constants import c, e
import yaml


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + '\n')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: native(v) for k, v in value.items()}
    if isinstance(value, list):
        return [native(v) for v in value]
    return value


def prepare(workspace, output, binary, build_record):
    workspace, output, binary = (Path(x).resolve() for x in (workspace, output, binary))
    output.mkdir(parents=True, exist_ok=False)
    old = workspace / 'results/sparc_sase_comparison_20261103'
    production = old / 'compact-production-r1'
    prior_plan = json.loads((production / 'plan.json').read_text())
    original = old / 'compact-preflight-r1/particles.h5'
    original_sha = sha(original)
    slices, nper, nbins = 231, 1024, 16
    total = prior_plan['active_slices']
    first = (total - slices) // 2
    last = first + slices
    inputs = output / 'inputs'
    inputs.mkdir()
    particles = inputs / 'particles.h5'
    with h5py.File(original, 'r') as f:
        source = f['particles']
        a = source['records'][first*nper:last*nper]
        attrs = dict(source.attrs)
    require(len(a) == slices*nper, 'Input crop was truncated')
    require(np.array_equal(a['source_id'], np.arange(first*nper+1, last*nper+1)), 'Source IDs changed')
    for field in ['proper_velocity', 'plane_position_m', 'macro_weight']:
        grouped = a[field].reshape((-1, nbins) + a[field].shape[1:])
        require(np.all(grouped == grouped[:, :1]), 'Phase family was split: ' + field)
    with h5py.File(particles, 'x') as f:
        group = f.create_group('particles')
        for key, value in attrs.items():
            group.attrs[key] = value
        group.attrs['selection'] = 'Original central slices 577:808; complete 16-clone families; no reweighting or phase shift'
        group.create_dataset('records', data=a, chunks=(min(65536, len(a)),))
    u = a['proper_velocity']
    gamma = np.sqrt(1 + np.sum(u*u, axis=1))
    beta = u[:, 2] / gamma
    weights = a['macro_weight'].astype(np.longdouble)
    times = a['arrival_time_offset_s'].astype(np.longdouble)
    ne = weights.sum()
    gb = prior_plan['physics']['boost_gamma']
    bb = np.sqrt(1-1/gb**2)
    dz, dt = prior_plan['geometry']['dz_m'], prior_plan['geometry']['dt_s']
    wavelength = prior_plan['physics']['wavelength_m']
    # The held trajectory is built from the lab longitudinal momentum only,
    # unlike the active particle's full box-frame velocity.
    light = (np.hypot(1., u[:, 2])+u[:, 2])/(gb+np.sqrt(gb-1)*np.sqrt(gb+1))
    held_uz = .5*(light-1/light)
    held = held_uz/np.hypot(1., held_uz)
    held_min, held_max = float(held.min()), float(held.max())
    mean_beta = np.sum(weights*beta)/ne
    mean_bt = np.sum(weights*beta*times)/ne
    tmin, tmax = float(times.min()), float(times.max())
    center = float(c*(mean_beta*np.longdouble(tmax)-mean_bt)+2*dz/gb)
    reference_time = float((np.longdouble(center)/c+mean_bt)/mean_beta)
    snapshot_z = beta*c*(reference_time-a['arrival_time_offset_s'])
    input_times_lab = -snapshot_z/(c*beta)
    origin_time = float(np.min(input_times_lab-2*dz/(gb*c*(bb+held))))
    initial_z = -gb*c*(bb+held)*(input_times_lab-origin_time)
    zmax, zmin = float(initial_z.max()), float(initial_z.min())
    translation = -.5*(zmin+zmax)
    last_z = .01 + 77*.028
    span = tmax-tmin+77*wavelength/c+(1-bb)*last_z/c+4e-15
    steps = int(np.ceil(span*gb/dt))
    active = int(np.ceil(last_z/(gb*dz))) + 64+128+2
    extent = max((1+held_max)*(steps+2)+66-translation/dz,
                 translation/dz-64+active+2)
    full = int(np.ceil(2*extent/256))*256
    # The crop has its own frame/grid phase. Off/on have identical geometry;
    # no claim is made that a cropped bunch reproduces the original phase.
    # Changing mesh.center does not align it: boxTranslation moves with center.
    timing_reference = json.loads((old / 'energy_critical_mode_20261105/onset/RESULT.json').read_text())
    old_fraction = timing_reference['plane_event_support']['entry']['window_cell_at_step0'] % 1
    default_fraction = (translation/dz+full/2) % 1
    mesh_center_z = 0.
    duration = (steps+.5)*dt
    stop = gb*gb*(center+bb*c*origin_time)+np.sqrt(gb*gb-1)*dz*steps
    geometry = dict(nx=80, ny=80, reference_nz=full, active_nz=active,
        dx_m=prior_plan['geometry']['dx_m'], dy_m=prior_plan['geometry']['dy_m'], dz_m=dz, dt_box_s=dt,
        boost_gamma=gb, boost_beta=bb, steps=steps, duration_box_s=duration,
        first_magnet_z_m=.01, last_magnet_z_m=last_z,
        center_lab_z_m=center, reference_time_offset_lab_s=reference_time,
        predicted_frame_time_origin_lab_s=origin_time,
        predicted_box_translation_m=translation, mesh_center_z_m=mesh_center_z,
        original_entry_fraction=old_fraction, predicted_entry_fraction=default_fraction,
        phase_alignment_claim=False, paired_grid_geometry_identical=True,
        held_beta_min=held_min, held_beta_max=held_max,
        field_cells_times_steps_ratio_to_old_full=active*steps/(prior_plan['geometry']['active_nz']*prior_plan['geometry']['steps']),
        field_raw_upper_GiB_per_case=3*(steps+1)*80*80*48/2**30)
    selection = dict(original=str(original), original_sha256=original_sha,
        first_original_slice=first, end_original_slice_exclusive=last,
        slices=slices, particles=len(a), represented_electrons=float(ne), charge_C=float(ne)*e,
        input_time_span_s=tmax-tmin, lab_length_m=(tmax-tmin)*c,
        nominal_slippage_m=77*wavelength, slippage_multiple=3,
        macro_weight_min=float(weights.min()), macro_weight_max=float(weights.max()),
        times_rephased=False, weights_renormalized=False, source_ids_preserved=True,
        core_first_slice_inclusive=first+77, core_end_slice_exclusive=last-77,
        limitations=['Hard crop creates a new bunch edge; paired controls are comparable, but old full-beam output is not the same-input baseline.',
                    'Primary BF comparison must report central 77 slices as well as the full crop.',
                    'If the unfiltered crop does not reproduce fast startup, no conclusion about the original full-beam mechanism is accepted.'])
    save(output / 'GEOMETRY.json', geometry)
    save(output / 'INPUT_SELECTION.json', selection)
    template = Path(__file__).resolve().parent
    scripts = output / 'scripts'
    scripts.mkdir()
    for name in ['prepare.py', 'run_case.py', 'run_case.sbatch', 'analyze.py']:
        shutil.copy2(template/name, scripts/name)
    prior = template.parent
    for folder, name in [('space_charge_physics', 'launcher.py'),
                         ('source_coupling_integration', 'binding.py'),
                         ('source_coupling_integration', 'analyze_closed.py')]:
        shutil.copy2(prior/folder/name, scripts/name)
    for name in ['report.py', 'review_discrete.py']:
        shutil.copy2(prior/'discrete_energy_ledger'/name, scripts/('report_discrete.py' if name == 'report.py' else name))
    base = yaml.safe_load((production/'energy/run.yaml').read_text())
    cases = []
    for stage in ['full']:
        for mode in ['off', 'on']:
            name = ('timing-' if stage == 'timing' else '')+mode
            directory = output/name
            directory.mkdir()
            cfg = copy.deepcopy(base)
            cfg['output']['manifest'] = str(directory/'run-manifest.yaml')
            cfg['mesh'].update(cells=[80, 80, full], center=[0., 0., mesh_center_z],
                duration=duration, longitudinal_shape='linear' if mode == 'off' else 'linear-notch')
            cfg['mesh']['moving_window'].update(cells=active, particle_collar_cells=8,
                reference_directory=str(directory/'reference'))
            cfg['beam']['reference']['initial_center_z'] = center
            cfg['beam']['input'].update(file=str(particles), electrons=float(ne))
            cfg['sources']['magnetic_elements'] = [cfg['sources']['magnetic_elements'][0]]
            cfg['detectors'].update(directory=str(directory/'detectors'), substep_particle_planes=True,
                field_planes=[dict(name=n, z=z, rhythm=dt/gb, buffer_samples=8, compression=1,
                                   particle_background={'enabled': False})
                              for n,z in [('entry',0.), ('und-entry',.01), ('p1',last_z)]],
                particle_planes=[dict(name=n+'-particles', z=z, buffer_records=32768, compression=0)
                                 for n,z in [('entry',0.), ('und-entry',.01), ('mid',.01+77*.028/2), ('p1',last_z)]],
                energy_planes=[dict(name=n+'-transport', z=z, rhythm=dt/gb, buffer_samples=64, compression=0)
                               for n,z in [('entry',0.), ('und-entry',.01), ('p1',last_z)]])
            cfg['energy_ledger'].update(directory=str(directory/'ledger'), enabled=True,
                sample_interval_steps=10, discrete_operator_audit=True)
            cfg['runtime']['resource_monitor']['progress_interval_steps'] = 50
            cfg['trajectory']['enabled'] = False
            expected_steps = steps
            case_stop = stop
            cfg['stop'] = {'mode': 'reference-center-z', 'z': case_stop}
            (directory/'run.yaml').write_text(yaml.safe_dump(native(cfg), sort_keys=False))
            for folder in ['detectors', 'ledger']:
                (directory/folder).mkdir()
            plan = dict(name=name, mode=mode, stage=stage, directory=str(directory),
                binary=str(binary), build_record=str(Path(build_record).resolve()), nodes=1, ranks=16, threads=8,
                input=str(particles), input_sha256=sha(particles), expected_particles=len(a),
                expected_steps=expected_steps, configured_duration_step_upper=math.ceil(cfg['mesh']['duration']/dt),
                case_limit_seconds=3300,
                require_all_particles_through_planes=True, wavelength_m=wavelength,
                first_original_slice=first, end_original_slice_exclusive=last,
                core_first_slice_inclusive=first+77, core_end_slice_exclusive=last-77)
            save(directory/'run-plan.json', plan)
            cases.append(plan)
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(Path(build_record).resolve()), seed='Original physical-noise seed 20261103; unchanged selected events',
        cases=cases, geometry=geometry, input_selection=selection,
        smoke_record=str(output.parent/'smoke-r1/COMPLETE.json'),
        accepted_for_execution=False,
        timing_method='Observe the early progress of the valid full off case; no upstream-of-magnet stop card is used.',
        reason='Program/module build and smoke identity must pass before explicit full-control submission')
    save(output/'plan.json', plan)
    paths = [particles, output/'plan.json', output/'INPUT_SELECTION.json', output/'GEOMETRY.json', *scripts.iterdir()]
    for case in cases:
        paths.extend([Path(case['directory'])/'run.yaml', Path(case['directory'])/'run-plan.json'])
    save(output/'FROZEN_INPUTS.json', {str(p):sha(p) for p in paths if p.is_file()})
    require(sha(original) == original_sha, 'Original input changed')
    print(json.dumps(geometry, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('workspace', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--build-record', type=Path, required=True)
    a = p.parse_args()
    prepare(a.workspace, a.output, a.binary, a.build_record)
