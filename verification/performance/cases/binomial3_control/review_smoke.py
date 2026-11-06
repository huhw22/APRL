"""Independent closed-smoke gate; refuses to inspect H5 before normal completion."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np
import yaml


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def source_convergence(directory, ranks, steps):
    expected = {f'source-coupling-rank-{rank:05d}' for rank in range(ranks)}
    files = sorted((directory / 'reference').glob('source-coupling-rank-*.csv'))
    require({path.stem for path in files} == expected, 'Missing/extra source-coupling rank CSV')
    markers = sorted((directory / 'reference').glob('source-coupling-rank-*.complete'))
    require({path.stem for path in markers} == expected, 'Missing/extra source-coupling completion marker')
    metrics = {'ranks': ranks, 'steps_per_rank': steps, 'max_self_force_defect': 0.,
               'max_relative_fixed_point_work_defect': 0., 'max_iterations': 0}
    evidence = {}
    header = ['step', 'iterations', 'self_force_defect', 'relative_fixed_point_work_defect',
              'trial_pushes', 'committed_pushes']
    for path in files:
        marker = path.with_suffix('.complete')
        require(marker.read_text().split() == ['1', str(steps)], 'Wrong source-coupling completion step')
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle)
            require(reader.fieldnames == header, 'Unexpected source-coupling CSV layout')
            rows = list(reader)
        require(len(rows) == steps, 'Incomplete source-coupling step history')
        for step, row in enumerate(rows, 1):
            values = {key: float(row[key]) for key in header}
            require(all(math.isfinite(value) for value in values.values()), 'Nonfinite source-coupling history')
            require(values['step'] == step, 'Source-coupling step history is not consecutive')
            for key in ['step', 'iterations', 'trial_pushes', 'committed_pushes']:
                require(values[key] >= 0 and values[key].is_integer(), 'Invalid source-coupling integer counter')
            require(values['iterations'] >= 1, 'Zero converged source iterations')
            force = values['self_force_defect']
            work = values['relative_fixed_point_work_defect']
            require(0 <= force < 1e-13, 'Source fixed-point force did not converge')
            require(0 <= work < 1e-8, 'Source fixed-point work did not converge')
            metrics['max_self_force_defect'] = max(metrics['max_self_force_defect'], force)
            metrics['max_relative_fixed_point_work_defect'] = max(metrics['max_relative_fixed_point_work_defect'], work)
            metrics['max_iterations'] = max(metrics['max_iterations'], int(values['iterations']))
        evidence[path.name] = {'sha256': sha(path), 'marker_sha256': sha(marker)}
    return dict(metrics=metrics, files=evidence)


def detector_closure(directory, card, expected_particles):
    config = card['detectors']
    expected = {plane['name'] + '.h5' for plane in config['field_planes']}
    expected.update(plane['name'] + '-energy.h5' for plane in config['energy_planes'])
    expected.update(plane['name'] + '.h5' for plane in config['particle_planes'])
    expected.update(plane['name'] + '-substeps.h5' for plane in config['particle_planes'])
    require(len(expected) == 9 and config['substep_particle_planes'] is True,
            'Control must contain exactly nine detector files including three substep files')
    files = sorted((directory / 'detectors').glob('*.h5'))
    require({path.name for path in files} == expected, 'Missing/extra detector H5')
    evidence = {}
    for path in files:
        with h5py.File(path, 'r') as handle:
            require(len(handle) == 1, 'Unexpected detector root group count')
            group = handle[next(iter(handle))]
            require(int(group['complete'][()]) == 1, 'Unclosed detector H5')
            for marker, records in [('committed_records', 'records'), ('committed_samples', 'time_s')]:
                if marker in group:
                    require(int(group[marker][()]) == len(group[records]), 'Uncommitted detector records')
            if group.name.rsplit('/', 1)[-1] in ['particle_plane', 'particle_substep_plane']:
                require(len(group['records']) == expected_particles, 'Wrong particle crossing count')
            datasets = {}
            def inspect(name, dataset):
                if not isinstance(dataset, h5py.Dataset):
                    return
                if dataset.shape and not dataset.size:
                    datasets[name] = {'shape': list(dataset.shape), 'finite': True}
                    return
                stride = 8 if name == 'fields' else 65536
                for begin in range(0, len(dataset) if dataset.shape else 1, stride):
                    block = dataset[begin:begin+stride] if dataset.shape else dataset[()]
                    names = block.dtype.names
                    if names:
                        require(all(np.isfinite(block[field]).all() for field in names
                                    if block[field].dtype.kind in 'fc'), 'Nonfinite detector record')
                    elif block.dtype.kind in 'fc':
                        require(np.isfinite(block).all(), 'Nonfinite detector dataset')
                datasets[name] = {'shape': list(dataset.shape), 'finite': True}
            group.visititems(inspect)
            evidence[path.name] = {'complete': True, 'datasets': datasets}
    return evidence


def main(root, output):
    root = Path(root).resolve()
    output = Path(output).resolve() if output else root / 'INDEPENDENT_REVIEW.json'
    require(not output.exists(), 'Refusing to overwrite independent review')
    # Both collective normal shutdown and the closed-file parent acceptance
    # must exist before any H5 is opened here. This script never polls a job.
    completed = root / 'COMPLETE.json'
    require(completed.is_file(), 'Smoke has not published collective normal completion; no H5 inspected')
    result = json.loads(completed.read_text())
    require(result.get('program_integrity_pass') is True, 'Smoke program checks did not pass')
    require(result.get('physical_acceptance') is False, 'Unexpected physical acceptance claim')
    for key in ['linear_science_bitwise', 'notch_science_bitwise']:
        require(result.get(key) and all(value is True for value in result[key].values()),
                'Original numerical outputs are not bitwise preserved: ' + key)
    require(result.get('three_point_node_consistency'), 'Missing single/dual-node consistency acceptance')
    plan = json.loads((root / 'plan.json').read_text())
    build = Path(plan['build_record']) / 'COMPLETE.json'
    identity = json.loads(build.read_text())
    require(identity.get('passed') is True and sha(plan['binary']) == identity['binary_sha256'] ==
            result['binary_sha256'], 'Candidate binary identity mismatch')
    frozen = json.loads((root / 'FROZEN_INPUTS.json').read_text())
    require(all(sha(path) == digest for path, digest in frozen.items()), 'Smoke frozen inputs changed')
    reviews = {}
    for case in plan['cases']:
        if case['shape'] != 'linear-binomial3':
            continue
        directory = Path(case['directory'])
        local = json.loads((directory / 'RESULT.json').read_text())
        require(local.get('program_integrity_pass') is True, 'Case not closed')
        text = (directory / 'aprl.log').read_text()
        require(all(token in text for token in ['Configured stop reached:', 'phase=solve-complete',
                'phase=mpi-finalize-complete', '[source-coupling] linear-binomial3:',
                '[source-coupling] conservative-source;']), 'Missing normal stop or actual three-point startup')
        card = yaml.safe_load((directory / 'run.yaml').read_text())
        manifest_path = Path(card['output']['manifest'])
        manifest = yaml.safe_load(manifest_path.read_text())
        resolved = manifest['resolved_mesh']
        require(manifest['mpi_size'] == case['ranks'], 'Manifest rank count mismatch')
        require(Path(manifest['configuration_path']).resolve() == directory / 'run.yaml', 'Wrong manifest card')
        require(resolved['longitudinal_shape'] == card['mesh']['longitudinal_shape'] == 'linear-binomial3' and
                resolved['particle_coupling'] == card['mesh']['particle_coupling'] == 'conservative-source' and
                resolved['binomial3_coefficients'] == [.25, .5, .25] and
                resolved['binomial3_response'] == 'cos(kz*dz/2)^2', 'Actual manifest did not select matched three-point mode')
        metadata = (directory / 'reference' / 'metadata.txt').read_text().splitlines()
        require('longitudinal_source linear-binomial3-cos2-three-point-v1' in metadata,
                'Prepared reference source identity is not three-point')
        convergence = source_convergence(directory, case['ranks'], plan['expected_steps'])
        detectors = detector_closure(directory, card, plan['expected_particles'])
        reviews[case['name']] = dict(program_integrity_pass=True, physical_acceptance=False,
            convergence=convergence, detector_files=detectors, manifest_sha256=sha(manifest_path),
            reference_metadata_sha256=sha(directory / 'reference' / 'metadata.txt'),
            case_result_sha256=sha(directory / 'RESULT.json'))
    require(set(reviews) == {'binomial3-single', 'binomial3-dual'}, 'Missing requested three-point cases')
    report = dict(passed=True, physical_acceptance=False, reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
        scope='Independent closure/convergence/actual-mode gate; does not establish long-FEL physical accuracy.',
        binary_sha256=identity['binary_sha256'], parent_complete_sha256=sha(completed),
        frozen_input_files=len(frozen), original_linear_and_notch_science_bitwise=True,
        single_dual_science_comparison='Tolerance comparison in parent COMPLETE.json, not bitwise.', cases=reviews)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    main(args.root, args.output)
