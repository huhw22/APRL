#!/usr/bin/env python3
"""Bounded offline audit of independent held-source reference cases; no physics pass gate."""
import argparse
import importlib.util
import json
from pathlib import Path
import struct

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args();root = args.root.resolve()
    helper_path = Path(__file__).with_name('oracle_review.py')
    spec = importlib.util.spec_from_file_location('oracle_review_helpers', helper_path)
    h = importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    require = h.require
    require((root/'COMPLETED.json').is_file(), 'Campaign is not complete')
    plan = json.loads((root/'plan.json').read_text())
    inputs = {Path(__file__).resolve(), helper_path, root/'COMPLETED.json', root/'plan.json',
              root/'source-files.json', root/'frozen.json'}
    cases = {}
    for item in plan['cases']:
        name = item['name'];directory = root/name
        metrics = json.loads((directory/'metrics.json').read_text())
        require(metrics['completed'] and metrics['mode'] == item['mode'] and
                metrics['beta_box'] == .23 and metrics['steps'] == 176 and
                metrics['cpml_cells'] == [0, 0, 20] and metrics['history_policy'] == 'world-position',
                'Unexpected fixture: '+name)
        initial, release, events = (h.table(directory/file, True) for file in h.CSV_NAMES[:3])
        plane = h.table(directory/'field-plane.csv')
        require(len(initial) == len(release) == len(events) == 1024 and len(plane) == 18000,
                'Incomplete cohort/field plane: '+name)
        summary, u, gamma, weights = h.particle_state(events, initial)
        radius = lambda t: float(np.sqrt(np.sum(weights*(t['x_m']**2+t['y_m']**2))))
        summary.update(initial_rms_radius_m=radius(initial), detector_rms_radius_m=radius(events))
        cases[name] = dict(metrics=metrics, initial=initial, release=release, events=events,
                           plane=plane, summary=summary, u=u, gamma=gamma, weights=weights)
        inputs.update(directory/file for file in h.CSV_NAMES)
        inputs.add(directory/'metrics.json')
    full_name = next(name for name, case in cases.items() if case['metrics']['mode'] == 'full')
    full = cases[full_name]
    donors = {case['metrics']['guard_cells']: name for name, case in cases.items()
              if case['metrics']['mode'] == 'held-record'}
    donor_reviews = []
    first_release = float(np.min(full['release']['step']+full['release']['fraction']))
    for guard, name in donors.items():
        case = cases[name];metrics = case['metrics']
        require(metrics['prescribed_reference_donor'] and not metrics['prescribed_actual_control'] and
                metrics['reference_provenance_kind'] == 'prescribed-uniform-initial-source', 'Wrong donor provenance')
        require(np.array_equal(case['initial'], full['initial']) and np.array_equal(case['release'], full['release']),
                'Donor initial or prescribed release differs from actual full')
        for key in ('x_m', 'y_m', 'ux', 'uy', 'uz', 'charge_C'):
            require(np.array_equal(case['events'][key], case['initial'][key]), 'Donor changed seed '+key)
        require(np.all(case['events']['delta_energy_box_eV'] == 0) and
                np.all(case['events']['delta_energy_lab_eV'] == 0), 'Donor gained particle energy')
        mask = full['plane']['step'] < first_release
        for key in full['plane'].dtype.names:
            require(np.array_equal(case['plane'][key][mask].view(np.uint64),
                                   full['plane'][key][mask].view(np.uint64)),
                    'Donor and actual full differ before first release: '+key)
        directory = root/name/'native-oracle'
        manifest = json.loads((directory/'manifest.json').read_text())
        require(manifest['completed'] and manifest['source_kind'] == 1 and manifest['frames'] == 177,
                'Donor manifest is not completed source kind 1')
        with (directory/'header.bin').open('rb') as stream:
            require(stream.read(8) == b'APRLOH02', 'Donor is not V2')
            integers = struct.unpack('<8Q', stream.read(64))
            require(integers == (24, 24, 424, 176, guard, 710271, 6, 1), 'Donor source kind/config mismatch')
        inputs.update(directory/file for file in ('header.bin', 'manifest.json', 'complete.bin'))
        donor_reviews.append({'case': name, 'source_kind': 1, 'V2_header_verified': True,
                              'seed_xy_momentum_charge_preserved': True, 'particle_delta_energy_zero': True,
                              'release_events_identical_to_actual_full': True,
                              'pre_release_field_samples_bitwise_equal': int(mask.sum()),
                              'last_pre_release_sample_step': int(full['plane']['step'][mask].max()),
                              'first_release_continuous_step': first_release,
                              'overall_wall_s': metrics['overall_wall_s'], 'scope': 'Full-grid prescribed-source donor remains O(L^2)'})
    comparisons = []
    for name, case in cases.items():
        mode = case['metrics']['mode']
        if mode not in ('residual-held-reference', 'residual-held-zero', 'residual'):
            continue
        zero = mode == 'residual-held-zero'
        base_name = donors[case['metrics']['guard_cells']] if zero else full_name
        base = cases[base_name];metrics = case['metrics']
        require(np.array_equal(case['initial'], base['initial']) and np.array_equal(case['release'], base['release']),
                'Candidate changed initial/prescribed release states')
        if mode.startswith('residual-held'):
            require(metrics['prescribed_reference_input'] and metrics['oracle_feed_channel'] and metrics['oracle_cpml_channel'] and
                    metrics['reference_provenance_kind'] == 'prescribed-uniform-initial-source' and
                    Path(metrics['oracle_directory']).resolve() == (root/donors[metrics['guard_cells']]/'native-oracle').resolve(),
                    'Candidate consumed the wrong reference source')
            require(metrics['prescribed_actual_control'] is zero and not metrics['prescribed_reference_donor'],
                    'Actual PIC/control mode confusion')
        nonzero_energy = bool(np.any(case['events']['delta_energy_lab_eV'] != 0))
        momentum_changed = bool(np.any(case['u'] != h.vector(case['initial'], ('ux', 'uy', 'uz'))))
        if zero:
            require(not nonzero_energy and not momentum_changed, 'Zero-response control is not prescribed')
        else:
            require(nonzero_energy and momentum_changed and
                    case['summary']['detector_rms_radius_m'] != case['summary']['initial_rms_radius_m'],
                    'Actual candidate accidentally retained prescribed particle motion')
        d = case['u']-base['u']
        dg = np.sum(d*(case['u']+base['u']), axis=1)/(case['gamma']+base['gamma'])
        energy_error = h.BOOST*(dg+h.BOOST_BETA*d[:, 2])*h.REST
        comparisons.append({'case': name, 'reference_case': base_name, 'zero_response_control_only': zero,
                            'actual_particle_energy_nonzero': nonzero_energy, 'actual_momentum_changed': momentum_changed,
                            'fields': h.plane_comparison(base['plane'], case['plane'], base['events'], case['events']),
                            'particles': {**case['summary'],
                                'mean_delta_difference_eV': case['summary']['mean_delta_energy_lab_eV']-base['summary']['mean_delta_energy_lab_eV'],
                                'sigma_difference_eV': case['summary']['sigma_energy_lab_eV']-base['summary']['sigma_energy_lab_eV'],
                                'particle_energy_rms_error_eV': float(np.sqrt(np.sum(case['weights']*energy_error**2))),
                                'particle_energy_max_abs_error_eV': float(np.max(np.abs(energy_error)))},
                            'cpml_history_peak_abs': metrics['cpml_history_peak_abs'],
                            'candidate_wall_s': metrics['overall_wall_s']})
    replay = []
    pairs = []
    donor_names = list(donors.values())
    pairs += [(root/donor_names[0], root/name) for name in donor_names[1:]]
    if 'b023-residual-g24' in cases:
        pairs.append((root/'b023-residual-g24', root.parent/'r1/b023-residual-g24'))
    pairs.append((root/full_name, Path(plan['legacy_replay'][full_name]['directory'])))
    for a, b in pairs:
        hashes = {}
        for file in h.CSV_NAMES:
            require((a/file).read_bytes() == (b/file).read_bytes(), 'Scientific CSV mismatch: '+str(a/file)+' vs '+str(b/file))
            hashes[file] = h.sha256(a/file);inputs.update((a/file, b/file))
        replay.append({'left': str(a), 'right': str(b), 'all_five_csv_byte_identical': True, 'sha256': hashes})
    current_sources = json.loads((root/'source-files.json').read_text())
    prior_path = root.parent/'r1/source-files.json';inputs.add(prior_path)
    prior_sources = json.loads(prior_path.read_text())
    core = {path: sha for path, sha in current_sources.items() if '/source/aprl/src/' in path}
    require(core and all(prior_sources.get(path) == sha for path, sha in core.items()), 'Core source changed from r1')
    checks = {}
    for folder in ('checks', 'oracle-checks'):
        path = root/folder/'metrics.json';inputs.add(path)
        checks[folder] = json.loads(path.read_text())
        require(checks[folder]['completed'], 'Incomplete operator/recording checks')
    result = {'format': 'aprl-independent-held-reference-review-v1', 'campaign': str(root),
              'data_contract_checks': 'PASS', 'physical_acceptance': 'Not assigned; no physical threshold imposed',
              'binary_sha256': plan['binary_sha256'], 'unchanged_core_source_files_vs_r1': len(core),
              'build_provenance_scope': 'Use SOURCE_FINAL.json for r4/r5 metadata-only rebuild provenance; this review does not claim r4 source snapshots exactly identify the compiled metadata strings.',
              'full_particle_moments': full['summary'], 'donors': donor_reviews, 'comparisons': comparisons,
              'csv_replay': replay, 'operator_checks': checks,
              'input_sha256': {str(path): h.sha256(path) for path in sorted(inputs)},
              'limitations': ['The independent held donor uses only seed trajectories, but its full-domain evolution retains O(L^2) cost.',
                              'Donor generation plus candidate cost must both be counted; no speedup or production claim.',
                              'Zero-response results compare only with the held donor and do not validate physical particle motion.',
                              'F-Fref is not established to be pure radiation; MPI, long-time and production integration remain unaccepted.',
                              'Native E(n), B(n-1/2) comparison; before-last-detection includes earlier outgoing particles.']}
    output = root/'HELD_REVIEW.json'
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(output)
    for row in comparisons:
        whole = row['fields']['intervals']['whole']
        print(row['case'], 'vs', row['reference_case'], 'E', whole['E']['relative_l2'],
              'B', whole['B']['relative_l2'], 'particle_RMS_eV', row['particles']['particle_energy_rms_error_eV'])


if __name__ == '__main__':
    main()
