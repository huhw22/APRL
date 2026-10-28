#!/usr/bin/env python3
"""Freeze the two authorized full-input moving-field-window production cards."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import yaml
from geometry import derive

HERE = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''): h.update(b)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    args = p.parse_args()
    w, root = args.workspace.resolve(), args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root / 'plan.json').exists(): raise ValueError('Campaign already prepared')
    source = w / 'source/aprl'
    accepted = json.loads((w / 'results/moving_window_integration_20261028/SOURCE_FINAL.json').read_text())
    binary = Path(accepted['binary'])
    binary_sha = digest(binary)
    assert os.access(binary, os.X_OK) and binary_sha == accepted['binary_sha256']
    for name, expected in accepted['source_files_sha256'].items():
        assert digest(source / name) == expected, 'Tested core changed: ' + name
    revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    scripts = root / 'scripts'; scripts.mkdir()
    for f in HERE.iterdir():
        if f.is_file(): shutil.copy2(f, scripts / f.name)
    cases = {}
    for name, old_name, plane in (
        ('undulator_4p5m', 'elegant_field_only_20261027', 4.5),
        ('free_drift_0p5m', 'elegant_free_drift_20261027', .5)):
        parent = w / 'results' / old_name / 'production'
        target = root / name; target.mkdir()
        geo = derive(parent, ranks=16)
        old = yaml.safe_load((parent / 'run.yaml').read_text())
        card = copy.deepcopy(old)
        assert card['detectors']['field_planes'][0]['z'] == plane
        card['output'] = {'overwrite': False, 'manifest': str(target / 'run-manifest.yaml')}
        card['mesh']['cells'] = [160, 160, geo['reference_cells_z']]
        card['mesh'].pop('longitudinal', None)
        card['mesh']['center'] = [0., 0., 0.]
        card['mesh']['duration'] = geo['duration_box_s']
        card['mesh']['moving_window'] = {'enabled': True, 'lower_guard_cells': 64,
            'upper_guard_cells': 128, 'block_layers': 32,
            'reference_directory': str(target / 'reference')}
        card['boundary']['cells'] = [10, 10, 0]
        card['initial_self_field']['model'] = 'cowan-relativistic-poisson'
        card.pop('checkpoint', None)
        card['detectors']['directory'] = str(target / 'detectors')
        card['energy_ledger']['directory'] = str(target / 'ledger')
        assert card['beam'] == old['beam'] and card['sources'] == old['sources']
        assert card['stop'] == old['stop'] == {'mode': 'after-last-element'}
        assert not card['detectors']['particle_planes']
        assert not card['particle_retirement']['enabled']
        assert card['detectors']['field_planes'] == old['detectors']['field_planes']
        assert card['mesh']['cell_size'] == old['mesh']['cell_size']
        inp = Path(card['beam']['input']['file']).resolve(); st = inp.stat()
        plan = {'format': 'aprl-moving-production-v1', 'case': name,
            'binary': str(binary), 'binary_sha256': binary_sha,
            'compiled_revision': accepted['compiled_revision'], 'core_commit': accepted['source_commit'],
            'preparation_revision': revision, 'source_run': str(parent),
            'source_config_sha256': digest(parent / 'run.yaml'),
            'nodes': 1, 'ranks': 16, 'threads_per_rank': 8,
            'input_file': str(inp), 'input_size_bytes': st.st_size,
            'input_mtime_ns': st.st_mtime_ns, 'input_sha256': geo['source_input_sha256'],
            'records': geo['source_records'], 'field_plane_z_m': plane,
            'geometry': geo, 'campaign': str(root),
            'preflight_pass': str(root / 'preflight/PASS.json'),
            'physics': 'Original full input, release positions/momenta, magnet/field plane and fine resolution; fast blocked Cowan reference plus moving E/B window.',
            'changes': ['Compact moving E/B window with left reference and causal right crop.',
                'Accepted blocked Cowan reference replaces the long-domain Poisson/held field preparation.',
                'Only transverse four CPML faces; z guard cells remain ordinary physical cells.',
                'Explicit disk-backed reference domain covers all inlet queries; no complete long E/B donor.',
                'Hard duration uses the existing 1.1 ballistic planning bound; physical stop remains after-last-element.',
                'No checkpoint; no particle-plane output or extra particle-retirement layer.'],
            'limits': ['First full-scale run of this moving-window version; weak-tail field accuracy remains to compare.',
                'All pending particles remain on the inlet rank; 16x8 is a resource choice, not a measured optimum.',
                'This is an enabled moving-window experiment, not a default-binary replacement.']}
        (target / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
        save(target / 'run-plan.json', plan)
        save(target / 'geometry-plan.json', geo)
        shutil.copy2(scripts / 'run.sbatch', target / 'run.sbatch')
        frozen = [target / 'run.yaml', target / 'run-plan.json', target / 'geometry-plan.json', target / 'run.sbatch']
        frozen += [f for f in scripts.iterdir() if f.is_file()]
        frozen += [binary, w / 'build/cluster-tools/binding_probe']
        frozen += [f for f in (binary.parent / 'runtime-libs').iterdir() if f.is_file()]
        save(target / 'frozen.json', {str(f): digest(f) for f in frozen})
        cases[name] = str(target)
    save(root / 'plan.json', {'format': 'aprl-moving-production-campaign-v1',
        'workspace': str(w), 'binary': str(binary), 'binary_sha256': binary_sha,
        'binding_probe': str(w / 'build/cluster-tools/binding_probe'),
        'nodes_per_case': 1, 'ranks': 16, 'threads': 8, 'cases': cases,
        'preparation_revision': revision, 'core_source_manifest_verified': True})
    print(json.dumps(cases, indent=2))


if __name__ == '__main__': main()
