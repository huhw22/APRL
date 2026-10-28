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


def accepted_pending(path, source):
    """Reuse completed acceptance; this does not run a new preflight."""
    path = path.resolve()
    accepted = json.loads(path.read_text())
    assert accepted['pass'] and accepted['core_source_sha256']
    binary = Path(accepted['binary']).resolve()
    binary_sha = digest(binary)
    assert os.access(binary, os.X_OK) and binary_sha == accepted['binary_sha256']
    for name, expected in accepted['core_source_sha256'].items():
        assert digest(source / name) == expected, 'Tested core changed: ' + name
    evidence = {str(path): digest(path)}
    for name, expected in accepted['proof_sha256'].items():
        proof = path.parent / name
        assert digest(proof) == expected, 'Accepted proof changed: ' + name
        evidence[str(proof)] = expected
    for name in ('pair-r1/SCIENCE_PASS.json', 'mpi-r1/MPI_PASS.json',
                 'long-queue-r1/LONG_QUEUE_PASS.json', 'full-r2/TIMING_REVIEW.json'):
        assert name in accepted['proof_sha256'], 'Missing accepted stage: ' + name
        proof = json.loads((path.parent / name).read_text())
        sha = proof['binary_sha256']
        assert proof['pass'] and (sha.get('candidate') if isinstance(sha, dict) else sha) == binary_sha
    for name in ('unit-r1/PASS', 'unit-r1/required.log'):
        proof = path.parent / name
        evidence[str(proof)] = digest(proof)
    assert (path.parent / 'unit-r1/PASS').read_text().strip() == 'PASS'
    assert '100% tests passed, 0 tests failed out of 76' in (path.parent / 'unit-r1/required.log').read_text()
    manifest_path = path.parent / 'pair-r1/science/undulator-16x8-candidate/run-manifest.yaml'
    manifest = yaml.safe_load(manifest_path.read_text())
    compiled = manifest['source_revision']
    stamp = accepted['compiled_stamp'].split(';', 1)[0].strip()
    assert compiled.split('-dirty', 1)[0].startswith(stamp.split('-dirty', 1)[0])
    assert compiled.endswith('-dirty') == stamp.endswith('-dirty')
    evidence[str(manifest_path)] = digest(manifest_path)
    return {'binary': str(binary), 'binary_sha256': binary_sha,
            'compiled_revision': compiled, 'source_commit': accepted['core_commit']}, evidence


def reused_case(campaign, name, target, binary_sha):
    """Require the unchanged cancelled card and its completed input geometry."""
    old_root = campaign / name
    old_plan = json.loads((old_root / 'run-plan.json').read_text())
    frozen = json.loads((old_root / 'frozen.json').read_text())
    for filename, expected in frozen.items():
        assert digest(filename) == expected, 'Old frozen file changed: ' + filename
    passed = json.loads((campaign / 'preflight/PASS.json').read_text())
    assert passed['pass'] and passed['binary_sha256'] == old_plan['binary_sha256']
    assert passed['card_sha256'][str(old_root / 'run.yaml')] == digest(old_root / 'run.yaml')
    geometry_path = campaign / 'preflight' / (name + '-geometry.json')
    checked = json.loads(geometry_path.read_text())
    assert checked['pass'] and checked['input_unchanged']
    assert checked['records'] == old_plan['records'] == 41732775
    assert checked['input_sha256'] == passed['input_sha256'] == old_plan['input_sha256']
    inp = Path(old_plan['input_file']); stat = inp.stat()
    assert (stat.st_size, stat.st_mtime_ns) == (old_plan['input_size_bytes'], old_plan['input_mtime_ns'])
    geo = json.loads((old_root / 'geometry-plan.json').read_text())
    assert geo == old_plan['geometry'] and geo['source_input_sha256'] == old_plan['input_sha256']
    old = yaml.safe_load((old_root / 'run.yaml').read_text())
    card = copy.deepcopy(old)
    replacements = ((('output', 'manifest'), 'run-manifest.yaml'),
                    (('detectors', 'directory'), 'detectors'),
                    (('energy_ledger', 'directory'), 'ledger'),
                    (('mesh', 'moving_window', 'reference_directory'), 'reference'))
    for keys, relative in replacements:
        current = card
        for key in keys[:-1]: current = current[key]
        current[keys[-1]] = str(target / relative)
    normalized = copy.deepcopy(card)
    for keys, _ in replacements:
        current, previous = normalized, old
        for key in keys[:-1]: current, previous = current[key], previous[key]
        current[keys[-1]] = previous[keys[-1]]
    assert normalized == old, 'Reused card changed beyond output/reference paths'
    assert (old_plan['nodes'], old_plan['ranks'], old_plan['threads_per_rank']) == (1, 16, 8)
    assert binary_sha != old_plan['binary_sha256'], 'Reuse mode expects the accepted pending replacement'
    files = [old_root / f for f in ('run.yaml', 'run-plan.json', 'geometry-plan.json', 'frozen.json')]
    files += [campaign / 'preflight/PASS.json', geometry_path]
    return old_root, old, card, geo, {str(f): digest(f) for f in files}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--accepted-source', type=Path, help='Completed pending SOURCE_FINAL.json; requires --reuse-campaign')
    p.add_argument('--reuse-campaign', type=Path, help='Unchanged prior moving campaign with completed geometry preflight')
    args = p.parse_args()
    if bool(args.accepted_source) != bool(args.reuse_campaign):
        p.error('--accepted-source and --reuse-campaign must be provided together')
    w, root = args.workspace.resolve(), args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root / 'plan.json').exists(): raise ValueError('Campaign already prepared')
    source = w / 'source/aprl'
    evidence = {}
    if args.accepted_source:
        accepted, evidence = accepted_pending(args.accepted_source, source)
    else:
        accepted = json.loads((w / 'results/moving_window_integration_20261028/SOURCE_FINAL.json').read_text())
    binary = Path(accepted['binary'])
    binary_sha = digest(binary)
    assert os.access(binary, os.X_OK) and binary_sha == accepted['binary_sha256']
    for name, expected in accepted.get('source_files_sha256', {}).items():
        assert digest(source / name) == expected, 'Tested core changed: ' + name
    revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    scripts = root / 'scripts'; scripts.mkdir()
    for f in HERE.iterdir():
        if f.is_file(): shutil.copy2(f, scripts / f.name)
    cases = {}
    reuse_checks = {}
    for name, old_name, plane in (
        ('undulator_4p5m', 'elegant_field_only_20261027', 4.5),
        ('free_drift_0p5m', 'elegant_free_drift_20261027', .5)):
        parent = w / 'results' / old_name / 'production'
        target = root / name; target.mkdir()
        if args.reuse_campaign:
            parent, old, card, geo, case_evidence = reused_case(args.reuse_campaign.resolve(), name, target, binary_sha)
            evidence.update(case_evidence)
            reuse_checks[name] = {'source_card': str(parent / 'run.yaml'),
                'source_card_sha256': digest(parent / 'run.yaml'), 'normalized_card_equal': True,
                'geometry_pass_reused': str(args.reuse_campaign.resolve() / 'preflight' / (name + '-geometry.json'))}
        else:
            geo = derive(parent, ranks=16)
            old = yaml.safe_load((parent / 'run.yaml').read_text())
            card = copy.deepcopy(old)
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
        assert card['detectors']['field_planes'][0]['z'] == plane
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
            'expected_pending_queue': bool(args.accepted_source),
            'mpi_environment': {'OMPI_MCA_btl': 'self,vader,tcp', 'OMPI_MCA_btl_openib_warn_no_device_params_found': '0'},
            'physics': 'Original full input, release positions/momenta, magnet/field plane and fine resolution; fast blocked Cowan reference plus moving E/B window.',
            'changes': ['Compact moving E/B window with left reference and causal right crop.',
                'Accepted blocked Cowan reference replaces the long-domain Poisson/held field preparation.',
                'Only transverse four CPML faces; z guard cells remain ordinary physical cells.',
                'Explicit disk-backed reference domain covers all inlet queries; no complete long E/B donor.',
                'Hard duration uses the existing 1.1 ballistic planning bound; physical stop remains after-last-element.',
                'No checkpoint; no particle-plane output or extra particle-retirement layer.'],
            'limits': ['Full-scale weak-tail field accuracy remains to compare with the old completed pair.',
                'All particles remain in memory; 16x8 is a resource choice, not a measured optimum.',
                'This is an enabled moving-window experiment, not a default-binary replacement.']}
        if args.accepted_source:
            plan['changes'] = ['Same physical and geometric card as the prior cancelled moving run; output/reference paths retargeted.',
                'Accepted pending admission schedule avoids per-step scans of particles waiting outside the active window.',
                'Pending axial ledger statistics retain the accepted affine diagnostic approximation.',
                'Original per-job MPI self,vader,tcp configuration restored; no system MPI change.']
            plan['accepted_source_manifest'] = str(args.accepted_source.resolve())
            plan['acceptance_reuse'] = reuse_checks[name]
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
    if args.accepted_source:
        preflight = root / 'preflight'; preflight.mkdir()
        save(preflight / 'ACCEPTANCE_REUSE.json', {'pass': True, 'binary_sha256': binary_sha,
            'new_simulation_preflight_run': False, 'cases': reuse_checks, 'evidence_sha256': evidence,
            'scope': 'Existing completed acceptance and full-input geometry reused after identity/card/stat checks; no new scientific or geometry simulation was run.'})
        input_hashes = {json.loads((Path(path) / 'run-plan.json').read_text())['input_sha256'] for path in cases.values()}
        assert len(input_hashes) == 1
        save(preflight / 'PASS.json', {'pass': True, 'binary_sha256': binary_sha,
            'input_sha256': next(iter(input_hashes)),
            'card_sha256': {str(Path(path) / 'run.yaml'): digest(Path(path) / 'run.yaml') for path in cases.values()},
            'accepted_evidence_reuse': True, 'new_simulation_preflight_run': False,
            'acceptance_reuse': str(preflight / 'ACCEPTANCE_REUSE.json'),
            'scope': 'Reused accepted pending binary tests plus unchanged full-input geometry; not a newly executed preflight, full-physics comparison or speedup guarantee.'})
        for path in cases.values():
            frozen = Path(path) / 'frozen.json'; hashes = json.loads(frozen.read_text())
            hashes.update(evidence)
            for f in (preflight / 'PASS.json', preflight / 'ACCEPTANCE_REUSE.json'): hashes[str(f)] = digest(f)
            save(frozen, hashes)
    save(root / 'plan.json', {'format': 'aprl-moving-production-campaign-v1',
        'workspace': str(w), 'binary': str(binary), 'binary_sha256': binary_sha,
        'binding_probe': str(w / 'build/cluster-tools/binding_probe'),
        'nodes_per_case': 1, 'ranks': 16, 'threads': 8, 'cases': cases,
        'preparation_revision': revision, 'core_source_manifest_verified': True,
        'accepted_evidence_reuse': bool(args.accepted_source)})
    print(json.dumps(cases, indent=2))


if __name__ == '__main__': main()
