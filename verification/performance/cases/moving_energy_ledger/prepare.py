#!/usr/bin/env python3
"""Freeze small moving-volume ledger validation; never submit or run MPI."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
BASELINE_SHA = 'cde740f2426fd87c4e0605759610f7d71423fe03d86925866037d8e0f2b7d116'
VARIANTS = ('old-on', 'new-on', 'new-off', 'old-off')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def canonical(card):
    result = copy.deepcopy(card)
    result['beam']['input']['file'] = '<INPUT>'
    result['output']['manifest'] = '<OUTPUT>'
    result['mesh']['moving_window']['reference_directory'] = '<REFERENCE>'
    for name in ('detectors', 'trajectory', 'energy_ledger'):
        result[name]['directory'] = '<OUTPUT>'
    result['energy_ledger']['enabled'] = '<ON_OR_OFF>'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    workspace, root = args.workspace.resolve(), args.root.resolve()
    require(not root.exists(), 'Use a new campaign directory')
    require(3 <= args.repeats <= 5, 'Use 3 to 5 predeclared performance repetitions')
    binaries = {'old': workspace/'build/pending-admission-20261028/aprl', 'new': args.candidate.resolve()}
    require(sha(binaries['old']) == BASELINE_SHA, 'Accepted pending-admission baseline changed')
    for binary in binaries.values():
        require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable: ' + str(binary))
        require((binary.parent/'runtime-libs').is_dir(), 'Need frozen runtime-libs beside binary')
    source = workspace/'results/moving_window_integration_20261028/dense-two-node-r3'
    original = yaml.safe_load((source/'moving/run.yaml').read_text())
    require(original['mesh']['cells'] == [160, 160, 896], 'Source fixture mesh changed')
    require(original['boundary']['cells'] == [10, 10, 0], 'Source fixture CPML changed')
    require(original['mesh']['moving_window']['lower_guard_cells'] == 64 and
            original['mesh']['moving_window']['upper_guard_cells'] == 128, 'Source guards changed')
    require(original['energy_ledger']['sample_interval_steps'] == 1, 'Science ledger cadence changed')
    root.mkdir(parents=True)
    scripts, inputs = root/'scripts', root/'inputs'
    scripts.mkdir(); inputs.mkdir()
    for path in HERE.iterdir():
        if path.is_file():
            shutil.copy2(path, scripts/path.name)
    shutil.copy2(source/'inputs/held.h5', inputs/'held.h5')
    cases = {'science': [], 'performance': [], 'mpi': []}
    for stage in cases:
        repetitions = args.repeats if stage == 'performance' else 1
        for repeat in range(repetitions):
            for magnet in ('undulator', 'free-drift'):
                order = list(VARIANTS if stage != 'mpi' else ('old-on', 'new-on'))
                # Pair order alternates; no winner-dependent retry or pruning.
                if repeat % 2:
                    order.reverse()
                for variant in order:
                    name = f'{magnet}-r{repeat + 1}-{variant}'
                    directory = root/stage/name
                    directory.mkdir(parents=True)
                    card = copy.deepcopy(original)
                    card['beam']['input']['file'] = str(inputs/'held.h5')
                    card['output']['manifest'] = str(directory/'run-manifest.yaml')
                    card['mesh']['moving_window']['reference_directory'] = str(directory/'reference')
                    for section, folder in (('detectors', 'detectors'), ('trajectory', 'trajectory'), ('energy_ledger', 'ledger')):
                        card[section]['directory'] = str(directory/folder)
                    card['energy_ledger']['enabled'] = variant.endswith('-on')
                    if magnet == 'free-drift':
                        card['sources']['magnetic_elements'] = []
                    if stage == 'performance':
                        card['energy_ledger']['sample_interval_steps'] = 10
                        card['trajectory']['enabled'] = False
                        card['detectors']['particle_planes'] = []
                        card['detectors']['field_planes'][0]['rhythm'] *= 10
                    card['runtime']['resource_monitor']['progress_interval_steps'] = 20
                    (directory/'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
                    cases[stage].append({'name': name, 'directory': str(directory), 'variant': variant,
                        'binary': variant.split('-')[0], 'ledger': variant.endswith('-on'), 'magnet': magnet,
                        'repeat': repeat + 1, 'nodes': 2 if stage == 'mpi' else 1, 'ranks': 16,
                        'threads': 8, 'steps': 257, 'timeout_s': 180, 'profile_interval': 0})
        # Diagnose only ledger enable/disable: all other same-stage card fields match.
        for magnet in ('undulator', 'free-drift'):
            selected = [x for x in cases[stage] if x['magnet'] == magnet]
            first = canonical(yaml.safe_load((Path(selected[0]['directory'])/'run.yaml').read_text()))
            require(all(canonical(yaml.safe_load((Path(x['directory'])/'run.yaml').read_text())) == first
                        for x in selected), 'Physical card mismatch within stage')
    plan = {'format': 'aprl-moving-energy-ledger-v1', 'workspace': str(workspace),
        'binaries': {name: str(path) for name, path in binaries.items()},
        'binary_sha256': {name: sha(path) for name, path in binaries.items()},
        'source_head': subprocess.check_output(['git', '-C', str(workspace/'source/aprl'), 'rev-parse', 'HEAD'], text=True).strip(),
        'binding_probe': str(workspace/'build/cluster-tools/binding_probe'),
        'input_file': str(inputs/'held.h5'), 'input_sha256': sha(inputs/'held.h5'), 'particles': 31675,
        'source_card': str(source/'moving/run.yaml'), 'source_card_sha256': sha(source/'moving/run.yaml'),
        'cases': cases, 'repeats': args.repeats, 'stage_budget_s': {'science': 1500, 'performance': 2700, 'mpi': 840},
        'mpi_environment': {'OMPI_MCA_btl': 'self,vader,tcp', 'OMPI_MCA_btl_openib_warn_no_device_params_found': '0'},
        'acceptance': {'physical_output': 'Same layout: exact numeric bytes for fields, time vectors, sorted particle planes and per-rank trajectory records.',
            'energy_balance': 'Independently recompute recorded-term identities; report physical residual and component curves without assigning an invented physics tolerance.',
            'lifecycle': '257 complete steps; all 31675 particles released and exited; closed complete=1 output required.',
            'performance': 'Profile disabled; 40-240 step interval, repeat medians/ranges and every pair retained; no performance PASS threshold.'},
        'limits': ['Small one-period source, not production energy/radiation convergence.',
            'Performance cards differ from science only in documented diagnostic cadence/outputs; no cross-stage speed ratio.',
            'One-node and two-node timings are never mixed.', 'No default binary, physics kernel, or production card changes.']}
    save(root/'plan.json', plan)
    frozen = [root/'plan.json', inputs/'held.h5', Path(plan['binding_probe'])]
    frozen += [p for p in scripts.iterdir() if p.is_file()]
    frozen += [Path(x['directory'])/'run.yaml' for group in cases.values() for x in group]
    for binary in binaries.values():
        frozen += [binary] + [p for p in (binary.parent/'runtime-libs').iterdir() if p.is_file()]
    save(root/'frozen.json', {str(path): sha(path) for path in frozen})
    print(root/'plan.json')


if __name__ == '__main__':
    main()
