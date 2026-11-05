"""Prepare four original short drift cards; build identity is sealed before launch."""
import argparse
import copy
import hashlib
import json
import shutil
from pathlib import Path

import yaml


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def prepare(workspace, build_record, output):
    workspace, output = Path(workspace).resolve(), Path(output).resolve()
    build_record = Path(build_record).resolve()
    output.mkdir(parents=True, exist_ok=False)
    scripts = output / 'scripts'
    scripts.mkdir()
    template = Path(__file__).resolve().parent
    for name in ['flow_prepare.py', 'flow_run.py', 'flow_analyze.py', 'flow.sbatch']:
        shutil.copy2(template / name, scripts / name)
    shutil.copy2(template / 'flow_prepare.py', scripts / 'prepare.py')
    cases_root = template.parent
    for folder, name in [('space_charge_physics', 'launcher.py'),
                         ('source_coupling_integration', 'binding.py'),
                         ('source_coupling_integration', 'analyze_closed.py')]:
        shutil.copy2(cases_root / folder / name, scripts / name)
    for name in ['report.py', 'review_discrete.py']:
        shutil.copy2(cases_root / 'discrete_energy_ledger' / name,
                     scripts / ('report_discrete.py' if name == 'report.py' else name))
    previous = workspace / 'results/longitudinal_main_flow_20261103/flow-r2'
    old = workspace / 'build/sparc-sase-r1-20261103/aprl'
    new = workspace / 'build/quintic-admission-current-20261105/aprl'
    require(old.is_file(), 'Missing baseline program')
    cases = []
    for shape in ['linear', 'quintic']:
        original = previous / ('drift-' + shape)
        require(json.loads((original / 'RESULT.json').read_text())['program_integrity_pass'],
                'Original short flow is not closed')
        for version, binary in [('old', old), ('new', new)]:
            name = shape + '-' + version
            directory = output / name
            directory.mkdir()
            card = yaml.safe_load((original / 'run.yaml').read_text())
            beam = copy.deepcopy(card['beam'])
            def rewrite(value):
                if isinstance(value, dict):
                    return {k: rewrite(v) for k, v in value.items()}
                if isinstance(value, list):
                    return [rewrite(v) for v in value]
                return value.replace(str(original), str(directory)) if isinstance(value, str) else value
            card = rewrite(card)
            card['beam'] = beam
            for folder in ['detectors', 'trajectory', 'ledger']:
                (directory / folder).mkdir()
            (directory / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
            cases.append(dict(name=name, shape=shape, version=version,
                              directory=str(directory), binary=str(binary),
                              ranks=16, threads=8, nodes=1, limit_seconds=300,
                              original_card=str(original / 'run.yaml')))
    plan = dict(workspace=str(workspace), new_simulation='Short regression only; no production',
                baseline_binary=str(old), baseline_binary_sha256=sha(old),
                candidate_binary=str(new), particles=31675, expected_steps=257,
                candidate_build_record=str(build_record),
                binding_probe=str(workspace / 'build/cluster-tools/binding_probe'), cases=cases,
                candidate_identity='Must pass build COMPLETE.json and freeze its binary hash before running')
    save(output / 'plan.json', plan)
    paths = [old, output / 'plan.json', *scripts.iterdir()]
    for case in cases:
        card_path = Path(case['directory']) / 'run.yaml'
        paths.extend([card_path, Path(yaml.safe_load(card_path.read_text())['beam']['input']['file'])])
    save(output / 'FROZEN_INPUTS.json', {str(p): sha(p) for p in paths if p.is_file()})
    save(output / 'PREPARED.json', dict(prepared=True, submitted=False,
                                      candidate_program_not_yet_frozen=True, cases=cases))
    print(output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('workspace', type=Path)
    p.add_argument('build_record', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    prepare(a.workspace, a.build_record, a.output)
