#!/usr/bin/env python3
"""Freeze one bounded, diagnostic-only replay of the accepted free-drift case."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import yaml


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--field-stages', action='store_true')
    a = p.parse_args()
    ws, build, root = a.workspace.resolve(), a.build.resolve(), a.root.resolve()
    root.mkdir()  # Never reuse an executed campaign.
    baseline = ws/'results/energy_transport_audit_20261029/campaign-r3/science/free-drift-r1-new-on'
    prior = ws/'results/energy_transport_audit_20261029/campaign-r3/plan.json'
    oldplan = json.loads(prior.read_text())
    assert json.loads((baseline/'CASE_REVIEW.json').read_text())['checks_pass']
    case = root/'free-drift'; case.mkdir()
    (case/'mesh-work').mkdir()
    if a.field_stages: (case/'field-stages').mkdir()
    card = yaml.safe_load((baseline/'run.yaml').read_text())
    assert not card['sources']['magnetic_elements'] and not card['sources']['incident_waves']
    assert card['mesh']['current_filter'] == 'none'
    card['output']['manifest'] = str(case/'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory'] = str(case/'reference')
    for name in ('detectors', 'trajectory', 'energy_ledger'):
        card[name]['directory'] = str(case/name if name != 'energy_ledger' else case/'ledger')
    (case/'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
    scripts = root/'scripts'; scripts.mkdir()
    here = Path(__file__).resolve().parent
    for name in ('run.py', 'run.sbatch', 'prepare_run.py', 'deposition_work_probe.cpp', 'grid_work_observer_checks.cpp'):
        shutil.copy2(here/name, scripts/name)
    transport = scripts/'transport'; transport.mkdir()
    for name in ('review.py', 'prepare.py', 'run.py'):
        shutil.copy2(here.parent/'energy_transport_audit'/name, transport/name)
    plan = {'root': str(root), 'case': str(case), 'baseline': str(baseline),
            'build': str(build), 'binary': str(build/'aprl'),
            'binding_probe': oldplan['binding_probe'], 'particles': oldplan['particles'],
            'probe': str(build/'deposition_work_probe'), 'steps': 257,
            'observer_checks': str(build/'grid_work_observer_checks'),
            'physics_changed': False, 'purpose': 'Independently measure grid J dot E; free drift only',
            'binary_sha256': sha(build/'aprl'), 'limits_s': {'probe': 90, 'case': 600},
            'source_parent': '7d130d8', 'baseline_build': str(ws/'build/energy-transport-audit-r2-20261029')}
    plan['field_stages']=a.field_stages
    save(root/'plan.json', plan)
    frozen = [root/'plan.json', case/'run.yaml', Path(card['beam']['input']['file']),
              Path(plan['binding_probe']), build/'aprl', build/'deposition_work_probe', build/'grid_work_observer_checks',
              baseline/'run.yaml', baseline/'CASE_REVIEW.json', prior]
    frozen += [x for x in scripts.rglob('*') if x.is_file()]
    save(root/'frozen.json', {str(x): sha(x) for x in frozen})
    print(root)


if __name__ == '__main__':
    main()
