#!/usr/bin/env python3
"""Prepare a trial-only continuation with an immutable copy of a checked control."""
import argparse
import json
from pathlib import Path
import shutil

import run


def prepare(root, previous, workspace):
    old = json.loads((previous / 'plan.json').read_text())
    submission = json.loads((previous / 'SUBMITTED.json').read_text())
    run.require(submission['nodes'] == ['node1', 'node8'], 'Reuse requires the original node pair')
    run.require(not (previous / 'trial/aprl.log').exists(), 'Trial already started; do not duplicate it')
    run.require((previous / 'JOB_EXIT_CODE').is_file(), 'Original allocation must have ended')
    for line in (previous / 'frozen.sha256').read_text().splitlines():
        sha, file = line.split('  ', 1)
        run.require(run.digest(file) == sha, 'Original frozen evidence changed: ' + file)
    run.prepare(root, workspace)
    plan = json.loads((root / 'plan.json').read_text())
    for key in ('variants', 'steps', 'dt_box_s', 'expected_mesh', 'input_sha256', 'ranks', 'threads_per_rank'):
        run.require(plan[key] == old[key], 'Recovery comparison changed: ' + key)
    # Copy completed outputs, never change or rerun the historical control.
    # Its configuration/manifest retain the original output paths as evidence.
    shutil.copytree(previous / 'control', root / 'control', dirs_exist_ok=True)
    plan['total_limit_s'] = 5400
    plan['reused_control'] = {
        'directory': str(previous / 'control'), 'job_id': submission['job_id'],
        'nodes': submission['nodes'], 'method': 'independent copy of closed completed-600-step outputs',
        'control_rerun': False, 'original_failed_check_preserved': str(previous / 'control/check.log')}
    plan['limitations'].append('The same-node control/trial pair spans two allocations after a checker-only failure; it is not a repeated measurement.')
    run.save(root / 'plan.json', plan)
    # Validate the reusable evidence before allocating any new CPU time.
    # Copy only the completed small checks for this offline revalidation;
    # the recovery allocation will independently run its new small cards.
    for name in ('preflight-result.json', 'PREFLIGHT_PASS'):
        shutil.copy2(previous / name, root / name)
    run.check_run(root, 'control')
    (root / 'PREFLIGHT_PASS').unlink()
    (root / 'preflight-result.json').rename(root / 'prior-preflight-result.json')
    (root / 'control/check.log').rename(root / 'control/original-failed-check.log')
    frozen = sorted(p for p in root.rglob('*') if p.is_file() and p.name != 'frozen.sha256')
    frozen += [Path(v['binary']) for v in plan['variants'].values()]
    frozen += [Path(v['toolchain']['cache']) for v in plan['variants'].values()]
    frozen += [workspace / 'build/cluster-tools/binding_probe']
    frozen += [Path(p) for p in plan['required_evidence']]
    (root / 'frozen.sha256').write_text(''.join(f'{run.digest(p)}  {p}\n' for p in frozen))
    print(root)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--previous', required=True, type=Path)
    parser.add_argument('--workspace', default=run.WORKSPACE, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.previous.resolve(), args.workspace.resolve())
