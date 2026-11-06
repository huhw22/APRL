"""Prepare only the completed-smoke-design FEL off/on cards on two nodes.

This reads frozen cards and plans, never running HDF5. It does not modify or
restart the source campaign and never launches a job.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import yaml


WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def rewrite(value, original, destination):
    if isinstance(value, dict):
        return {key: rewrite(item, original, destination) for key, item in value.items()}
    if isinstance(value, list):
        return [rewrite(item, original, destination) for item in value]
    return value.replace(str(original), str(destination)) if isinstance(value, str) else value


def main(workspace, output, binary, build_record, source_campaign=None, nodes=2, ranks=32):
    workspace, output, binary, build_record = [Path(value).resolve()
        for value in [workspace, output, binary, build_record]]
    source_campaign = (Path(source_campaign).resolve() if source_campaign else
                       workspace / 'results/mpi_particle_balance_20261106/smoke-r1')
    assert nodes == 2 and ranks > 0 and ranks % nodes == 0
    parent = json.loads((source_campaign / 'plan.json').read_text())
    selected = {case['mode']:case for case in parent['cases']
                if case['family'] == 'fel7-w64' and case['mode'] in ['off', 'on']}
    assert set(selected) == {'off', 'on'}
    # Validate every input and source card before creating destination cases.
    source_inputs = {}
    for case in selected.values():
        card_path = Path(case['directory']) / 'run.yaml'
        assert sha(card_path) == case['card_sha256']
        source = Path(case['input'])
        source_inputs[str(source)] = sha(source)
        assert source_inputs[str(source)] == case['input_sha256']
        assert case['expected_particles'] == 178725 and case['expected_steps'] == 379
        assert case['active_grid'] == [80, 80, 342]
        assert 342 >= 8 * ranks
    output.mkdir(parents=True, exist_ok=False)
    cases, cards, canonical_hashes = [], {}, []
    for mode in ['off', 'on']:
        old = selected[mode]
        original = Path(old['directory'])
        name = 'fel7-w64-' + mode
        directory = output / name
        directory.mkdir()
        base = yaml.safe_load((original / 'run.yaml').read_text())
        card = rewrite(copy.deepcopy(base), original, directory)
        card['beam'] = copy.deepcopy(base['beam'])
        assert card['mesh']['particle_coupling'] == 'conservative-source'
        assert card['mesh']['longitudinal_shape'] == 'linear-binomial3'
        assert bool(card['mesh'].get('mpi_particle_balance', {}).get('enabled', False)) == (mode == 'on')
        for folder in ['detectors', 'ledger']:
            (directory / folder).mkdir()
        if card.get('trajectory', {}).get('enabled'):
            (directory / 'trajectory').mkdir()
        assert not (directory / 'reference').exists()
        path = directory / 'run.yaml'
        path.write_text(yaml.safe_dump(card, sort_keys=False))
        digest = sha(path)
        cards[str(path)] = digest
        canonical = rewrite(copy.deepcopy(card), directory, '@CASE')
        canonical['mesh'].pop('mpi_particle_balance', None)
        canonical_sha = hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()
        canonical_hashes.append(canonical_sha)
        case = copy.deepcopy(old)
        case.update(name=name, family='dual-fel7-w64', mode=mode, directory=str(directory),
            binary=str(binary), build_record=str(build_record), nodes=nodes, ranks=ranks,
            threads=1, limit_seconds=600, case_limit_seconds=600, profile_interval=0,
            card_sha256=digest, normalized_physical_card_sha256=canonical_sha,
            copied_card=str(original / 'run.yaml'), copied_card_sha256=old['card_sha256'])
        save(directory / 'run-plan.json', case)
        cases.append(case)
    assert len(set(canonical_hashes)) == 1
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(build_record), accepted_for_execution=False, cases=cases,
        nodes=nodes, ranks=ranks, threads=1, recommended_slurm_minutes=30, profile_interval=0,
        source_inputs=source_inputs, frozen_cards=cards,
        source_campaign=str(source_campaign), source_campaign_plan_sha256=sha(source_campaign / 'plan.json'),
        source_cards={str(Path(case['directory']) / 'run.yaml'):case['card_sha256']
                      for case in selected.values()},
        physical_acceptance=False, submitted=False,
        scope='Two-node single-thread off/on ownership regression, using exactly the smoke-design seven-period w64 physical cards. No old program rerun.',
        unchanged=['Source events/weights/charge', 'All source-card physical settings and diagnostics',
                   'Nx/Ny80/reference1024/active342', 'Spacing/dt/boost/stop379 steps'],
        controlled_changes=['Separate output directories', 'r3 candidate program identity',
                            'Two-node32MPIx1 resource layout shared by both versions'])
    save(output / 'plan.json', plan)
    save(output / 'FROZEN_CARDS.json', cards)
    save(output / 'PREPARATION_COMPLETE.json', dict(complete=True, submitted=False,
        accepted_for_execution=False, cases=2, nodes=nodes, ranks=ranks, threads=1,
        source_campaign_unchanged=True, source_inputs_unchanged=True,
        physical_cards_identical_except_balance=True, reference_directories_absent=True))
    (output / 'README.md').write_text(
        '# Two-node MPI ownership regression\n\n'
        'Prepared only. The two cases copy the original smoke-design fel7-w64 off/on '
        'cards, preserving input, geometry, timestep,379-step stop and every diagnostic. '
        'Only output paths, r3 program identity and the common two-node32MPIx1 layout '
        'change. No historical program is rerun. The source campaign is not modified, '
        'and no running HDF5 is read.600s per case and30minutes total are proposed.\n')
    print(json.dumps(dict(complete=True, output=str(output), cases=2,
                          nodes=nodes, ranks=ranks, expected_steps=379), indent=2))
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--out', '--output', dest='output', type=Path,
                        default=WORKSPACE / 'results/mpi_particle_balance_20261106/two-nodes-r1')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--build-record', '--buildrecord', dest='build_record', type=Path, required=True)
    parser.add_argument('--source-campaign', type=Path)
    parser.add_argument('--nodes', type=int, default=2)
    parser.add_argument('--ranks', type=int, default=32)
    args = parser.parse_args()
    main(args.workspace, args.output, args.binary, args.build_record,
         args.source_campaign, args.nodes, args.ranks)
