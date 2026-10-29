#!/usr/bin/env python3
"""Prepare/finalize a bounded energy-transport audit; never submit jobs.

Create: --workspace WORKSPACE --root NEW_ROOT [--candidate NEW_BINARY]
Freeze after the binary is ready: --root ROOT --finalize --candidate NEW_BINARY
An unfrozen plan is reviewable but cannot be executed.
"""
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
BASELINE_SHA = '490e4e6663418efcdce0bb84fb8c77248fbe0f7c38a9c8d220a5d4874a124fd1'
VARIANTS = ('old-on', 'new-on', 'new-off', 'new-audit-on')


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
    """Only optional diagnostics and output identities may differ within a pair."""
    result = copy.deepcopy(card)
    result['beam']['input']['file'] = '<INPUT>'
    result['output']['manifest'] = '<OUTPUT>'
    result['mesh']['moving_window']['reference_directory'] = '<REFERENCE>'
    for name in ('detectors', 'trajectory', 'energy_ledger'):
        result[name]['directory'] = '<OUTPUT>'
    result['energy_ledger']['enabled'] = '<ON_OR_OFF>'
    result['energy_ledger'].pop('discrete_operator_audit', None)
    result['detectors'].pop('energy_planes', None)
    return result


def finalize(root, candidate):
    plan = json.loads((root / 'plan.json').read_text())
    require(not plan['frozen'] and not (root / 'frozen.json').exists(), 'Plan already frozen')
    require(candidate is not None, '--candidate is required to finalize')
    candidate = candidate.resolve()
    plan['binaries']['new'] = str(candidate)
    require(sha(plan['binaries']['old']) == BASELINE_SHA, 'Accepted v4 baseline changed')
    frozen = [root / 'plan.json', Path(plan['input_file']), Path(plan['binding_probe'])]
    frozen += [p for p in (root / 'scripts').iterdir() if p.is_file()]
    frozen += [Path(case['directory']) / 'run.yaml'
               for group in plan['cases'].values() for case in group]
    for case in plan['cases']['science']:
        if case.get('reused_from'):
            frozen += [Path(path) for path in case['reused_evidence_sha256']]
            frozen += [Path(case['reused_from'])/relative for relative in case['reused_output_sha256']]
    for name, value in plan['binaries'].items():
        binary = Path(value)
        require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable: ' + value)
        require((binary.parent / 'runtime-libs').is_dir(), 'Missing frozen runtime-libs: ' + value)
        plan['binary_sha256'][name] = sha(binary)
        frozen += [binary] + [p for p in (binary.parent / 'runtime-libs').iterdir() if p.is_file()]
    plan['source_head_at_freeze'] = subprocess.check_output(
        ['git', '-C', str(Path(plan['workspace']) / 'source/aprl'), 'rev-parse', 'HEAD'], text=True).strip()
    plan['frozen'] = True
    save(root / 'plan.json', plan)
    save(root / 'frozen.json', {str(p): sha(p) for p in frozen})
    print(root / 'frozen.json')



def prepare_reuse(plan, old_root):
    import review
    old_plan = json.loads((old_root/'plan.json').read_text())
    require(old_plan.get('frozen') is True, 'Reuse source must be a frozen campaign')
    require(plan['binaries']['new'] is not None, 'Reuse requires an explicit candidate binary')
    expected = {key:sha(path) for key,path in plan['binaries'].items()}
    require(expected == old_plan['binary_sha256'], 'Reuse source uses different binaries')
    require(plan['input_sha256'] == old_plan['input_sha256'], 'Reuse input contents differ')
    binding_path = old_root/'science/execution/BINDING.json'
    binding = json.loads(binding_path.read_text())
    require(binding['pass'] and binding['ranks']==16 and binding['threads']==8 and binding['physical_cores']==128,
            'Reuse source lacks fixed16x8 binding evidence')
    require(len(binding['nodes'])==1, 'Only completed one-node science cases may be reused')
    reused = []
    for item in plan['cases']['science']:
        if item['variant'] not in ('old-on','new-on','new-off'):
            continue
        candidates = [old for old in old_plan['cases']['science'] if
            old['magnet']==item['magnet'] and old['variant']==item['variant'] and old['repeat']==item['repeat']]
        if not candidates:
            continue
        require(len(candidates)==1, 'Ambiguous reused case')
        old_item = candidates[0]; directory = Path(old_item['directory'])
        report_path = directory/'CASE_REVIEW.json'
        if not report_path.is_file():
            continue
        old_report = json.loads(report_path.read_text())
        require(old_report['checks_pass'], 'Cannot reuse a failed case')
        for key in ('binary','ledger','energy_planes','magnet','repeat','nodes','ranks','threads','steps'):
            require(old_item[key]==item[key], 'Reuse case differs: '+key)
        require(canonical(yaml.safe_load((directory/'run.yaml').read_text())) ==
                canonical(yaml.safe_load((Path(item['directory'])/'run.yaml').read_text())),
                'Reuse physics card differs')
        old_card = yaml.safe_load((directory/'run.yaml').read_text())
        require(sha(old_card['beam']['input']['file']) == plan['input_sha256'], 'Reused input changed')
        inspected = review.inspect(directory, old_item, old_plan, export_csv=False)
        require(inspected['checks_pass'], 'Reused closed output failed independent inspection')
        item['reused_from'] = str(directory)
        item['reused_campaign_root'] = str(old_root)
        item['reused_original_nodes'] = binding['nodes']
        item['reused_original_binding'] = binding
        evidence = [old_root/'plan.json', old_root/'frozen.json', binding_path,
                    report_path, directory/'run.yaml', directory/'aprl.log', directory/'aprl.log.launch.json']
        item['reused_evidence_sha256'] = {str(path):sha(path) for path in evidence}
        item['reused_output_sha256'] = {relative:sha(directory/relative) for relative in inspected['files']}
        item['reused_binary_sha256'] = expected
        reused.append(item['name'])
    require(reused, 'No eligible completed science cases to reuse')
    plan['reused_science'] = {'source_root':str(old_root), 'cases':reused,
        'scope':'Same binary/input/physics/layout and closed numeric outputs verified; original node recorded. No performance samples reused.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--finalize', action='store_true')
    parser.add_argument('--plan-only', action='store_true', help='Keep plan unfrozen even when candidate exists')
    parser.add_argument('--reuse-science', type=Path, help='Read-only reuse of complete same-binary science cases')
    args = parser.parse_args()
    root = args.root.resolve()
    if args.finalize:
        finalize(root, args.candidate)
        return
    require(args.workspace is not None, '--workspace is required when creating a plan')
    workspace = args.workspace.resolve()
    require(not root.exists(), 'Use a new campaign directory')
    require(args.repeats == 3, 'This bounded campaign declares exactly three repetitions')
    baseline = workspace / 'build/moving-energy-ledger-r2-20261029/aprl'
    require(sha(baseline) == BASELINE_SHA, 'Accepted v4 baseline changed')
    source = workspace / 'results/moving_window_integration_20261028/dense-two-node-r3'
    original = yaml.safe_load((source / 'moving/run.yaml').read_text())
    require(original['mesh']['cells'] == [160, 160, 896], 'Source fixture mesh changed')
    require(original['boundary']['cells'] == [10, 10, 0], 'Source fixture CPML changed')
    moving = original['mesh']['moving_window']
    require(moving['lower_guard_cells'] == 64 and moving['upper_guard_cells'] == 128,
            'Source guard geometry changed')
    require(original['energy_ledger']['sample_interval_steps'] == 1, 'Science cadence changed')
    cadence = original['detectors']['field_planes'][0]['rhythm']
    planes = [dict(name='audit-entry', z=-0.05, rhythm=cadence, buffer_samples=16, compression=0),
              dict(name='audit-exit', z=0.06, rhythm=cadence, buffer_samples=16, compression=0)]
    root.mkdir(parents=True)
    scripts, inputs = root / 'scripts', root / 'inputs'
    scripts.mkdir(); inputs.mkdir()
    for path in HERE.iterdir():
        if path.is_file():
            shutil.copy2(path, scripts / path.name)
    shutil.copy2(source / 'inputs/held.h5', inputs / 'held.h5')
    cases = {'science': [], 'performance': [], 'mpi': []}
    for stage in cases:
        repeats = args.repeats if stage == 'performance' else 1
        for repeat in range(repeats):
            for magnet in ('undulator', 'free-drift'):
                order = list(VARIANTS if stage != 'mpi' else ('new-on', 'new-audit-on'))
                if repeat % 2:
                    order.reverse()
                for variant in order:
                    name = f'{magnet}-r{repeat + 1}-{variant}'
                    directory = root / stage / name
                    directory.mkdir(parents=True)
                    card = copy.deepcopy(original)
                    card['beam']['input']['file'] = str(inputs / 'held.h5')
                    card['output']['manifest'] = str(directory / 'run-manifest.yaml')
                    card['mesh']['moving_window']['reference_directory'] = str(directory / 'reference')
                    for section, folder in (('detectors', 'detectors'), ('trajectory', 'trajectory'), ('energy_ledger', 'ledger')):
                        card[section]['directory'] = str(directory / folder)
                    ledger = variant != 'new-off'
                    audit = variant == 'new-audit-on'
                    card['energy_ledger']['enabled'] = ledger
                    # Omit the new key entirely for the accepted old parser.
                    if audit:
                        card['detectors']['energy_planes'] = copy.deepcopy(planes)
                    detector_names = [entry['name'] for key in ('field_planes','particle_planes','energy_planes')
                                      for entry in card['detectors'].get(key,[])]
                    require(len(detector_names)==len(set(detector_names)), 'Detector names must be globally unique')
                    if magnet == 'free-drift':
                        card['sources']['magnetic_elements'] = []
                    if stage == 'performance':
                        card['energy_ledger']['sample_interval_steps'] = 10
                        card['trajectory']['enabled'] = False
                        card['detectors']['particle_planes'] = []
                        card['detectors']['field_planes'][0]['rhythm'] *= 10
                        # Energy planes retain full cadence. This measures the
                        # real incremental cost of the requested online audit.
                    card['runtime']['resource_monitor']['progress_interval_steps'] = 20
                    (directory / 'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
                    cases[stage].append(dict(name=name, directory=str(directory), variant=variant,
                        binary='old' if variant.startswith('old-') else 'new', ledger=ledger,
                        energy_planes=audit, magnet=magnet, repeat=repeat + 1,
                        nodes=2 if stage == 'mpi' else 1, ranks=16, threads=8, steps=257,
                        timeout_s=180, profile_interval=0))
        for magnet in ('undulator', 'free-drift'):
            selected = [case for case in cases[stage] if case['magnet'] == magnet]
            reference = canonical(yaml.safe_load((Path(selected[0]['directory']) / 'run.yaml').read_text()))
            require(all(canonical(yaml.safe_load((Path(case['directory']) / 'run.yaml').read_text())) == reference
                        for case in selected), 'Physics card differs within stage')
    plan = dict(format='aprl-energy-transport-audit-v1', frozen=False, workspace=str(workspace),
        binaries={'old': str(baseline), 'new': str(args.candidate.resolve()) if args.candidate else None},
        binary_sha256={'old': BASELINE_SHA},
        source_head_at_prepare=subprocess.check_output(['git', '-C', str(workspace / 'source/aprl'),
            'rev-parse', 'HEAD'], text=True).strip(),
        binding_probe=str(workspace / 'build/cluster-tools/binding_probe'),
        input_file=str(inputs / 'held.h5'), input_sha256=sha(inputs / 'held.h5'), particles=31675,
        source_card=str(source / 'moving/run.yaml'), source_card_sha256=sha(source / 'moving/run.yaml'),
        cases=cases, repeats=args.repeats, energy_planes=planes,
        energy_plane_files={'entry':'audit-entry-energy.h5','exit':'audit-exit-energy.h5'},
        stage_budget_s={'science': 1500, 'performance': 2700, 'mpi': 840},
        mpi_environment={'OMPI_MCA_btl': 'self,vader,tcp', 'OMPI_MCA_btl_openib_warn_no_device_params_found': '0'},
        acceptance={
            'physical_output': 'Exact legacy E/B, time vectors, ID-sorted particle planes and per-rank trajectories within the same layout; audit files checked separately.',
            'energy': 'Independent recorded-term and pusher-work identities; report physical residual without residual-derived work or invented closure threshold.',
            'lifecycle': '257 physical steps, complete/committed output, all 31675 particles released/exited; all audit planes cover the same cohort.',
            'performance': 'Profile off; three repeats; 40-240 step throughput; all pairs and extrema retained.'},
        limits=[
            'This campaign validates diagnostics, not full-production radiation or space-charge convergence.',
            'All field-plane quantities are total self-consistent E/B, not isolated near fields.',
            'Boost moving-boundary and fixed-lab-plane budgets have separate frames and sampling surfaces.',
            'No cross-node or cross-stage timing comparison is used as implementation speedup.',
            'No new BF physics fixture; no default program, original input or production-card mutation.',
            'Input cohort is a short production-like fixture, not the original 41-million-particle beam.'])
    if args.reuse_science:
        prepare_reuse(plan, args.reuse_science.resolve())
    save(root / 'plan.json', plan)
    print(root / 'plan.json')
    if args.candidate and args.candidate.is_file() and not args.plan_only:
        finalize(root, args.candidate)


if __name__ == '__main__':
    main()
