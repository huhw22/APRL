#!/usr/bin/env python3
"""Check checkpoint-disabled cost against the accepted three-optimization build."""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil

import yaml

WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
CONTROL_BUILD = Path('build/compute-followup-20261026')
CONTROL_SHA = '5edd096eec8d02fef03613a4f74731d90eba5574cb6c109fe26c7fb4f16a9483'
CONTROL_REVISION = '8adce330f646701123ef5b312b914b37152af874'
CASES = ('grid128', 'particles262144')
LAYOUT = (4, 4, 'compact')
REPEATS = 5
HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def helper_path():
    frozen = HERE / 'fixed_driver.py'
    return frozen if frozen.exists() else HERE.parent / 'performance/cases/fixed_mpi_followup/run.py'


fixed = module('checkpoint_disabled_fixed_helpers', helper_path())
fixed.CASE_NAMES = CASES
fixed.REPEATS = REPEATS
memory, save, digest = fixed.memory, fixed.save, fixed.digest


def checkpoint_evidence(pass_file, trial, control):
    """A timing campaign may start only after exact-binary restart correctness."""
    directory = pass_file.parent
    paths = [pass_file, directory / 'plan.json', directory / 'RESULT.json', directory / 'frozen.json']
    if not all(path.is_file() for path in paths):
        raise ValueError('Missing completed checkpoint correctness evidence')
    plan = json.loads(paths[1].read_text())
    result = json.loads(paths[2].read_text())
    if plan.get('format') != 'aprl-checkpoint-acceptance-v1' or result.get('status') != 'PASS':
        raise ValueError('Checkpoint correctness did not pass')
    if Path(plan['binary']).resolve() != trial or plan['binary_sha256'] != digest(trial):
        raise ValueError('Checkpoint correctness used a different candidate executable')
    if Path(plan['baseline']).resolve() != control or plan['baseline_sha256'] != CONTROL_SHA:
        raise ValueError('Checkpoint correctness used a different accepted baseline')
    if (directory / 'FAIL').exists() or (directory / 'FAILED').exists():
        raise ValueError('Checkpoint correctness contains a failure marker')
    memory.verify_frozen(directory)
    return dict(pass_path=str(pass_file), files={str(path): digest(path) for path in paths},
                trial_sha256=digest(trial), control_sha256=CONTROL_SHA)


def prepare(args):
    workspace, root = args.workspace.resolve(), args.root.resolve()
    control, trial = args.control.resolve(), args.trial.resolve()
    if root.exists() or control == trial:
        raise ValueError('Use a fresh result directory and a separate candidate executable')
    if digest(control) != CONTROL_SHA:
        raise ValueError('Control must be the accepted 8adce33 three-optimization executable')
    revisions = {label: fixed.compiled_revision_evidence(binary)
                 for label, binary in (('control', control), ('trial', trial))}
    if revisions['control']['compiled_revision'] != CONTROL_REVISION or \
            revisions['trial']['compiled_revision'] != args.trial_revision:
        raise ValueError('Compiled revision does not match the requested full SHA')
    builds = {label: fixed.build_metadata(binary)
              for label, binary in (('control', control), ('trial', trial))}
    if builds['control']['selected'] != builds['trial']['selected']:
        raise ValueError('Control/candidate compiler, flags, MPI or HDF5 toolchain differs')
    if builds['control']['selected']['APRL_ENABLE_OPENMP'] != 'ON' or \
            builds['control']['selected']['APRL_REQUIRE_PARALLEL_HDF5'] != 'ON':
        raise ValueError('OpenMP and required parallel HDF5 must be enabled')
    required = fixed.required_evidence(args.required_pass.resolve(), trial)
    correctness = checkpoint_evidence(args.checkpoint_pass.resolve(), trial, control)
    source = args.source_campaign.resolve()
    if not (source / 'PASS').is_file():
        raise ValueError('The source workload campaign must have passed')
    source_plan = json.loads((source / 'plan.json').read_text())
    source_frozen = json.loads((source / 'frozen.json').read_text())
    expected = {item['name']: item for item in source_plan['cases']}
    original_files = [source / 'cards' / (case + '.yaml') for case in CASES]
    original_files += [source / 'scripts' / name
                       for name in ('generic_inspector.py', 'basic_inspector.py')]
    for file in original_files:
        if digest(file) != source_frozen[str(file)]:
            raise ValueError('Source card or inspector changed: ' + str(file))
    for case in CASES:
        file = source / 'cards' / (case + '.yaml')
        card = yaml.safe_load(file.read_text())
        if digest(file) != expected[case]['sha256'] or 'checkpoint' in card:
            raise ValueError('Input must be the unchanged checkpoint-disabled card: ' + case)
        if card['beam']['input']['type'] != 'generated-gaussian':
            raise ValueError('This benchmark must use the generated-gaussian cards')
    runtime_dir = control.parent / 'runtime-libs'
    runtime_files = {str(runtime_dir / name): sha
                     for name, sha in fixed.CONTROL_RUNTIME_LIBS.items()}
    for file, sha in runtime_files.items():
        if digest(file) != sha:
            raise ValueError('Accepted runtime library changed: ' + file)
    probe = workspace / 'build/node-survey/binding_probe'
    if not probe.is_file():
        raise ValueError('The binding probe is missing')
    root.mkdir(parents=True)
    scripts, cards = root / 'scripts', root / 'cards'
    scripts.mkdir(); cards.mkdir()
    for name in ('run.py', 'run.sh', 'README.md'):
        shutil.copy2(HERE / name, scripts / name)
    shutil.copy2(helper_path(), scripts / 'fixed_driver.py')
    shutil.copy2(fixed.helper_path(), scripts / 'memory_driver.py')
    for name in ('generic_inspector.py', 'basic_inspector.py'):
        shutil.copy2(source / 'scripts' / name, scripts / name)
    cases = []
    for case in CASES:
        copied = cards / (case + '.yaml')
        shutil.copy2(source / 'cards' / copied.name, copied)
        record = dict(expected[case]); record.update(card=str(copied), sha256=digest(copied))
        cases.append(record)
    binaries = {label: dict(path=str(binary), sha256=digest(binary),
                           compiled_revision=revisions[label]['compiled_revision'], build=builds[label])
                for label, binary in (('control', control), ('trial', trial))}
    plan = dict(format='aprl-checkpoint-disabled-performance-v1', root=str(root), workspace=str(workspace),
                prepared_at=datetime.now(timezone.utc).isoformat(), source_campaign=str(source),
                binaries=binaries, compiled_revisions=revisions, required_regression=required,
                checkpoint_correctness=correctness,
                runtime_libraries=dict(directory=str(runtime_dir), files=runtime_files),
                binding_probe=str(probe), binding_probe_sha256=digest(probe),
                layout=dict(ranks=4, threads=4, placement='compact', cpus=16),
                expected_host='gnode1', physical_core_budget=16, smt_workers=False,
                cases=cases, repeats=REPEATS, order_seed=fixed.SEED,
                run_timeout_s=180, total_timeout_s=7200, primary_throughput_interval=[20, 80],
                profile_interval=0, checkpoint_enabled=False, expected_runs=28,
                limits=['Same card and fixed 4x4 for each comparison; different cases are not speedups.',
                        'Four pilots, four warmups, twenty throughput runs; every run has profiling off.',
                        'Control is the accepted combined 8ad build; production remains untouched.',
                        'PASS means complete runs and strict scientific equivalence, not zero-cost proof.',
                        'Report paired timings and uncertainty; this does not measure snapshot write bandwidth.'])
    save(root / 'plan.json', plan)
    frozen = [root / 'plan.json', *sorted(scripts.iterdir()), *sorted(cards.iterdir()), control, trial, probe]
    frozen += [Path(item['flags']) for item in revisions.values()]
    frozen += [Path(item['cache']) for item in builds.values()]
    frozen += [Path(path) for path in runtime_files]
    frozen += [Path(path) for path in required['files']]
    frozen += [Path(path) for path in correctness['files']]
    save(root / 'frozen.json', {str(path): digest(path) for path in frozen})
    print(json.dumps(dict(prepared=str(root), layout='4 MPI x 4 threads', cases=list(CASES),
                          repeats=REPEATS, profile_interval=0, expected_runs=28), indent=2))


class Campaign(fixed.Campaign):
    def __init__(self, root):
        super().__init__(root)
        if self.plan['format'] != 'aprl-checkpoint-disabled-performance-v1' or \
                set(self.cases) != set(CASES) or self.plan['profile_interval'] != 0:
            raise ValueError('Unexpected disabled-performance plan')
        self.report['format'] = 'aprl-checkpoint-disabled-performance-result-v1'

    def one(self, case_name, variant, phase, repeat=0):
        if phase not in ('pilot', 'warmup', 'throughput'):
            raise ValueError('Profiling and snapshot-writing phases are forbidden in this campaign')
        if 'checkpoint' in yaml.safe_load(Path(self.cases[case_name]['card']).read_text()):
            raise ValueError('Checkpoint configuration appeared in a throughput card')
        super().one(case_name, variant, phase, repeat)
        row = self.report['runs'][-1]
        directory = Path(row['directory'])
        if row['profile_interval'] != 0 or '[checkpoint]' in (directory / 'run.log').read_text():
            row.update(valid=False, status='failed'); self.checkpoint()
            raise ValueError('Checkpoint activity or profiling appeared in a disabled run')
        if (directory / 'checkpoints').exists():
            row.update(valid=False, status='failed'); self.checkpoint()
            raise ValueError('Disabled run unexpectedly created a checkpoint directory')

    def pilot(self):
        required = fixed.required_evidence(Path(self.plan['required_regression']['pass_path']),
                                           Path(self.plan['binaries']['trial']['path']))
        if required != self.plan['required_regression']:
            raise ValueError('Required regression evidence changed')
        checked = checkpoint_evidence(Path(self.plan['checkpoint_correctness']['pass_path']),
                                      Path(self.plan['binaries']['trial']['path']),
                                      Path(self.plan['binaries']['control']['path']))
        if checked != self.plan['checkpoint_correctness']:
            raise ValueError('Checkpoint correctness evidence changed')
        if (self.root / 'PILOT_PASS').exists():
            return
        self.binding()
        for case in CASES:
            for variant in ('control', 'trial'):
                self.one(case, variant, 'pilot')
        self.report.update(pilot_pass=True, status='pilot_passed'); self.checkpoint()
        (self.root / 'PILOT_PASS').write_text('Fixed 4x4 binding and both disabled-case paired pilots passed.\n')

    def matrix(self):
        if (self.root / 'PASS').exists():
            return
        if not (self.root / 'PILOT_PASS').exists():
            raise ValueError('Pilot must pass first')
        self.report['status'] = 'running_throughput'; self.checkpoint()
        for case, variant in fixed.paired_order(CASES, 0):
            self.one(case, variant, 'warmup')
        for repeat in range(1, REPEATS + 1):
            for case, variant in fixed.paired_order(CASES, repeat):
                self.one(case, variant, 'throughput', repeat)
        self.summarize()
        if len(self.report['runs']) != 28 or not all(row['valid'] for row in self.report['runs']):
            raise ValueError('Incomplete disabled-performance campaign')
        self.report.update(status='passed', pass_=True); self.checkpoint()
        (self.root / 'PASS').write_text('All 28 disabled runs passed scientific checks; consult paired timing summary for overhead.\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=WORKSPACE)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', choices=('prepare', 'pilot', 'matrix', 'all'), required=True)
    parser.add_argument('--control', type=Path, default=WORKSPACE / CONTROL_BUILD / 'aprl')
    parser.add_argument('--trial', type=Path)
    parser.add_argument('--trial-revision')
    parser.add_argument('--required-pass', type=Path)
    parser.add_argument('--checkpoint-pass', type=Path)
    parser.add_argument('--source-campaign', type=Path,
                        default=WORKSPACE / 'results/memory_scaling/campaign-20261025')
    args = parser.parse_args()
    if args.stage == 'prepare':
        if any(value is None for value in (args.trial, args.trial_revision, args.required_pass, args.checkpoint_pass)):
            parser.error('prepare requires --trial, --trial-revision, --required-pass and --checkpoint-pass')
        prepare(args); return 0
    if platform.node() != 'gnode1' or os.environ.get('APRL_ALLOW_UNSCHEDULED_CPU_TESTS') != '1' or os.environ.get('SLURM_JOB_ID'):
        parser.error('Run only on the authorized gnode1 CPU host, without fabricated Slurm variables')
    root = args.root.resolve(); memory.verify_frozen(root)
    campaign = Campaign(root)
    if (root / 'PASS').exists() or (args.stage == 'pilot' and (root / 'PILOT_PASS').exists()):
        print('Requested stage already complete; no duplicate launches.'); return 0
    if campaign.report['status'] == 'failed':
        raise ValueError('Prepare a new result directory after a failed campaign')
    if not campaign.report['runs']:
        campaign.report['initial_resources'] = fixed.idle_host_snapshot(); campaign.checkpoint()
    try:
        campaign.report['status'] = 'running'; campaign.checkpoint()
        if args.stage in ('pilot', 'all'):
            campaign.pilot()
        if args.stage in ('matrix', 'all'):
            campaign.matrix()
        memory.verify_frozen(root)
        return 0
    except Exception as error:
        campaign.report.update(status='failed', pass_=False, error=str(error))
        campaign.summarize(); campaign.checkpoint()
        (root / 'FAILED').write_text(str(error) + '\n')
        print('FAILED:', error, flush=True)
        return 1
    finally:
        campaign.checkpoint()


if __name__ == '__main__':
    raise SystemExit(main())
