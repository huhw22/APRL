"""Validate closed duration-bounded prefixes without promoting physical completion.

The exit contract is preserved from pending_admission_optimization/run.py.
Its check_duration_exit helper is frozen from full_scale_followup/run.py as
duration_review.py. Every HDF5 read follows the all-rank post-close proof.
"""
import csv
import json
import math
from pathlib import Path
import re

import h5py
import numpy as np
import yaml

from common import require
import report_discrete


def bounded_duration_status(text, record, steps, duration_helper):
    """Classify only the expected post-close MPI abort, never physical success.

    Simulation::solve places the duration-readiness MPI_Barrier after every
    rank's local output closes. A root *complete* trace proves all peers entered
    that barrier; sixteen trace lines are neither enabled nor required. main
    reports the expected exception and calls MPI_Abort(...,1), after which srun
    --kill-on-bad-exit can emit a generic STEP CANCELLED and kill peer ranks.
    """
    status = duration_helper.check_duration_exit(text, record['exit_code'], steps=steps)
    require(math.isfinite(record['seconds']) and 0 < record['seconds'] < record['limit_s'],
            'Launcher reached its timeout or has invalid elapsed time')
    phases = ('trajectory-close-complete', 'detector-close-complete', 'ledger-close-complete',
              'particle-boundary-report-complete', 'final-resource-report-complete',
              'duration-exhaustion-readiness-sync-enter', 'duration-exhaustion-readiness-sync-complete')
    markers = []
    for phase in phases:
        found = list(re.finditer(r'\[teardown\] rank=0\b[^\n]* phase=' + phase + r'\b', text))
        require(len(found) == 1, 'Missing/duplicate root close/readiness proof: ' + phase)
        markers.append(found[0].start())
    require(markers == sorted(markers), 'Root close/readiness phases are out of order')
    fatal = r'DUE TO TIME LIMIT|Timed out|Segmentation fault|Bus error|Floating point exception|Killed by signal|Signal:\s*\w+|oom[_ -]?kill|out.of.memory|\[shutdown-watchdog\]|Direct E/B run stopped after a complete field step'
    require(not re.search(fatal, text, re.I), 'Unexpected timeout/signal/crash/watchdog termination')
    require('Configured stop reached:' not in text, 'Bounded timing run unexpectedly claims physical completion')
    final = re.findall(r'\[resource\] phase=final step=(\d+)\b', text)
    require(final == [str(steps)], 'Final complete-step report is missing/duplicate/wrong')
    expected = duration_helper.EXHAUSTED
    failure = re.search(r'Simulation failure on MPI rank \d+:\s*' + re.escape(expected), text)
    aborts = list(re.finditer(r'MPI_ABORT was invoked on rank (\d+) in communicator MPI_COMM_WORLD\s+with errorcode (\d+)\.', text))
    require(aborts and all(match[2] == '1' for match in aborts), 'Missing expected MPI_Abort(errorcode=1) evidence')
    require(failure is not None, 'Missing complete expected duration exception')
    ready = markers[-1]
    expected_cancel = re.compile(r'slurmstepd: error: \*\*\* STEP (\d+\.\d+) ON (\S+) CANCELLED AT [0-9T:+.\-]+ \*\*\*')
    expected_exit = re.compile(r'srun: error: \S+: tasks? [0-9,\-]+: (?:Killed|Exited with exit code 1)')
    cancellations, peer_cleanup = [], []
    offset = 0
    for line in text.splitlines(keepends=True):
        stripped = line.strip()
        if re.search(r'CANCELLED', line, re.I):
            match = expected_cancel.fullmatch(stripped)
            require(match is not None and offset > ready and offset > failure.start() and
                    any(offset > abort.end() for abort in aborts),
                    'Cancellation is not the generic STEP cleanup after proven duration abort')
            cancellations.append({'step_id': match[1], 'host': match[2], 'line': stripped})
        elif stripped.startswith(('slurmstepd: error:', 'srun: error:')):
            require(expected_exit.fullmatch(stripped) is not None and offset > ready and
                    any(offset > abort.end() for abort in aborts),
                    'Unexpected scheduler/runtime error: ' + stripped)
            peer_cleanup.append(stripped)
        offset += len(line)
    status.update(readiness_complete_marker_observed=True,
        all_rank_readiness_proof='Rank-zero completion of the collective placed after every local output close; trace is root-only by default.',
        mpi_abort_errorcode=1, scheduler_step_cancellations_after_duration_abort=cancellations,
        scheduler_peer_cleanup_lines=peer_cleanup, physical_simulation_complete=False,
        contract='Post-close expected duration abort only; HDF5 committed extents/finite values are independently required by full_review.')
    return status


def duration_status(directory, case, runtime=None):
    require(case.get('bounded_duration_prefix') is True, 'Duration prefix was not explicitly authorized by its card plan')
    directory = Path(directory)
    require(not (directory / 'PIC_COMPLETE.json').exists(), 'A duration prefix must not have a PIC_COMPLETE marker')
    runtime = runtime if runtime is not None else json.loads((directory / 'aprl.log.launch.json').read_text())
    import duration_review
    return bounded_duration_status((directory / 'aprl.log').read_text(), runtime,
                                    case['expected_steps'], duration_review)


def load_prefix_ledger(directory, case):
    """The unchanged v5/discrete loader, with explicit expected physical flag0."""
    directory = Path(directory)
    duration_status(directory, case)
    with h5py.File(directory / 'ledger/energy-ledger.h5', 'r') as handle:
        group = handle['energy_ledger']
        rows = group['records'][:]
        require(int(group.attrs['format_version']) == 5, 'Need unchanged v5 physical ledger')
        require(int(group['complete'][()]) == 0 and
                int(group['committed_records'][()]) == len(rows), 'Prefix ledger transaction/physical flag differs')
        ranks = int(group.attrs['mpi_size'])
        frame = group.attrs['frame']
    require(b'boosted computational frame' in frame if isinstance(frame, bytes) else
            'boosted computational frame' in frame, 'Wrong ledger frame')
    require(len(rows) and int(rows['step'][-1]) == case['expected_steps'] and ranks == case['ranks'],
            'Prefix ledger final step/rank differs')
    paths = sorted((directory / 'ledger').glob('discrete-field-rank-*.csv'))
    require([path.name for path in paths] ==
            [f'discrete-field-rank-{rank:05d}.csv' for rank in range(ranks)], 'Missing/extra discrete rank files')
    steps = int(rows['step'][-1])
    sums = {}
    for path in paths:
        require(path.with_suffix('.complete').read_text().split() == ['0', str(steps)],
                'Discrete prefix history must close with physical_complete0 and exact steps')
        with path.open() as stream:
            reader = csv.DictReader(stream)
            records = list(reader)
            names = reader.fieldnames
        require([int(row['step']) for row in records] == list(range(1, steps + 1)),
                'Missing/repeated discrete prefix step')
        values = {name:np.array([report_discrete.LD(row[name]) for row in records])
                  for name in names if name != 'step'}
        require(all(np.isfinite(value).all() for value in values.values()), 'Nonfinite discrete prefix')
        report_discrete.validate_columns(values)
        for name, value in values.items():
            if name not in sums:
                sums[name] = np.zeros_like(value)
            sums[name] = np.maximum(sums[name], value) if name.startswith('max_') or '_prediction_scale_' in name else sums[name] + value
    return rows, sums, ranks, steps


def finite(data, where):
    if data.dtype.names:
        require(all(np.isfinite(data[name]).all() for name in data.dtype.names
                    if data[name].dtype.kind in 'fc'), 'Nonfinite prefix data: ' + str(where))
    elif data.dtype.kind in 'fc':
        require(np.isfinite(data).all(), 'Nonfinite prefix data: ' + str(where))


def verify_prefix_case(case, runtime=None):
    directory = Path(case['directory'])
    status = duration_status(directory, case, runtime)
    card = yaml.safe_load((directory / 'run.yaml').read_text())
    cfg = card['detectors']
    expected = {plane['name'] + '.h5' for key in ['field_planes', 'particle_planes'] for plane in cfg.get(key, [])}
    expected |= {plane['name'] + '-energy.h5' for plane in cfg.get('energy_planes', [])}
    if cfg.get('substep_particle_planes'):
        expected |= {plane['name'] + '-substeps.h5' for plane in cfg.get('particle_planes', [])}
    require({path.name for path in (directory / 'detectors').glob('*.h5')} == expected, 'Missing/extra prefix detectors')
    files = sorted((directory / 'detectors').glob('*.h5')) + sorted((directory / 'trajectory').glob('*.h5'))
    files.append(directory / 'ledger/energy-ledger.h5')
    scientific = {}
    counts = {}
    require(case.get('require_all_particles_through_planes') is False and case.get('input_source_ids_contiguous') is True,
            'Prefix reader requires explicit partial-cohort/source-ID declarations')
    first, last = case['expected_source_id_first'], case['expected_source_id_last']
    require(last - first + 1 == case['expected_particles'], 'Source range differs from input count')
    for path in files:
        with h5py.File(path, 'r') as handle:
            group = handle[next(iter(handle))]
            require(int(group['complete'][()]) == 0, 'Prefix HDF5 physical completion must remain0')
            for counter, values in [('committed_records', 'records'), ('committed_samples', 'time_s')]:
                if counter in group:
                    require(int(group[counter][()]) == len(group[values]), 'Uncommitted prefix extent: ' + str(path))
            def visit(name, dataset):
                if not isinstance(dataset, h5py.Dataset):
                    return
                if not dataset.shape:
                    finite(np.asarray(dataset[()]), path)
                    return
                block = 4 if name == 'fields' else 131072
                for start in range(0, len(dataset), block):
                    finite(dataset[start:start + block], path)
            group.visititems(visit)
            if 'fields' in group:
                require(len(group['fields']) == len(group['time_s']), 'Field prefix extents disagree')
                require(np.all(np.diff(group['time_s'][:]) > 0), 'Nonmonotonic prefix field time')
            if 'records' in group and 'source_id' in (group['records'].dtype.names or ()):
                records = group['records'][:]
                require(len(records) <= case['expected_particles'], 'Partial prefix exceeds the full input')
                require(np.all((records['source_id'] >= first) & (records['source_id'] <= last)), 'Prefix contains unknown source IDs')
                require(len(np.unique(records['source_id'])) == len(records) and
                        len(np.unique(records['particle_id'])) == len(records), 'Repeated prefix particle identity')
                require(np.all(records['weight'] == 1.), 'Physical-electron prefix weight changed')
                counts[path.name] = len(records)
        scientific[str(path.relative_to(directory))] = True
    rows, sums, ranks, steps = load_prefix_ledger(directory, case)
    cadence = int(card['energy_ledger']['sample_interval_steps'])
    expected_steps = list(range(0, steps + 1, cadence))
    if expected_steps[-1] != steps:
        expected_steps.append(steps)
    require(np.array_equal(rows['step'], expected_steps), 'Prefix ledger sampling/final flush differs')
    _, _, ledger = report_discrete.calculate(rows, sums, steps)
    histories = sorted((directory / 'reference').glob('source-coupling-rank-*.csv'))
    require([path.name for path in histories] ==
            [f'source-coupling-rank-{rank:05d}.csv' for rank in range(ranks)], 'Missing prefix coupling rank history')
    force = work = 0.
    for path in histories:
        # This marker denotes closed coupling records, not physical completion.
        require(path.with_suffix('.complete').read_text().split() == ['1', str(steps)], 'Unclosed coupling prefix')
        data = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        require(np.array_equal(data['step'], np.arange(1, steps + 1)), 'Incomplete coupling prefix steps')
        finite(data, path)
        force = max(force, float(np.max(data['self_force_defect'])))
        work = max(work, float(np.max(data['relative_fixed_point_work_defect'])))
    require(force < 1e-13 and work < 1e-8, 'Prefix fixed-point iteration failed')
    events = (directory / 'aprl.log').read_text().count('[particle-balance] step=')
    require((events == 0) if case['mode'] == 'off' else (events > 0), 'Wrong static/balanced prefix event count')
    return dict(passed=True, program_integrity_pass=True, prefix_complete=True,
        physical_complete=False, physical_acceptance=False, steps=steps, ranks=ranks,
        duration_exit=status, scientific_files=scientific, strict_ledger=ledger,
        max_force_defect=force, max_work_defect=work, balance_events=events,
        particle_plane_counts=counts, complete_source_through_all_planes=False,
        verification_scope='closed150-step prefix, not full beam/FEL completion')
