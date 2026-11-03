import re,math,sys
from common import need as require
EXHAUSTED="mesh.duration was exhausted before the configured physical stop"
STEPS=120

def check_duration_exit(log, code, steps=STEPS):
    """Accept only the expected bounded stop after root closed all output.

    MPI_Abort may truncate a peer's error prefix after the readiness barrier.
    Keep incomplete prefixes as evidence, but never discard a fully formed
    rank error or let a following expected error conceal its message.
    """
    require(code == 1 and EXHAUSTED in log,
            'Expected clean duration-exhaustion exit code 1')
    expected_ranks = []
    truncated = []
    for position in re.finditer(r'Simulation failure', log):
        tail = log[position.start():]
        segment = re.split(r'\n|(?=Simulation failure)|(?=srun:)', tail[len('Simulation failure'):],
                           maxsplit=1)[0]
        segment = 'Simulation failure' + segment
        full = re.match(r'^Simulation failure on MPI rank ([0-9]+):[ \t]*', segment)
        if full:
            message = segment[full.end():]
            require(message.startswith(EXHAUSTED),
                    'Unexpected complete rank failure: ' + segment)
            expected_ranks.append(int(full.group(1)))
        else:
            require(re.fullmatch(r'Simulation failure(?: on MPI rank(?: [0-9]*)?)?[ \t]*', segment),
                    'Unexpected malformed simulation failure: ' + segment)
            truncated.append(segment)
    require(expected_ranks, 'Missing a complete expected duration failure message')
    for phase in ('detector-close-complete', 'ledger-close-complete',
                  'final-resource-report-complete', 'duration-exhaustion-readiness-sync-enter'):
        require(re.search(r'\[teardown\] rank=0\b[^\n]* phase=' + phase + r'\b', log),
                'Missing root output-close/readiness marker: ' + phase)
    require(re.search(r'\[resource\] phase=final step=' + str(steps) + r'\b', log),
            'Missing final complete-step report')
    return {'expected_complete_rank_messages': len(expected_ranks),
            'truncated_rank_prefixes': truncated,
            'root_outputs_closed_before_readiness_barrier': True,
            'readiness_complete_marker_observed': bool(re.search(
                r'\[teardown\] rank=0\b[^\n]* phase=duration-exhaustion-readiness-sync-complete\b', log))}

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
