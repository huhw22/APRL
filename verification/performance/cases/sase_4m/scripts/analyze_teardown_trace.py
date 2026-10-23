#!/usr/bin/env python3
"""Locate an APRL shutdown hang from flushed per-rank teardown markers.

Exit status is 0 when at least one valid marker was analyzed, 1 when the
input contains no valid teardown marker, and 2 for a command/input error.
With --require-clean-exit, an incomplete detailed trace instead returns 3.
With --require-safe-finalize-hang, anything except a detailed, all-rank trace
that completed the solve and entered MPI_Finalize before hanging returns 4.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence, TextIO


SOLVE_TAIL_OPERATIONS = (
    ("final-stop-sync", "final stop synchronization"),
    ("final-energy-sample", "final energy sample"),
    ("trajectory-close", "trajectory close"),
    ("detector-close", "detector close"),
    ("ledger-close", "energy-ledger close"),
    ("particle-boundary-report", "particle-boundary report"),
    ("final-resource-report", "final resource report"),
)
SOLVE_TAIL_PHASES = tuple(
    phase
    for operation, _label in SOLVE_TAIL_OPERATIONS
    for phase in (f"{operation}-enter", f"{operation}-complete")
) + ("solve-complete",)
LEGACY_TEARDOWN_PHASES = (
    "simulation-destructor-enter",
    "simulation-destructor-complete",
    "mpi-finalize-enter",
    "mpi-finalize-complete",
)
FINALIZE_READINESS_PHASES = (
    "finalize-readiness-sync-enter",
    "finalize-readiness-sync-complete",
)
TEARDOWN_PHASES = (
    "simulation-destructor-enter",
    "simulation-destructor-complete",
) + FINALIZE_READINESS_PHASES + (
    "mpi-finalize-enter",
    "mpi-finalize-complete",
)
FULL_PHASES = SOLVE_TAIL_PHASES + TEARDOWN_PHASES
PHASE_INDEX = {phase: index for index, phase in enumerate(FULL_PHASES)}
FIELD_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")


@dataclass(frozen=True)
class Event:
    rank: int
    host: str
    pid: int
    elapsed_s: float
    phase: str
    source: str
    line_number: int


@dataclass
class ParsedTrace:
    events: list[Event]
    malformed: list[str]


def parse_lines(lines: Iterable[str], source: str) -> ParsedTrace:
    events: list[Event] = []
    malformed: list[str] = []
    for line_number, line in enumerate(lines, 1):
        if "[teardown]" not in line:
            continue
        fields = dict(FIELD_RE.findall(line.split("[teardown]", 1)[1]))
        try:
            rank = int(fields["rank"])
            host = fields["host"]
            pid = int(fields["pid"])
            elapsed_s = float(fields["elapsed_s"])
            phase = fields["phase"]
            if (rank < 0 or pid <= 0 or not host or
                    not math.isfinite(elapsed_s) or elapsed_s < 0):
                raise ValueError("invalid numeric/host field")
            if phase not in PHASE_INDEX:
                raise ValueError("unknown phase")
        except (KeyError, ValueError) as error:
            malformed.append(f"{source}:{line_number}: {error}: {line.strip()}")
            continue
        events.append(Event(
            rank=rank,
            host=host,
            pid=pid,
            elapsed_s=elapsed_s,
            phase=phase,
            source=source,
            line_number=line_number,
        ))
    return ParsedTrace(events=events, malformed=malformed)


def compact_ranks(ranks: Iterable[int]) -> str:
    values = sorted(set(ranks))
    if not values:
        return "none"
    ranges: list[str] = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def classify_teardown(expected: set[int],
                      ranks_by_phase: dict[str, set[int]]) -> tuple[str, str]:
    missing_enter = expected - ranks_by_phase["simulation-destructor-enter"]
    missing_destructor = expected - ranks_by_phase["simulation-destructor-complete"]
    missing_readiness_enter = (
        expected - ranks_by_phase["finalize-readiness-sync-enter"]
    )
    missing_readiness_complete = (
        expected - ranks_by_phase["finalize-readiness-sync-complete"]
    )
    missing_finalize_enter = expected - ranks_by_phase["mpi-finalize-enter"]
    missing_finalize_complete = expected - ranks_by_phase["mpi-finalize-complete"]

    if missing_enter:
        return (
            "INCOMPLETE_BEFORE_TEARDOWN",
            "Not every expected rank reached the Simulation destructor marker. "
            "The run stopped earlier than teardown, APRL_TEARDOWN_TRACE=all was not "
            "active, or the collected log is incomplete; this trace cannot yet "
            "attribute the hang to destructor, MPI_Finalize, or Slurm cleanup.",
        )
    if missing_destructor:
        return (
            "SIMULATION_DESTRUCTOR",
            "Every expected rank entered Simulation destruction, but one or more did "
            "not complete it. The hang is in Simulation/member destruction, before "
            "MPI_Finalize.",
        )
    # Preserve useful diagnosis for logs from the older APRL build, which had
    # MPI_Finalize markers but predated the readiness barrier.  Such a trace is
    # diagnostic only and still fails both strict acceptance predicates below.
    legacy_finalize_sequence = (
        not ranks_by_phase["finalize-readiness-sync-enter"]
        and not ranks_by_phase["finalize-readiness-sync-complete"]
        and bool(ranks_by_phase["mpi-finalize-enter"])
    )
    if not legacy_finalize_sequence and missing_readiness_enter:
        return (
            "BETWEEN_DESTRUCTOR_AND_FINALIZE_READINESS",
            "Every expected rank completed Simulation destruction, but one or more "
            "did not enter the finalize-readiness barrier. The short "
            "MPI_Finalize timeout was not armed; suspect rank failure or an "
            "incomplete log in the gap before readiness synchronization.",
        )
    if not legacy_finalize_sequence and missing_readiness_complete:
        return (
            "FINALIZE_READINESS_SYNC",
            "Every expected rank entered the finalize-readiness barrier, but one or "
            "more did not complete it. The hang is in the readiness collective "
            "guarded by the general teardown timeout, before any rank arms the "
            "short MPI_Finalize timeout.",
        )
    if missing_finalize_enter:
        if legacy_finalize_sequence:
            return (
                "BETWEEN_DESTRUCTOR_AND_MPI_FINALIZE",
                "This legacy trace predates finalize readiness. Every expected rank "
                "completed Simulation destruction, but one or more did not emit "
                "mpi-finalize-enter; suspect process failure or an incomplete log "
                "in the gap before MPI_Finalize.",
            )
        return (
            "BETWEEN_FINALIZE_READINESS_AND_MPI_FINALIZE",
            "Every expected rank completed finalize readiness, but one or more did "
            "not emit mpi-finalize-enter. These markers are adjacent in main(), so "
            "suspect process failure or an incomplete/corrupted log rather than a "
            "normal blocking operation.",
        )
    if missing_finalize_complete:
        return (
            "MPI_FINALIZE_OR_PMIX",
            "Every expected rank entered MPI_Finalize, but one or more did not return. "
            "The hang is in the MPI_Finalize/OpenMPI/PMIx teardown path, not in the "
            "APRL Simulation destructor.",
        )
    return (
        "POST_MPI_FINALIZE",
        "Every expected rank completed Simulation destruction and returned from "
        "MPI_Finalize. Combine this trace with a live APRL PID check: remaining APRL "
        "processes indicate a C/C++ process-exit, static-destructor, or runtime hook "
        "hang after MPI_Finalize; no APRL processes while the step remains RUNNING "
        "indicates srun/slurmstepd step cleanup.",
    )


def classify(expected: set[int], ranks_by_phase: dict[str, set[int]],
             detailed_solve_tail: bool) -> tuple[str, str]:
    if detailed_solve_tail:
        previous_operation = "solve-body"
        for operation, label in SOLVE_TAIL_OPERATIONS:
            enter_phase = f"{operation}-enter"
            complete_phase = f"{operation}-complete"
            missing_enter = expected - ranks_by_phase[enter_phase]
            if missing_enter:
                code = "SOLVE_TAIL_BEFORE_" + operation.upper().replace("-", "_")
                return (
                    code,
                    "Not every expected rank reached %s. All earlier recorded stages "
                    "ended at %s; suspect the gap before this operation, an earlier "
                    "rank failure, or an incomplete log. Missing ranks: %s."
                    % (enter_phase, previous_operation, compact_ranks(missing_enter)),
                )
            missing_complete = expected - ranks_by_phase[complete_phase]
            if missing_complete:
                code = "SOLVE_TAIL_" + operation.upper().replace("-", "_")
                return (
                    code,
                    "Every expected rank entered %s, but one or more did not complete "
                    "it. The hang is inside the solve-tail %s operation. Missing "
                    "completion ranks: %s."
                    % (enter_phase, label, compact_ranks(missing_complete)),
                )
            previous_operation = complete_phase

        missing_solve_complete = expected - ranks_by_phase["solve-complete"]
        if missing_solve_complete:
            return (
                "SOLVE_TAIL_AFTER_FINAL_RESOURCE_REPORT",
                "All expected ranks completed the final resource report, but not all "
                "emitted solve-complete. The two markers have no intended collective "
                "between them; suspect rank failure or incomplete logging. Missing "
                "ranks: %s." % compact_ranks(missing_solve_complete),
            )

    return classify_teardown(expected, ranks_by_phase)


def phase_rank_sets(events: Sequence[Event]) -> dict[str, set[int]]:
    result = {phase: set() for phase in FULL_PHASES}
    for event in events:
        result[event.phase].add(event.rank)
    return result


def has_detailed_solve_tail(events: Sequence[Event]) -> bool:
    return any(event.phase in SOLVE_TAIL_PHASES for event in events)


def has_finalize_readiness(events: Sequence[Event]) -> bool:
    return any(event.phase in FINALIZE_READINESS_PHASES for event in events)


def trace_has_strict_integrity(events: Sequence[Event],
                               expected_count: int) -> bool:
    """Reject concatenated, duplicated, conflicting, or reordered traces."""
    expected = set(range(expected_count))
    if {event.rank for event in events} != expected:
        return False
    keys = [(event.rank, event.phase) for event in events]
    if len(keys) != len(set(keys)):
        return False

    identity_by_rank: dict[int, set[tuple[str, int]]] = defaultdict(set)
    ranks_by_identity: dict[tuple[str, int], set[int]] = defaultdict(set)
    event_by_rank_phase: dict[tuple[int, str], Event] = {}
    for event in events:
        identity = (event.host, event.pid)
        identity_by_rank[event.rank].add(identity)
        ranks_by_identity[identity].add(event.rank)
        event_by_rank_phase[(event.rank, event.phase)] = event
    if any(len(identities) != 1 for identities in identity_by_rank.values()):
        return False
    if any(len(ranks) != 1 for ranks in ranks_by_identity.values()):
        return False

    for rank in expected:
        phase_times = [
            event_by_rank_phase[(rank, phase)].elapsed_s
            for phase in FULL_PHASES
            if (rank, phase) in event_by_rank_phase
        ]
        if any(later < earlier
               for earlier, later in zip(phase_times, phase_times[1:])):
            return False
    return True


def full_trace_is_clean(events: Sequence[Event], expected_count: int,
                        malformed: Sequence[str] = ()) -> bool:
    expected = set(range(expected_count))
    by_phase = phase_rank_sets(events)
    return (
        not malformed
        and trace_has_strict_integrity(events, expected_count)
        and has_detailed_solve_tail(events)
        and all(expected <= by_phase[phase] for phase in FULL_PHASES)
    )


def trace_is_safe_finalize_hang(events: Sequence[Event],
                                expected_count: int,
                                malformed: Sequence[str] = ()) -> bool:
    """Return true only for the campaign's narrowly accepted teardown fault.

    Every rank must have completed every solve-tail/HDF5-close phase,
    destroyed Simulation, completed the finalize-readiness barrier, and
    entered MPI_Finalize. At least one rank must not have returned from
    MPI_Finalize; zero successful returns are acceptable because the readiness
    barrier proves every rank was safe to terminate first. Requiring the detailed solve-tail
    format prevents a legacy four-marker trace from being accepted without
    evidence that the scientific outputs were closed first.
    """
    expected = set(range(expected_count))
    by_phase = phase_rank_sets(events)
    phases_before_finalize_return = FULL_PHASES[:-1]
    finalize_returned = expected & by_phase["mpi-finalize-complete"]
    return (
        not malformed
        and trace_has_strict_integrity(events, expected_count)
        and has_detailed_solve_tail(events)
        and all(expected <= by_phase[phase]
                for phase in phases_before_finalize_return)
        and len(finalize_returned) < expected_count
        and classify(expected, by_phase, True)[0] == "MPI_FINALIZE_OR_PMIX"
    )


def render_report(events: Sequence[Event], expected_count: int | None,
                  malformed: Sequence[str], streams: Sequence[str]) -> str:
    observed = {event.rank for event in events}
    if expected_count is None:
        expected = set(observed)
        expected_label = "inferred from observed markers"
    else:
        expected = set(range(expected_count))
        expected_label = "explicit ranks 0..%d" % (expected_count - 1)
    unexpected = observed - expected
    if unexpected:
        raise ValueError(
            "observed ranks outside --expected-ranks: %s" % compact_ranks(unexpected)
        )

    detailed_solve_tail = has_detailed_solve_tail(events)
    readiness_trace = has_finalize_readiness(events)
    report_phases = FULL_PHASES if detailed_solve_tail else TEARDOWN_PHASES
    by_phase = phase_rank_sets(events)
    unique_keys = {(event.rank, event.phase) for event in events}
    duplicates = len(events) - len(unique_keys)
    host_by_rank: dict[int, set[str]] = defaultdict(set)
    for event in events:
        host_by_rank[event.rank].add(event.host)

    lines = [
        "APRL teardown trace analysis",
        "sources: %s" % ", ".join(streams),
        "valid marker records: %d (unique rank/phase=%d, duplicates=%d)"
        % (len(events), len(unique_keys), duplicates),
        "observed ranks: %d [%s]" % (len(observed), compact_ranks(observed)),
        "expected ranks: %d [%s] (%s)"
        % (len(expected), compact_ranks(expected), expected_label),
        "trace format: %s" % (
            "detailed solve-tail + teardown" if detailed_solve_tail
            else (
                "readiness-aware six-stage teardown" if readiness_trace
                else "legacy four-stage teardown"
            )
        ),
    ]

    conflicting_hosts = {
        rank: hosts for rank, hosts in host_by_rank.items() if len(hosts) != 1
    }
    host_counts = Counter(
        next(iter(hosts)) for hosts in host_by_rank.values() if len(hosts) == 1
    )
    lines.append("rank placement by host: " + (
        ", ".join(f"{host}={count}" for host, count in sorted(host_counts.items()))
        if host_counts else "none"
    ))
    if conflicting_hosts:
        lines.append(
            "WARNING: rank host changed: "
            + "; ".join(
                f"rank {rank}={','.join(sorted(hosts))}"
                for rank, hosts in sorted(conflicting_hosts.items())
            )
        )

    lines.append("")
    lines.append("Phase coverage (counts are unique ranks):")
    for phase in report_phases:
        phase_events = [event for event in events if event.phase == phase]
        phase_ranks = by_phase[phase]
        missing = expected - phase_ranks
        per_host = Counter()
        for rank in phase_ranks:
            phase_hosts = {event.host for event in phase_events if event.rank == rank}
            for host in phase_hosts:
                per_host[host] += 1
        missing_by_host = Counter()
        missing_host_unknown = 0
        for rank in missing:
            known_hosts = host_by_rank.get(rank, set())
            if len(known_hosts) == 1:
                missing_by_host[next(iter(known_hosts))] += 1
            else:
                missing_host_unknown += 1
        missing_host_parts = [
            f"{host}={count}" for host, count in sorted(missing_by_host.items())
        ]
        if missing_host_unknown:
            missing_host_parts.append(f"unknown={missing_host_unknown}")
        elapsed = [event.elapsed_s for event in phase_events]
        lines.extend([
            f"  {phase}",
            "    reached: %d/%d" % (len(phase_ranks), len(expected)),
            "    missing ranks: %s" % compact_ranks(missing),
            "    hosts: %s" % (
                ", ".join(f"{host}={count}" for host, count in sorted(per_host.items()))
                if per_host else "none"
            ),
            "    missing by known host: %s" % (
                ", ".join(missing_host_parts) if missing_host_parts else "none"
            ),
            "    elapsed_s range: %s" % (
                "%.9g .. %.9g" % (min(elapsed), max(elapsed)) if elapsed else "none"
            ),
        ])

    ordering_violations: list[str] = []
    for rank in sorted(observed):
        first_by_phase: dict[str, Event] = {}
        for event in events:
            if event.rank == rank and event.phase not in first_by_phase:
                first_by_phase[event.phase] = event
        sequence = [first_by_phase[phase].elapsed_s for phase in report_phases
                    if phase in first_by_phase]
        phase_sequence = [phase for phase in report_phases if phase in first_by_phase]
        if any(later < earlier for earlier, later in zip(sequence, sequence[1:])):
            ordering_violations.append(
                f"rank {rank}: " + ", ".join(
                    f"{phase}={first_by_phase[phase].elapsed_s:.9g}"
                    for phase in phase_sequence
                )
            )
    if ordering_violations:
        lines.append("")
        lines.append("WARNING: non-monotonic first-marker timestamps:")
        lines.extend("  " + item for item in ordering_violations)

    diagnosis, explanation = classify(expected, by_phase, detailed_solve_tail)
    lines.extend(["", "Diagnosis: " + diagnosis, explanation])
    if expected_count is None:
        lines.append(
            "NOTE: --expected-ranks was omitted, so completely silent ranks cannot "
            "be detected. Use --expected-ranks 256 for this two-node campaign."
        )
    if malformed:
        lines.append("")
        lines.append("Malformed teardown-like lines: %d" % len(malformed))
        lines.extend("  " + item for item in malformed[:10])
        if len(malformed) > 10:
            lines.append("  ... %d more" % (len(malformed) - 10))
    return "\n".join(lines)


def read_source(name: str) -> tuple[TextIO, bool]:
    if name == "-":
        return sys.stdin, False
    return Path(name).open("r", encoding="utf-8", errors="replace"), True


def self_test() -> None:
    def event_lines(phases: Sequence[str], completed: dict[int, int]) -> list[str]:
        lines = []
        for rank, phase_count in completed.items():
            for index, phase in enumerate(phases[:phase_count]):
                lines.append(
                    "[teardown] rank=%d host=node%d pid=%d elapsed_s=%.3f phase=%s\n"
                    % (rank, rank // 2, 1000 + rank, index + rank / 10.0, phase)
                )
        return lines

    expected = set(range(4))
    cases = (
        ({0: 2, 1: 2, 2: 1, 3: 2}, "SIMULATION_DESTRUCTOR"),
        ({0: 4, 1: 4, 2: 3, 3: 4}, "FINALIZE_READINESS_SYNC"),
        ({0: 6, 1: 6, 2: 5, 3: 6}, "MPI_FINALIZE_OR_PMIX"),
        ({0: 6, 1: 6, 2: 6, 3: 6}, "POST_MPI_FINALIZE"),
    )
    for completion, wanted in cases:
        parsed = parse_lines(event_lines(TEARDOWN_PHASES, completion), "self-test")
        actual, _ = classify(expected, phase_rank_sets(parsed.events), False)
        assert actual == wanted, (actual, wanted)

    for operation, _label in SOLVE_TAIL_OPERATIONS:
        complete_phase = f"{operation}-complete"
        rank_two_count = FULL_PHASES.index(complete_phase)
        completion = {rank: len(FULL_PHASES) for rank in expected}
        completion[2] = rank_two_count
        parsed = parse_lines(event_lines(FULL_PHASES, completion), "self-test")
        actual, _ = classify(expected, phase_rank_sets(parsed.events), True)
        wanted = "SOLVE_TAIL_" + operation.upper().replace("-", "_")
        assert actual == wanted, (actual, wanted)

    complete_counts = {rank: len(FULL_PHASES) for rank in expected}
    complete_trace = parse_lines(
        event_lines(FULL_PHASES, complete_counts), "self-test"
    )
    assert full_trace_is_clean(complete_trace.events, len(expected))
    assert not trace_is_safe_finalize_hang(complete_trace.events, len(expected))
    assert classify(expected, phase_rank_sets(complete_trace.events), True)[0] == (
        "POST_MPI_FINALIZE"
    )
    legacy_trace = parse_lines(
        event_lines(
            LEGACY_TEARDOWN_PHASES,
            {rank: len(LEGACY_TEARDOWN_PHASES) for rank in expected},
        ),
        "self-test",
    )
    assert not full_trace_is_clean(legacy_trace.events, len(expected))
    assert not trace_is_safe_finalize_hang(legacy_trace.events, len(expected))
    assert classify(expected, phase_rank_sets(legacy_trace.events), False)[0] == (
        "POST_MPI_FINALIZE"
    )

    finalize_hang_counts = {rank: len(FULL_PHASES) for rank in expected}
    finalize_hang_counts[2] -= 1
    finalize_hang_trace = parse_lines(
        event_lines(FULL_PHASES, finalize_hang_counts), "self-test"
    )
    assert trace_is_safe_finalize_hang(
        finalize_hang_trace.events, len(expected)
    )
    assert not trace_is_safe_finalize_hang(
        finalize_hang_trace.events + finalize_hang_trace.events[:1],
        len(expected),
    )
    assert not trace_is_safe_finalize_hang(
        finalize_hang_trace.events, len(expected), ["malformed marker"]
    )

    conflicting_identity = list(finalize_hang_trace.events)
    conflicting_identity[-2] = replace(
        conflicting_identity[-2], pid=conflicting_identity[-2].pid + 10000
    )
    assert not trace_is_safe_finalize_hang(conflicting_identity, len(expected))

    reversed_time = list(finalize_hang_trace.events)
    finalize_enter_index = next(
        index for index, event in enumerate(reversed_time)
        if event.rank == 0 and event.phase == "mpi-finalize-enter"
    )
    reversed_time[finalize_enter_index] = replace(
        reversed_time[finalize_enter_index], elapsed_s=0.0
    )
    assert not trace_is_safe_finalize_hang(reversed_time, len(expected))

    no_finalize_return_counts = {
        rank: len(FULL_PHASES) - 1 for rank in expected
    }
    no_finalize_return_trace = parse_lines(
        event_lines(FULL_PHASES, no_finalize_return_counts), "self-test"
    )
    assert trace_is_safe_finalize_hang(
        no_finalize_return_trace.events, len(expected)
    )

    readiness_hang_counts = {rank: len(FULL_PHASES) for rank in expected}
    readiness_hang_counts[2] = FULL_PHASES.index(
        "finalize-readiness-sync-complete"
    )
    readiness_hang_trace = parse_lines(
        event_lines(FULL_PHASES, readiness_hang_counts), "self-test"
    )
    assert classify(
        expected, phase_rank_sets(readiness_hang_trace.events), True
    )[0] == "FINALIZE_READINESS_SYNC"
    assert not trace_is_safe_finalize_hang(
        readiness_hang_trace.events, len(expected)
    )

    early_hang_counts = {rank: len(FULL_PHASES) for rank in expected}
    early_hang_counts[2] = FULL_PHASES.index("detector-close-complete")
    early_hang_trace = parse_lines(
        event_lines(FULL_PHASES, early_hang_counts), "self-test"
    )
    assert not trace_is_safe_finalize_hang(early_hang_trace.events, len(expected))
    assert compact_ranks([0, 1, 2, 4, 6, 7]) == "0-2,4,6-7"
    malformed = parse_lines(["[teardown] rank=0 phase=bad\n"], "self-test")
    assert not malformed.events and len(malformed.malformed) == 1
    print("analyze_teardown_trace self-test passed")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Parse flushed APRL per-rank solve-tail/teardown markers and locate a "
            "final solve operation, Simulation destructor, MPI_Finalize/PMIx, or "
            "post-MPI shutdown hang."
        )
    )
    parser.add_argument(
        "logs", nargs="*", metavar="LOG",
        help="APRL/sbatch log to read; use - for standard input",
    )
    parser.add_argument(
        "--expected-ranks", type=int,
        help="expected MPI size; ranks are assumed to be 0..N-1",
    )
    parser.add_argument(
        "--require-clean-exit", action="store_true",
        help=(
            "return 3 unless every expected rank has every detailed solve-tail and "
            "teardown marker through mpi-finalize-complete"
        ),
    )
    parser.add_argument(
        "--require-safe-finalize-hang", action="store_true",
        help=(
            "return 4 unless every expected rank completed the detailed solve-tail "
            "and finalize-readiness barrier, then entered MPI_Finalize, with at "
            "least one rank not returning from MPI_Finalize; zero returns are "
            "allowed after readiness is proven"
        ),
    )
    parser.add_argument(
        "--self-test", action="store_true", help="run parser/classifier self-tests",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if not args.logs:
        parser.error("at least one LOG is required unless --self-test is used")
    if args.expected_ranks is not None and args.expected_ranks <= 0:
        parser.error("--expected-ranks must be a positive integer")
    if args.require_clean_exit and args.expected_ranks is None:
        parser.error("--require-clean-exit requires --expected-ranks")
    if args.require_safe_finalize_hang and args.expected_ranks is None:
        parser.error("--require-safe-finalize-hang requires --expected-ranks")
    if args.require_clean_exit and args.require_safe_finalize_hang:
        parser.error(
            "--require-clean-exit and --require-safe-finalize-hang are mutually "
            "exclusive"
        )
    if args.logs.count("-") > 1:
        parser.error("standard input may be specified only once")

    events: list[Event] = []
    malformed: list[str] = []
    for name in args.logs:
        try:
            stream, close_after = read_source(name)
            try:
                parsed = parse_lines(stream, name)
            finally:
                if close_after:
                    stream.close()
        except OSError as error:
            print(f"ERROR: cannot read {name}: {error}", file=sys.stderr)
            return 2
        events.extend(parsed.events)
        malformed.extend(parsed.malformed)

    if not events:
        print(
            "ERROR: no valid [teardown] markers were found. Rebuild the diagnostic "
            "APRL and submit with APRL_TEARDOWN_TRACE=all.",
            file=sys.stderr,
        )
        for item in malformed[:10]:
            print("  " + item, file=sys.stderr)
        return 1
    try:
        print(render_report(events, args.expected_ranks, malformed, args.logs))
    except ValueError as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 2
    if (args.require_clean_exit and
            not full_trace_is_clean(events, args.expected_ranks, malformed)):
        print(
            "ERROR: --require-clean-exit failed: not every expected rank has the "
            "complete detailed solve-tail and teardown marker sequence through "
            "mpi-finalize-complete.",
            file=sys.stderr,
        )
        return 3
    if (args.require_safe_finalize_hang and
            not trace_is_safe_finalize_hang(
                events, args.expected_ranks, malformed)):
        print(
            "ERROR: --require-safe-finalize-hang failed: the trace is not the "
            "narrowly accepted all-rank MPI_Finalize/PMIx hang after every "
            "detailed solve-tail/output-close phase and the finalize-readiness "
            "barrier completed.",
            file=sys.stderr,
        )
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
