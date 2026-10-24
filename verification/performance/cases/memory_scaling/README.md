# gnode1 workload, rank, and thread experiment

This candidate-only experiment investigates the accepted `ac58e67` executable.
It does not rebuild a binary, use the production ELEGANT input, launch an old
version, or modify the currently running production job.

The user authorized bounded CPU-only SSH tests on `gnode1`, a 16-physical-core,
32-logical-CPU host with two sockets and eight NUMA domains. Execution rejects
other hosts and requires `APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1`. It holds the same
workspace lock as the previous CPU calibration, checks initial load, and has
one two-hour deadline across pilot and matrix stages. Each process launch has
a 180-second timeout. This is not a scheduler reservation against other users.

## Cases and controls

All cards use a fixed-seed generated Gaussian, 100000 represented electrons,
sigma `(0.25,0.25,0.25)` micrometres, unchanged boost and cell spacing, self-field
feedback, initial Poisson accuracy, CPML, a fixed physical endpoint, and full
diagnostics. The reduced sigma fits the smallest physical aperture; the pilot
must confirm zero loss and full diagnostic cohorts before expansion.

- Grid series: `16²×128`, `32²×128`, `64²×128`, `128²×128`; 16384 particles.
- Particle series at `32²×128`: 1024, 4096, 16384, 65536, 262144 particles.
- There are eight unique cards; the shared central case is counted once.
- Compact layouts: `1×1`, `2×1`, `4×1`, `8×1`, `16×1`, `8×2`, `4×4`, `2×8`, `1×16`.
- An additional `8×1` layout maps one rank per NUMA domain, contrasting with
  compact placement. Binding probes record the actual CPU and NUMA assignment.

Each card is byte-identical across layouts and repetitions. Different cards
change work and sometimes boundary/discretization physics; their time differences
are not implementation speedups. The beam remains narrow in z: rank load balance
and serial deposition/ledger work are confounders and are recorded by profiling.
Fewer compact cores also use fewer NUMA domains and memory channels; this is not
an isolated pure-core-count experiment.

The pilot runs smallest and largest cards at `16×1`, `1×1`, and `1×16`. It first
checks binding for every requested layout. Full matrix execution requires
`PILOT_PASS`. The matrix performs one warmup per card/layout, three throughput
repetitions in reproducibly shuffled order, and separate profiling runs at
`16×1`, `4×4`, and `1×16`. Front anchoring makes the actual step count depend on
the generated head position. Each case's first completed `16×1` control supplies
its expected step count; all layouts/repeats must match it. No endpoint or duration
is extended for timing. Primary throughput uses steps 20–80, excluding the first 20 steps and final output closure; complete
loop and process wall times remain available. Profiling is disabled in throughput.

Checks reuse the established scientific inspector: normal physical stop,
complete/finite HDF5 products, unique source IDs, all generated particles at the
detector, no losses, valid ledger closure, fixed binary/mesh identity, and strict
bitwise output equality against a candidate single-thread control at the same
MPI rank count. Cross-rank differences are reported at the existing tolerance
without relaxing it or misrepresenting those comparisons as passed.

No cache-miss or hardware bandwidth counters are assumed available. RSS and
configured cells per rank are reported as proxies, not measured cache working
sets. One optional `/proc/<pid>/numa_maps` snapshot is collected only during a
profile trial; it includes runtime/shared pages and is not a remote-access count.
Typed profile summaries record rank-mean/rank-maximum stage times, the slowest
rank IDs, and min/max/mean particles per rank over the same 20–80-step interval.
Communication stages include waits and are nested in other stages, so they are
not additive network time. Particle preparation has 2048-particle batches and
uses a serial path below 32 particles; this can also explain small-case changes.

An independently prepared first-touch/array-traversal microbenchmark may provide
conditional NUMA or cache-capacity evidence. Its effective byte rate is neither
a hardware DRAM counter nor an APRL throughput/speedup measurement; keep those
results separate from this scientific simulation matrix.

## Preparation and execution

Prepare a new result directory on the admin host; this only writes metadata:

    python run.py --stage prepare --root /absolute/new/result/directory

The prepared directory contains frozen scripts, case cards, binary/probe hashes,
and commands. On the already authorized `gnode1` host, use the frozen launcher:

    APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 bash ROOT/scripts/run.sh --root ROOT --stage pilot
    APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 bash ROOT/scripts/run.sh --root ROOT --stage matrix

`--stage all` runs pilot and matrix sequentially with the same gates and deadline.
Never resume a failed campaign silently; preserve it and prepare a new reviewed
directory. `report.json`, per-run comparison files, `summary.csv`, and binding
records preserve the evidence. Only `throughput` rows enter timing summaries.
