Checkpoint/restart correctness campaign (draft until first measured PASS)

run.py prepare writes only tiny test inputs/configuration. It never starts MPI
or submits Slurm. run.py run uses an externally supplied launcher argument array;
there is no shell interpolation. Use a real Slurm allocation, or the separately
authorized CPU-only gnode1 exception with --allow-gnode1. Do not run on admin.
The requested test layout is limited to 16 CPUs; a Slurm allocation itself may
remain an exclusive full node according to workspace policy. Binding validation
and the required suite must precede this campaign and be recorded by its caller.

Example preparation (replace each /absolute/... with an actual path):
python3 run.py prepare --root /absolute/new-result-directory \
  --binary /absolute/candidate/aprl --baseline /absolute/accepted/aprl \
  --ranks 4 --threads 1 --timeout 180 \
  --launcher-json '["srun","--mpi=pmix_v3","-n","{ranks}","-c","{threads}","--cpu-bind=cores"]'
Inside the allocated job: python3 run.py run --root /absolute/new-result-directory
For a pre-bound gnode1 run, prepare a suitable mpiexec launcher and add
--allow-gnode1 to run. No production scripts, binaries or results are modified.
Prepare a separate directory for every 1/2/4/8 MPI test. The campaign checks
same-rank continuation; it intentionally rejects restart with different MPI size.

Cases:
1. generated-gaussian CPML: 64 particles, 16x16x256 mesh, nonzero self field,
   snapshots at steps0,7,17,55 and requested boost times0,13.7fs. Ledger cadence5
   and trajectory/field periods differ from snapshot steps. This checks default
   disabled output against the previous accepted binary, then checks enabling
   snapshots leaves all continuous-run scientific data bitwise identical.
2. held+retirement: the same small mesh plus a reproducible 64-record HDF5-v4
   auxiliary input (held cannot use generated-gaussian by design). Checkpoints
   at1,7,17,55 cross held activation and current-carrier retirement. The campaign
   requires all64 real particles removed before55 and checks continuation after
   that point. This fixture is a correctness probe, not a new physical model or
   a change to any production window, layer or physical stopping condition.
3. ballistic-validation: existing small generated fixture, validation enabled;
   snapshots1,7,17 must retain reference-plane histories and output cadence.

Each run must reach its configured physical stop and return0. A duration-exhausted
exit1 is not silently accepted for these fixtures; enlarge a fixture's guard only
after inspecting why its intended physical stop was not reached. The unchanged
production behavior still permits expected duration1 in separately bounded tests.
The HDF5 checks require complete1, committed lengths, finite data, identical types,
scientific attributes and per-field bit patterns. Compound padding is excluded.
Restart files are compared to the exact scientific suffix beginning at committed
prefix counts read from every rank's metadata/diagnostic_prefixes snapshot block;
ledger step rows must additionally equal continuous rows with step>snapshot.
Restart parent/step/time/output-scope provenance is checked before being excluded
from scientific attribute equality. Source prefix records remain in parent files.

Negative gates use an incomplete archive header, actual different MPI count,
changed mesh spacing, and a post-retirement layer modification even with explicit
allow_retirement_change. A timeout or successful physical run cannot satisfy a
negative gate. Logs, cards and launch.json are retained per run; RESULT.json and
PASS are emitted only after all gates, FAIL.json identifies any first failure.
Negative steps/times, nonfinite time, zero timeout and empty restart path are also
rejected. A second-generation restart tests prefixes relative to an already
restarted parent. Every held restart runs with the isolated original input renamed
out of the way and restores it afterwards, proving no original input-file lookup.
Explicit gates require nonzero CPML history at17, frozen particles at1, nonzero
retirement carriers at17, no real particles/carriers at55, and unmatched ballistic
reference entries at7. These are fixture requirements to be checked in actual
results, not claims that the as-yet unexecuted fixture has passed.

Required implementation contract:
- filename: checkpoint.directory/step-000000000007.h5 (12-digit absolute step).
- HDF5 /header is uint64 array; completion is element2 (checkpoint archive v1).
- full state and diagnostics cadence restored without a duplicate initial sample.
- fresh restart output is the suffix, with restart provenance in manifest/files.
- [checkpoint] phase=write step=... bytes=... seconds=... gib_per_s=... path=...
  reports write pause and size. The script records these lines without interpreting
  small-file bandwidth as full-scale multi-node throughput.

Additional required acceptance beyond this first campaign:
- Positive reset_diagnostics and allowed pre-retirement layer branch, including
  changed diagnostic planes; preserve the parent files and label the new branch.
- Reject existing output collision, corrupt block shape and changed immutable
  physics beyond mesh spacing.
- Isolated unavailable-peer/write-failure watchdog tests; no system/MPI changes.
- Measure disabled-loop overhead with repeated profiling-off paired runs before
  any performance claim; then a bounded two-node moderate-state I/O measurement.
No scientific/radiation convergence or production speedup follows from this suite.
