Full-size incremental optimization comparison, 2026-10-26

This is one fresh, exclusive two-node Slurm allocation: 256 MPI ranks,
one thread per rank, no fixed node names. It measures the new combination
of Cowan row reuse, particle-boundary hit reuse and diagnostic-only sampling
against the previously accepted Cowan build. It does not re-count the earlier
Cowan interior-path gain and does not run the historical pre-optimization code.

Control: build/cowan-accepted-r2-20261025/aprl, compiled77f1ef9, SHAdc537da...
Trial: build/compute-followup-20261026/aprl, compiled8adce33, SHA5edd096...
Both exact full identities, real compile flags, toolchain selections, trial
72-required evidence and the accepted combination evidence are verified and
frozen. Both use the same two frozen, read-only FFTW libraries from accepted-r2.

Prepare with run.py prepare --root NEW_ABSOLUTE_RESULT_DIRECTORY. Preparation
is metadata-only: no simulator, Slurm, MPI initialization or particle loading.
The frozen reference is production18893's compact card and input identity.
Every full run retains all 41,732,775 particles, 160x160x46592 geometry, held
startup, original boost and resolution, 0.3m retirement layer, physical stop,
detectors and cadence. Only mesh.duration (600 complete steps) and output
paths change. The full input is shared read-only, and its SHA is rechecked
on the allocated nodes. Existing production18893 and cluster.env stay unchanged.

Submit run.sbatch NEW_ABSOLUTE_RESULT_DIRECTORY. The two-hour Slurm hard cap
is a maximum, not a runtime estimate. Each full run is separately capped at
50 minutes including initialization/close; each small preflight at10 minutes.
A consumed outer timeout is a failure, never a successful bounded simulation.
Actual allocated nodes are recorded, and no script submits a later job.

Sequence: frozen identity checks, per-node dependency resolution, 256x1 binding,
paired small held/retirement runs, strict small scientific-output comparison,
control600 then trial600, each checked before continuing, final strict scientific
output comparison and matching step-window timing. Profiling is disabled.

Duration exhaustion intentionally returns exit1 and marks output incomplete.
Acceptance still requires the exact exhaustion message, root detector/ledger/
final-report close evidence, root duration-readiness-barrier entry, all600 steps,
unchanged geometry, and finite readable committed output. A peer's incomplete
error prefix cut by MPI_Abort is recorded only when all other gates pass; every
fully formed rank-error message must be the expected duration error. A later
expected message cannot hide an earlier different error. Barrier completion
logging is not required because a peer can abort after the barrier returns.
The old18906 run remains rejected because close/final evidence and H5 are broken.

Timing intervals are50-100,100-400,400-600. There is only one pair, so this is
an early-production sample, not repeated production statistics or validation
of late retirement/far-field behavior. Existing cross-rank differences are not
reclassified or hidden by tolerance. All checked scientific data and attributes
must match bit-for-bit; fixed-lab time cadence retains the existing FP64 ULP/
accumulated-time bound, with its five offline positive/negative regressions.

Read slurm-JOBID.out for sequence progress, control/aprl.log and trial/aprl.log
for step timing, comparison.json/PASS only after full checking. The startup
confirmation and any decision to stop following the job belong to the parent
workflow; this script itself runs the bounded pair to completion or first failure.

Offline tests: python test_cadence.py; python test_duration_exit.py. The latter
reads the accepted200-step and old broken600-step logs and checks synthetic
unknown-error, truncated-prefix, timeout, missing-close and wrong-step cases.
