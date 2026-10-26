# Checkpoint-disabled performance check

This benchmark measures the cost of adding the snapshot implementation when
no snapshot is requested. It compares the previously accepted combined build
(`8adce330f646701123ef5b312b914b37152af874`, executable SHA256
`5edd096eec8d02fef03613a4f74731d90eba5574cb6c109fe26c7fb4f16a9483`)
with a separately built checkpoint candidate. The production default and the
older Cowan-only baseline are not used.

Every simulation uses gnode1, **4 MPI ranks × 4 OpenMP threads**, the same 16
physical cores, profiling disabled, and an unchanged generated-gaussian card.
The two cards are `grid128` (128³ cells, 16,384 particles) and
`particles262144` (32×32×128 cells, 262,144 particles), copied byte-for-byte
from the accepted memory-scaling campaign. Neither card contains a
`checkpoint` block. Particle weights, physical endpoint, field precision,
substeps, feedback, initial field settings, boundaries and diagnostics remain
identical within each pair. Different cards are different workloads, not a
source of implementation speedup.

The driver reuses frozen copies of the existing fixed-MPI execution, binding,
identity, scientific inspection and timing helpers. It does not edit those
templates. It overrides only the selected two cases, accepted baseline
identity, and the campaign sequence; no profiling phase is allowed.

Before preparation, the candidate must have a completed 72-required-test
suite and a completed `verification/checkpoint_restart` correctness campaign
for the exact same candidate and accepted executables. The script binds both
evidence sets to their recorded hashes. It also checks identical build
configuration, exact compiled revision, accepted FFTW runtime libraries, the
source card hashes and the binding probe. All prepared scripts, cards,
programs, caches and evidence files are frozen in the new result directory.

Preparation only inspects/copies metadata and does not launch MPI:

```bash
/public/home/wfollowm/micromamba/envs/physics/bin/python3 \
  /public/home/wfollowm/APRL_Workspace/source/aprl/verification/checkpoint_performance/run.py \
  --stage prepare \
  --root /public/home/wfollowm/APRL_Workspace/results/checkpoint_disabled/NEW_DIRECTORY \
  --trial /public/home/wfollowm/APRL_Workspace/build/CHECKPOINT_BUILD/aprl \
  --trial-revision FULL_COMPILED_40_CHARACTER_REVISION \
  --required-pass /public/home/wfollowm/APRL_Workspace/results/REQUIRED_SUITE/PASS \
  --checkpoint-pass /public/home/wfollowm/APRL_Workspace/results/CHECKPOINT_CORRECTNESS/PASS
```

After resource review, execute the frozen launcher only on the authorized
CPU-only gnode1 host:

```bash
APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 bash \
  /public/home/wfollowm/APRL_Workspace/results/checkpoint_disabled/NEW_DIRECTORY/scripts/run.sh \
  --root /public/home/wfollowm/APRL_Workspace/results/checkpoint_disabled/NEW_DIRECTORY \
  --stage all
```

The launcher holds the existing CPU-calibration lock, rejects Slurm variables,
and limits the campaign to two hours. Each process invocation is limited to
180 seconds. The initial resource check rejects other MPI/APRL work or less
than 95% aggregate CPU idleness; it does not wait, kill other work, or modify
system/MPI settings. The existing binding probe confirms distinct physical
cores and the expected NUMA placement before any simulation.

Sequence: four pilots (both executables on each case), four warmups, then
five paired repetitions per case (20 measured runs): **28 simulations total**.
All runs have profiling off. Case order is deterministically shuffled and
control/trial order alternates. Every trial and later control is compared
bitwise with its case's first completed control. Missing outputs, changed
threads/geometry, a different physical endpoint, timeout, nonzero exit, or
scientific mismatch fails the campaign. Unexpected checkpoint log entries or
a created checkpoint directory also fail.

The primary timing is `(elapsed(step 80) - elapsed(step 20)) / 60`, excluding
initialization, the first 20 steps and output close. `report.json`,
`summary.csv` and `comparison.csv` contain raw repeats, median/min/max timings,
paired ratios and memory peaks; per-run logs and comparisons remain available.
`PILOT_PASS` covers only the pilots. `PASS` means all requested runs and
scientific checks finished, **not proof of exactly zero overhead**. A repeatable
slowdown must be investigated rather than relabelled as success; small timing
differences must be read against the five-pair variability. This experiment
does not measure snapshot write bandwidth or the complete production job.

A failed or incomplete campaign cannot silently rerun duplicate names. Keep
its evidence and prepare a fresh result directory. This template has no
submission side effect; creating it or invoking `--stage prepare` does not
start a benchmark.
