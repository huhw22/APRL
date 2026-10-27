Cowan quasi-steady initializer: small full-program integration

run.py prepare creates immutable cards and freezes executable/input hashes. It
does not run APRL. Use a fresh directory for each of 1, 2, and 4 MPI ranks, in
that order, with one thread per rank. All computation is bounded to 180 seconds
per run and 2 hours per campaign on the explicitly authorized gnode1 CPU host.
The shared CPU lock and binding probe run before simulation. Production runs,
default binaries, MPI/network settings and large ELEGANT input remain untouched.

Basic cards use the built-in generated-gaussian input: 4096 macroparticles,
1 pC, seed 271027, 32 x 32 x 192 cells, dx=dy=2 um and dz=.5 um. There is a
static-in-boost PEC case and beta_box=.23 CPML case. All three variants of a
case preserve its particle seed, physical mesh, end position, diagnostics,
boundary and snapshot times:
  accepted: frozen old source, relativistic-poisson
  legacy:   new executable, relativistic-poisson
  cowan:    new executable, cowan-relativistic-poisson
The two binaries must be built with the same compiler/dependencies/build flags.
The runner freezes binaries but build-flag verification belongs to the build
record prepared by the caller.

--include-held adds an established 64-particle HDF5-v4 arrival-event fixture as
a supplementary integration test. This is necessary because generated-gaussian
is a common-time snapshot; the full program explicitly requires HDF5-v4 for
frozen_lab_plane_injection. This fixture is not the baseline performance case.
Its early snapshot contains unreleased particles, and all 64 must eventually
activate. Cross-rank activation pushes are recorded, not silently assumed.

Acceptance separates numerical software equivalence from physical approximation:
* Old executable versus new executable with the old model: all scientific H5
  datasets and field/particle/CPML/carrier snapshot blocks must match bitwise.
* Old/new initializer step-0 particles must match bitwise, preserving initial
  positions, momenta, IDs, weights and held flags.
* New initializer step-7 restart must match every post-snapshot diagnostic
  suffix bitwise. CPML history must become nonzero. Changing the initializer
  model on restart must explicitly reject the physics contract.
* REVIEW.json measures delta kinetic energy, mean lab energy shift, lab/boost
  energy spread, interior field energy, energy residual and snapshot RMS sizes.
  New-model physical output is not required to match old-model bitwise. A
  software PASS is not automatic physical acceptance of every parameter range.
* aggregate compares all six full-plane E/B arrays at identical sample times
  across 1/2/4 ranks. It reports both old and new cross-rank discrepancy scales,
  without relaxing tolerances or claiming prior cross-rank differences fixed.

Example (prepare locally; execution only in the authorized CPU session):
  "$PYTHON_BIN" run.py prepare --root /absolute/fresh/rank1 --ranks 1 --include-held
  "$PYTHON_BIN" /absolute/fresh/rank1/scripts/run.py run --root /absolute/fresh/rank1 --allow-gnode1
Set APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 and load the same cluster runtime before
execution. Repeat rank2/rank4 only after rank1 passes. No Slurm ID is invented.
  "$PYTHON_BIN" run.py aggregate --roots /absolute/rank1 /absolute/rank2 /absolute/rank4 --output /absolute/cross-rank.json

The small-case elapsed times include initialization, snapshots and diagnostics;
they are not a production throughput benchmark. These cases do not establish
pure radiation, radiation convergence or full-scale speedup.

Independent verification and regression:
  Configure APRL_BUILD_COWAN_INITIALIZER_CHECKS=ON, then build all targets.
  The optional aprl_cowan_initializer_checks target uses FFTW only as a small
  serial DST oracle. Production APRL does not acquire an FFTW dependency.
  The 1/2/4-rank oracle cases include the existing minimum of two z cells per
  rank, CPML guards and transverse padding. They check a residual-derived
  error bound, not bitwise equality between different elliptic solvers.
  On authorized idle gnode1, set APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 and run:
    bash run_checks.sh /absolute/build /absolute/fresh/check-results
  Use MPIEXEC_PREFLAGS="--map-by;slot:PE=2;--bind-to;core;--nooversubscribe"
  for this suite, whose existing hybrid test uses up to 4 ranks x 2 threads.
  Build-local verification/runtime-libs may provide FFTW when needed on the
  test host; record its hash. run_checks.sh saves failed CTest detail as well
  as successes. No MPI/network/system configuration changes are needed.
  The default required.initial-field.cowan test has no FFTW dependency and
  checks static-source invariance under the production Maxwell step.

After completion, analyze.py creates the particle comparison and held plots
and PHYSICS_REVIEW.json from rank1's original HDF5 outputs. It uses a fixed
particle cohort, masks incomplete held trajectory samples, and separately
reports ordinary ledger residuals and their validity flags. See --help for
the input/output arguments. The figures compare finite-domain initializers;
they do not identify a pure radiation field or establish a speedup.
