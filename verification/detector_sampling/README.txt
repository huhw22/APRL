MPI fixed-lab detector interpolation regression — 2026-10-27

Problem and scope
-----------------
The old field-plane sampler clamps Ez/Bx/By to a local z slab's outermost
half-cell. At internal MPI boundaries the true trilinear stencil also needs
the neighbouring slab. A fixed laboratory plane moves across these boundaries
in the boosted grid, producing periodic diagnostic errors amplified by the
Lorentz transform. Solver halo exchange was already present; this fix changes
only LabDetectorManager sampling, not Maxwell or particle advancement.

The new path uses the global Yee stencil and fetches one adjacent raw Ez/Bx/By
layer only for a due sample requiring it. Buffers are reused. At 160x160 this is
619528 bytes per communicating rank; no whole-volume copy is introduced. No
transfer or buffer allocation occurs when field planes are disabled, not due,
or have an entirely local stencil. No throughput acceleration is claimed.

All six field components and their complex spatial phase are retained. The
existing native leapfrog E/B times are unchanged and described by HDF5 metadata.
The additional spatial_sampling attribute identifies the corrected stencil.
The completed exchange adds no checkpoint state or snapshot schema change.

Tests
-----
aprl_detector_sampling_checks OUTPUT_DIRECTORY [--sweep-only]
  Construct known analytic Yee fields, sample the actual HDF5-writing detector,
  compare every pixel/component to an unsplit EBFieldGrid. This is an analytic
  diagnostic fixture, not a Maxwell solution. Global grids are fixed across
  layouts; 257 cells checks uneven slabs and 256 cells checks moving planes.
  CTest required.detector.sampling-{1,2,4} checks boundaries, positive/negative
  gamma=929.13 boosts, a large coordinate origin, simultaneous planes and cadence.
  The sweep contains 1016 samples, 6x5 pixels, a strong linear background and a
  1e-5 box-field tone at Fourier bin37 with spatially varying phase.
  Default mode still checks exact internal interfaces and physical endpoints.
  --sweep-only exists for old-code controls: their generic sampler can reject an
  exact interface because of local-origin cancellation before recording a trace.

report.py FIXED_DIRECTORY --control OLD_DIRECTORY
  Compare actual old executable output to the legacy local sampler, then check
  raw full-plane Fourier amplitudes, phases and all 30x30 Ex pixel-pair cross
  spectra. No temporal filtering, clipping, mean subtraction or window is used.
  The figure omits only the DC bin from its displayed frequency range.

check_dynamics.py prepare --root FRESH_DIRECTORY --old OLD_APRL --new NEW_APRL
  --launcher-json '["/absolute/mpi-launch.sh","-n","{ranks}"]'
check_dynamics.py run --root DIRECTORY [--allow-gnode1]
  Three same-4-rank runs: old continuous, fixed continuous, fixed restart.
  Fixed-seed generated Gaussian, 64 particles, 12x12x80 grid, Cowan, CPML and
  self field. Compare scientific snapshots, non-field output and restart tail.
  Exclude identity metadata and x87 padding only; no tolerance for state changes.

Cluster use
-----------
Build with the existing MPI/parallel-HDF5 toolchain in a separate directory.
Set MPIEXEC_EXECUTABLE to this directory's mpi-launch.sh for cluster CTest.
run_gnode.sh BUILD FRESH_RESULT requires authorized gnode1, its CPU-test flag,
the existing standalone binding probe and the accepted FFTW runtime-library
copies in BUILD/runtime-libs. It runs one bounded regression suite.

run.sbatch RESULT BUILD reserves two exclusive 128-CPU nodes for at most10min;
each actual test uses64 ranks,32/node,one thread. It checks binding, runs fixed
full tests and expects the old control to exit1 with a completed raw field file.
The control binary BUILD/aprl_detector_sampling_legacy_checks is the exact same
fixture object linked to the unchanged accepted pre-fix direct_eb_core archive.
The recorded link command and all hashes are in legacy-build-final.json.
Reuse the accepted runtime environment; no MPI/network/system tuning is tested.

Accepted evidence
-----------------
Workspace results/detector_sampling_20261027/RESULT_REVIEW.txt and PROVENANCE.json.
75/75 required tests passed; changed fixture also reran1/2/4 successfully.
Slurm18931: node7/node8,19s,batchCOMPLETED0. Its old-control step deliberately
exits1; this expected regression failure is not a failed acceptance job.
Corrected64-rank output matches all182880 reference scalars exactly. Weak tone
complex amplitude and full pixel-pair cross spectra also match. Old control
reproduces the comb and matches legacy interpolation exactly.
Gaussian dynamics:42 steps,43 saved frames,owner changes10/33,restart7; snapshots
at7/30 and particle/ledger outputs unchanged; restart tail matches continuously.

Initial failed attempts remain recorded: missing FFTW deployment libraries on
gnode1 were restored from the accepted runtime bundle; job18930 hit an old-
sampler exact-interface rounding exception in the fixture. The revised sweep
phase avoids that exception for both versions; no threshold was relaxed and
the fixed default still tests exact interfaces.

Limits
------
This establishes removal of the identified detector-boundary comb, not all
possible MPI-dependent numerical errors. Old18893 field history has no saved
neighbour-grid history and cannot be repaired exactly by deleting spectral
lines. A new run or a suitable earlier full-state snapshot is needed for new
corrected data. Old diagnostic prefixes retain their original sampling error.
No production18893 rerun or default-binary switch is part of this acceptance.
E/B temporal alignment, near-field subtraction and pure-radiation attribution
remain separate questions. Synthetic spectral amplitudes are not18893 energies.
