Moving-window full-program integration fixture (version 2, 2026-10-28)

Purpose
-------
This is a bounded functional comparison of the existing fixed full-window held
path and the new moving-window path. The production Cowan advance, deposition,
particle push, original boost/dz and transverse CPML parameters are preserved.
The fixture is not a production radiation calculation or an N-convergence test.
No script submits jobs, changes default programs, or edits production input.

The normal generated-gaussian mode cannot exercise formal held injection:
Simulation requires fixed-plane arrival events in HDF5 format 4. prepare creates
a small deterministic v4 file for this specific functional purpose. Its sampled
transverse/momentum moments follow the previously frozen production-moment audit.
The nominal longitudinal RMS is 8 boosted dz. Two explicit +/-4.5 sigma arrival
records make the source span about 72 cells, exercising entry from behind the
64-cell lower guard. These are boundary witnesses, not random-tail statistics.
Default N=64 is an interface fixture, with approximately 495 electrons per macro.
Optional --particles 31675 retains the same charge and has near-unit weights;
it is a separate campaign and requires independent review. Neither choice is a
faithful reproduction of the full ELEGANT distribution.

Geometry and physics
--------------------
Both variants use 160 x 160 x 896 reference cells, dx=dy=5 um,
dz=1.2543262415359034 um, boost gamma=929.1305492858543, and center=0.
The Cowan initialization tolerance remains 1e-8. Both use xy CPML 10/10 and
z CPML 0 so the boundary setup is not silently changed between the pair.
The fixed case invokes the existing CG Cowan initializer. The moving case uses
the blocked reference on the SAME finite source domain and a compact active
window with lower/upper guards 64/128, block_layers=32, and automatic cells.
The source domain is not silently enlarged if the query range is inadequate.

Held release is at lab z=-0.05 m. The undulator retains K=1.092714657005,
period=0.0234 m, polarization=pi/2, Gaussian fringe cutoff=1e-9, and the original
80 particle-steps-per-period target; only its length is reduced to one period.
Particle planes are at release, an upstream auxiliary plane 16*gamma*dz after
release, and the downstream field plane at +0.06 m. The auxiliary plane is used
for a full-cohort ballistic comparison with an explicit fringe-clearance bound.
The exit plane includes magnetic dynamics; its width change is not attributed
solely to space charge. Raw E/B retain their existing native B half-time offset.

The runtime derives an active window of about 289 cells, initial reference
offset 422, and right particle safe-exit lab z about +0.206364 m. The stop
reference position is the safe exit plus 6*gamma*dz, beyond the field outer edge,
so particle exit and subsequent current-carrier cropping are exercised.
prepare reconstructs the v4 held-frame time origin and estimates the initial
reference-centre boost worldline before deriving the stop step. It also applies
initializeParticles' explicit-grid centring translation: approximately 38.02694
dz. The release plane is at this translated box coordinate, not box z=0.
Runtime field-cell count, initial offset and safe exit must match the prepared
version2 contract; the analyzer rejects mismatches. Version1 frozen result
directories are never rewritten and their logged actual geometry takes priority.
The reference-centre worldline
does not start at initial_center_z on the t'=0 plane: for N=64 it starts near
-0.0870 m. For N=64 the corrected stop takes 258 steps; duration is 258.5 steps,
allowing at most 259 advances and runtime query bound ceil(duration/dt)+1=260.
For N=31675 the predicted stop takes 257 steps. Both retain source ID2 initially
10.02694 cells behind the inlet, with physical drift crossing the old lower
boundary during step9, before the window shifts. These witness positions and
crossing estimates are saved in plan.json and checked before preparation.
The conservative beta=.25 source lookup bound leaves at least 96 reference cells
on the left. Runtime validates its measured beta and actual range again.
plan.json records all estimates; this is explicit trajectory/query geometry,
not a generic femtosecond quiet-time padding rule.

Build and prepare
-----------------
Build the current core with APRL_ENABLE_MOVING_WINDOW=ON and the established
MPI/HDF5 toolchain. Existing FFTW cache variables APRL_VERIFY_FFTW_INCLUDE_DIR
and APRL_VERIFY_FFTW_LIBRARY are reused. This feature remains OFF by default.
The ON verification targets include aprl_moving_detector_checks,
aprl_moving_window_mpi_checks, and aprl_blocked_reference_reader_checks.
Run module/required checks before the full-program campaign. Preparation only
writes fixture files, cards, frozen hashes and source provenance; it never runs
APRL. A fresh absolute directory and an existing executable are required.

  python run.py prepare --root /absolute/new/rank1 --binary /absolute/build/aprl
  python run.py prepare --root /absolute/new/rank2 --binary /absolute/build/aprl --ranks 2
  python run.py prepare --root /absolute/new/two-node --binary /absolute/build/aprl --ranks 2 --nodes 2

For an exact layout replay of an existing frozen campaign, append
  --geometry-from /absolute/previous/plan.json
This verifies every input particle record, source normalization/seed, units and
all non-path card settings, then copies that plan's source geometry, duration
and physical stop unchanged. The new plan records its provenance and hashes.
This option is required when comparing against a version1 run: corrected future
geometry may otherwise change the stop by a step. It never rewrites the old
plan, and its original geometry estimates remain clearly labelled as estimates.
New dense tests without this option use the corrected version2 geometry.

Default steps use 16 OpenMP threads per rank (at most 32 CPUs here), inside the
user's exclusive 128-CPU-per-node reservation. The actual rank/thread and node
layout is saved separately from the reservation. Use the matching run.sbatch
or run-two-nodes.sbatch with the prepared directory argument. Nodes are not
pinned. The runner requires a real Slurm allocation and a successful binding
probe; no unrestricted local or fabricated Slurm execution is supported.
Every srun step explicitly receives --chdir pointing to its resolved case path.
Layouts should be executed sequentially: 1 rank, 2 ranks on one node, then
2 ranks with one rank on each of two nodes. Preparation permits any of these;
the reviewer controls the validation order and checks available resources.

Each case has a MAXIMUM of 1200 seconds. The entire job still has a 1500-second
budget. Later cases receive min(1200, remaining budget minus 10 seconds), and
do not start with less than 60 seconds remaining. Existing 600/60-second APRL
shutdown watchdogs remain enabled; the external campaign timeout is also
bounded. Failure stops the campaign; incomplete output is never accepted.
Default order is moving, then fixed, catching fast moving-path failures first.
Optional --baseline adds old-fixed last for a same-card disabled-path bitwise
replay. All three share the SAME 1500-second budget; two slow CG solves may not
fit. --variants moving,fixed or another explicit comma-separated subset/order
freezes just those cases. Use a new directory and the remaining subset after
a budget-limited campaign; never overwrite completed output or rerun it silently.
The earlier binary must use the same compiler/dependencies and support Cowan initial
fields; the runner hashes it but cannot establish build-flag equivalence.

Review
------
  python analyze.py --root /absolute/completed/rank1 --output /absolute/new/review1
  python analyze.py --root /absolute/completed/rank1 --other-root /absolute/completed/rank2 --output /absolute/new/cross-layout

If cases were completed in separate frozen continuations, --fixed-case,
--moving-case and --old-fixed-case select their directories explicitly. Input
arrays and same-layout rank counts are checked before combining them.

The analyzer requires complete, committed, finite HDF5 outputs, a clean physical
stop, unique full source-ID cohorts at every particle plane, unchanged identity
and weights, and actual right particle exits/current-carrier crops. It reports:
* First recorded per-ID trajectory positions/momenta and input-to-release
  consistency. Formal trajectory output omits held particles, so these first
  released samples are explicitly NOT labelled as the t_box=0 initial state.
* Per-ID stable kinetic-energy differences, weighted mean/RMS, energy spread,
  and transverse widths/covariances. Small energy changes are computed by
  rationalizing delta gamma, never normalized by the much larger total energy.
* Upstream auxiliary widths/covariances after subtracting an analytic ballistic
  drift from each particle's own release position/momentum. Every particle must
  cross, and the maximum full-step lab displacement must remain before fringe.
* Native full-plane E and B differences at common lab times without temporal
  resampling, separately for all output and the prefix through the full cohort's
  detector passage. The PNG/PDF shows fixed/moving traces and explicit differences.
  Region metrics cover all pixels, the physical xy aperture excluding CPML, and
  a static two-sigma ellipse around the release centroid. Each reports regional
  L2, amplitude-squared area sums, peak timing, component errors, the worst local
  pixel/time and its local vector relative error. The previously identified
  sparse maximum (sample89, x=-37.5um, y=-22.5um) is inspected in addition to each
  new case's own maxima, so an error that moves cannot disappear from the review.
  E and B are always reported separately in their own SI units. No unscaled
  mixed E/B norm is used; a future combined norm would have to use E and cB.
  Native B remains staggered, so E^2/B^2 sums are not Poynting power.
* Existing ledger values plus separately summed inserted/cropped moving-volume
  field energy and exit/crop charge. A stationary-volume ledger residual alone
  is not a conservation or radiation measure after changing the volume.
* Optional same-MPI old/new fixed-path replay requires bitwise scientific
  datasets. Cross-layout differences remain separately reported and do not
  redefine the fixed-to-moving physical tolerance.

The analyzer does not invent a physical PASS tolerance or claim production
speedup. Initialization plus short diagnostic-rich runtime is not steady-state
throughput. First inspect structural checks, then judge approximation scales
against the accepted held baseline. Actual outcomes belong in a separate
RESULTS file; this versioned template intentionally contains no unmeasured PASS.
Dense-vs-sparse review must retain both the fixed-site and global maximum error,
separate pre/post full-cohort intervals, and compare particle errors against the
measured held self-field response: ballistic-corrected widths/covariances,
per-ID delta-K mean/RMS, and changes of energy spread. Never normalize a weak
delta by total beam energy or total beam width. A smaller integrated field error
alone does not justify discarding a local field discrepancy. Density-dependent
change is evidence for investigation, not by itself a proof of the mechanism.
Current limitations include frozen held Cowan0, current_filter:none, no checkpoint
or restart, no particle-retirement layer, no incident wave, and no transverse
ballistic/source-padding extension. Particle storage is not optimized here.
