Held inlet / right causal cut / transverse CPML validation

User-approved 2026-10-28 replacement for the longitudinal moving-CPML research
path. These are optional verification executables; no production src file,
default executable, checkpoint or MPI Simulation integration is changed.

Build in a new directory with the workspace MPI/parallel-HDF5 environment:
  cmake -S source/aprl -B build/held-cut-NEW -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpic++ \
    -DCMAKE_PREFIX_PATH=/public/home/wfollowm/aprl_dev/aprl-runtime-lab \
    -DBUILD_TESTING=ON -DAPRL_ENABLE_OPENMP=ON -DAPRL_BUILD_HELD_CUT_PROBE=ON \
    -DAPRL_BUILD_POSTPROCESSORS=OFF -DAPRL_BUILD_SASE_BUNCH_GENERATOR=OFF
  cmake --build build/held-cut-NEW --target aprl_held_cut_probe \
    aprl_held_cut_boundary_checks aprl_held_cut_causality_checks -j8

Scientific constants are anchored to the real 18932 card: gamma929.1305493,
dx=dy5um,dz1.25432624um, cdt=dz,160x160 transverse mesh,10-cell xy CPML,
order3/kappa8/alpha0/target1e-8,64/128 normal longitudinal guards, original
undulator and 80-samples/period-derived Boris substeps. No0.8c substitute.
The generated Gaussian is shorter/fewer particles; transverse sizes use full
production statistics, lab momentum/correlations use the reproducible modest
input sample (sample_source.py). It is a statistical surrogate, not the complete
ELEGANT microstructure or a claim of radiation convergence. Charge scales with
longitudinal sigma to keep the Gaussian peak line density.

The only grid source adapter is HeldCutBoundary. It prescribes the lower E
halo at native state n and lower Ex/Ey endpoint at n+1. New slices receive E/B
and eight transverse CPML history slabs at the same native state. Existing
values move by world position. z CPML histories are absent, z PEC flags false;
the right halo is explicit, not an implicit mirror. Right integer faces and
the half-step sacrifice collar cannot be used by particle/field sampling.
Original Cowan/Boris/current-deposition and CPML equations remain unchanged.

Operator checks deliberately corrupt the right E/B terminal faces and all
upper halos. Twelve successive steps plus left cuts must leave every trusted
field/history value bitwise equal to a full domain. Missed-shift and missing
incoming-history negative controls must fail equality. Small transverse arrays
in these indexing checks are explicitly NOT scientific geometry substitutes.

Causality checks independently poison both remote endpoints, all six native
field components, eight transverse history slabs and electric halos for40steps
per side. Complete-step lower influence reaches native layer n+1 (first step
can reach2); upper influence grows by at most n. The component-staggered
recurrence proves this envelope for arbitrary n, with tests exercising both
fronts. This is dynamic propagation only, not independence of the initial
elliptic solve from its finite-domain boundary.

Probe CLI:
  aprl_held_cut_probe OUT STEPS PARTICLES PREWARM live|frozen \
    [WARM_HORIZON [INITIAL_MODEL [BACKEND [BLOCK_LAYERS]]]]
Defaults remain64 and relativistic-poisson; BACKEND defaults to cg and
BLOCK_LAYERS to16. The bounded optional horizon is
64..256; the only additional initializer is the already implemented
cowan-relativistic-poisson. Its selection is an explicit research variant,
not a change to the production card or propagation kernels. BACKEND is cg or
blocked; blocked requires the optional build described below, Cowan initial
model, frozen mode and PREWARM0. The probe bounds BLOCK_LAYERS to1..128.
All science cases contain BOTH actual full and cut PIC worlds, with the same
prepared initial fields and identical prescribed release states.
  live: an independent held-source full grid evolves the reference; it never
        reads future actual PIC data. Its cost is included, so this is an
        interface correctness anchor, not a cost-saving production method.
  frozen: store one short-held-prepared native field/history state, translate
          it at mean held velocity, and use that approximate inlet only.
          No full reference evolution is performed after preparation.

The actual production activation plane remains lab z=-.05m; the actual4.5m
detector is not reached. A separate auxiliary local plane16dz downstream at
equal boost time is used for interface diagnostics, not silently substituted
for the production detector. Original prescribed magnet fields are included
when the short trajectory reaches the fringe. Raw E/B records have different
native time levels; B is sampled at the E-plane's boost position, not at its
own fixed-lab event. They support matched component comparisons but are not
simultaneous laboratory E/B or a Poynting-flux measurement. Plane pixels are
sparse, not a full-plane power integration.

For paired prewarm tests the full geometry is selected from a shared horizon
max(90,STEPS) and the selected prewarm horizon (default64); the width remains208. At64/90steps
the full grid is624, initial window offset248. This exactly retains the earlier
64-step fixture while allowing90steps/64prewarm with proven dynamic clearance.
Prewarm variation leaves release positions/momenta and the full/cut geometry
fixed. Do not call it a scan of production physical guards or solver parameters.
With WARM_HORIZON256 the90step complete reference has1008 cells and offset440;
all64/128/256prewarm and Cowan0 variants share that geometry. This is reference
padding for the experiment, not expansion of the208cell candidate window.
prepared-reference.csv records sampled native fields at time0 for same-source,
same-phase preparation comparisons. It does not replace future-inlet diagnostics.

run.py prepare freezes binaries, source snapshots, production card/log,
moment sample and scripts into a NEW directory. run.py run refuses reuse,
checks hashes and physical-core binding, and executes cases serially. gnode1
mode is the authorized CPU-only1MPI<=16threads exception,180s/case/2h total.
Slurm mode reserves one exclusive128CPU node, uses1MPIx16threads per step,
600s/case/25min total; run.sbatch requests no specific node and no requeue.
The optional --case-timeout-s 1200 sets a MAXIMUM for each Slurm case, not a
separate1200s allocation per case. The total remains1500s/25min. Before each
case, run.py uses min(case maximum, remaining campaign time minus10s); the
second/later case therefore receives a shorter limit when necessary. A case
is not started with less than60s remaining. Binding/checks consume the same
campaign time; default600 and gnode1's180 remain unchanged. Do not expand the
total limit to guarantee completion of a pair; report an incomplete pair.
The single-rank production Poisson initializer is serial even when the field
steps use16threads. Do not relax its tolerance to fit the gnode1 time limit.

review.py audits production/source constants, initial/release identity,
particle/event membership, stable lab energy changes, centered beam sizes,
native component plane differences, and measured cost components. Source
velocity correlations can cause ordinary beam focusing; total RMS change
is not isolated space-charge expansion. All-source averages and released-only
averages are both reported. No arbitrary scientific PASS tolerance is imposed.
Compare field and per-particle differences before deciding usability.

Source/initial/prewarm/reference evolution/entry preparation/cut advance/remap
costs are separate. Diagnostics and output are research overhead. Never quote
only the cut advance as whole-method speedup or claim O(N) from a live full-grid
reference. Whole-undulator weak-radiation phase, true4.5m plane, long-time
particle exit, multi-rank exchange and production checkpoint remain separate
acceptance stages.

Particle-sampling limitation: the64particle shortened source preserves scaled
charge and Gaussian line density, but represents about494.92 electrons per
macro-particle; the production input represents1 electron per record. Matching
its weight at this shortened charge requires about31675 particles. The sparse
fixture can validate interfaces and compare reference models on identical
sources; it cannot alone accept or reject production quasi-steady accuracy.
Possible CIC granularity effects require a separate sampling convergence test.

The follow-up convergence study separates preparation sensitivity, subsequent
frozen prediction error, full-control sensitivity, and the extra effect of the
cut. Exact pointwise replay of one held history is not the sole physical goal;
particle energy/beam moments and persistent field disturbances must be measured.
No blanket physical-pass threshold is inferred from the existing held spread.

Results: /public/home/wfollowm/APRL_Workspace/results/held_cut_20261028/
Follow-up: /public/home/wfollowm/APRL_Workspace/results/held_quasisteady_20261028/
Its offline preparation/ballistic figures can be regenerated into a NEW directory:
  python3 convergence_review.py --root NEW_REVIEW_DIRECTORY \
    --campaign /public/home/wfollowm/APRL_Workspace/results/held_quasisteady_20261028/preflight \
    --campaign /public/home/wfollowm/APRL_Workspace/results/held_quasisteady_20261028/long-held \
    --campaign /public/home/wfollowm/APRL_Workspace/results/held_quasisteady_20261028/cowan
The default horizon256 excludes the legacy-geometry replay case. Optional
--geometry-control OLD_CASE_DIRECTORY reports finite-domain sensitivity
separately. This script reads existing data and never starts a simulation.

Fast-reference design follow-up (2026-10-28):
  FAST_REFERENCE_DESIGN.txt derives transverse DST plus longitudinal blocked
  elimination, preserving distant source contributions without full-grid CG.
  modal_reference_check.py is a fixed small algebra check, not PIC or a speed
  benchmark. It uses the previous case's spacing/common beta, a full 3-D DST
  oracle, block sizes1/7/16/33, and missing-far-source negative controls.
  On authorized idle gnode1 after binding inspection, within the existing CPU
  test lock, use one physical core/one BLAS worker and timeout180s:
    python3 modal_reference_check.py --output /absolute/NEW/algebra.json
  In this Python algebra prototype, checkpoints still occupy O(A*Nz/block)
  RAM; it does not implement disk-backed summaries or native E/B queries.
  The separate C++ generator below adds those operations. Neither prototype
  establishes particle streaming or production integration.
  Results: results/fast_reference_design_20261028/.

Optional blocked Cowan reference implementation (2026-10-28):

Build aprl_blocked_cowan_reference and aprl_blocked_cowan_reference_checks by
adding -DAPRL_BUILD_BLOCKED_COWAN_REFERENCE=ON to the NEW verification build
configuration above. Keep APRL_BUILD_HELD_CUT_PROBE=ON for the paired PIC
probe. FFTW3 double-precision headers/library are required only when the
optional FFTW-dependent targets are selected. Explicit dependency paths are:
  -DAPRL_VERIFY_FFTW_INCLUDE_DIR=/absolute/private/include
  -DAPRL_VERIFY_FFTW_LIBRARY=/absolute/private/lib/libfftw3.so.3
Use the same MPI/HDF5 toolchain as the existing probe, and freeze these private
dependencies with the new binary. No system package, MPI/network setting or
production executable needs changing. Build targets, after configuration:
  cmake --build build/held-cut-NEW --target aprl_held_cut_probe \
    aprl_blocked_cowan_reference_checks -j8
Leaving the new option OFF preserves the previous build and its CG-only CLI;
even with it ON, omitting BACKEND still selects the original CG path. The
default-CG twelve scientific CSV files must replay the matched legacy case.

All generator calls currently require one live MPI rank. Reference queries
are serial and not thread safe; existing field/PIC work may still use the
recorded OpenMP thread count. The generator supports the same finite
Dirichlet Cowan initial operator with one mass-weighted common longitudinal
beta. That solver beta and the frozen-reference translation beta have their
existing separate definitions; neither is a replacement for the near-c
production boost/window speed. The particle release records are unchanged.

Examples below describe executable arguments, not permission to run on the
login node. Use the existing bounded compute-node runner and binding checks:
  aprl_blocked_cowan_reference_checks /absolute/NEW_CHECK_DIRECTORY
  aprl_held_cut_probe /absolute/NEW_CG_CASE \
    90 64 0 frozen 256 cowan-relativistic-poisson cg
  aprl_held_cut_probe /absolute/NEW_BLOCKED_CASE \
    90 64 0 frozen 256 cowan-relativistic-poisson blocked 16
For run.py prepare, the paired --case specifications are:
  --case cg64:90:64:0:frozen:256:cowan-relativistic-poisson:cg
  --case blocked64:90:64:0:frozen:256:cowan-relativistic-poisson:blocked:16
Use --environment slurm --threads 16 with a fresh campaign directory.
The approved slower-initializer maximum
is --case-timeout-s 1200, with the unchanged total1500s rule above. A later
matched-weight short-bunch pair may use31675 particles and120steps in BOTH
backends; its geometry must match within that pair. It is a separate source
sampling case, so comparing its runtime with the64particle case does not
measure an implementation speedup.

The generator preserves all source contributions. It forms CIC contributions
in original-particle order, stably groups them by longitudinal source block,
then uses transverse DST and a longitudinal forward elimination/backward
reconstruction. It does not truncate remote charge at a guard distance. For
A transverse interior nodes, L longitudinal interior layers and block depth B,
rho.bin and phi.bin each store8*A*L bytes; checkpoints.bin stores two double
transverse arrays per block, about16*A*ceil(L/B) bytes. These files remain in
CASE/blocked-reference after the generator is destroyed. Their directory and
report must be new; existing output is refused rather than overwritten.
The files use native-endian doubles, not a portable snapshot/restart format.
Short reads, invalid lengths and nonfinite query data are checked; there is
no payload checksum, so finite-valued corruption is not generally detectable.

Solver workspace and query caches have fixed O(A*B) capacity when A and B are
fixed. Queries retain at most three potential blocks, one native E/B block
with an extra right plane for adjacent-plane interpolation, and a source
block only when scalar charge diagnostics request it. However, source
preparation still requires the caller's full particle input and O(P) temporary
CIC contribution indices, plus a conservative sorting-scratch allowance.
Those indices are freed before solve/query. This is not yet a fully streamed
particle reader, constant-memory whole simulation or distributed-I/O method.
Global preparation and scalar disk volume still grow with L; each inlet
query uses the prepared blocks instead of rescanning all particles.

In blocked mode, the runtime inlet queries the generator directly and has no
full-donor field/layout pointer. The probe nevertheless retains an independently
filled complete donor and complete full PIC world as comparison oracles, and
allocates a full rho array for the independent original-depositor Gauss check.
Their memory and elapsed time are verification costs, not evidence that the
whole probe has bounded memory. reference_field_bytes records the retained
donor. initialization_wall_s includes generator construction, donor fill and
the Gauss audit; blocked_oracle_donor_fill_wall_s and
blocked_independent_gauss_wall_s separately identify the latter two operations.
blocked_independent_gauss_bytes records that audit's extra rho array.

CASE/blocked-reference-report.json records preparationSeconds for generator
construction alone: geometry/modes, source grouping, forward/backward passes
and file checks. It excludes later native E/B queries, donor fill and Gauss.
sourceSeconds/forwardSeconds/reconstructionSeconds are nested subtimers, not
additional costs. querySeconds measures native-cache builds, excluding cache
hits/copies; it is not the complete inlet cost. The report is written after
the probe's queries, so binary I/O counters include diagnostics/oracle queries
as well as candidate use. Count source preparation, disk traffic, reconstruction,
and measured cut preparation/remap/feed/advance when discussing total cost.
Do not add overlapping timers or quote FFT kernel time alone as acceleration.

Byte fields report array capacities and binary I/O, not measured peak RSS.
sourceIndexBytes is temporary index capacity; sourceSortingAllowanceBytes is
a conservative extra allowance, not measured allocation. CacheCapacityBytes
fields retain their historical maxima; residentArrayBytes is the current
sum of persistent potential/source/E-B cache capacities. queryWorkspaceBytes
records the maximum separate temporary construction workspace. These belong
to different phases and must not simply be summed as a simultaneous measured
peak. FFTW/MPI/runtime private allocations are excluded; JSON I/O is excluded
from bytesRead/bytesWritten. No production memory saving follows from these
capacity counters alone.

Acceptance sequence and postprocessing:
  1. Independent serial generator checks compare original CIC, native six
     field components and Gauss against the existing initializer, including
     finite guards, block seams/endpoints, preserved particles and fault/cache
     controls. Small indexing geometries are not scientific production cases.
  2. Replay the legacy default CG science files, then compare same-source,
     same-geometry CG/blocked Cowan0 full/cut runs. Respect the finite CG
     tolerance; different solvers need not produce bitwise identical fields.
  3. Inspect prepared/native-plane differences and two separate comparisons:
     blocked minus CG at fixed full/cut world; cut minus its own matched full
     world within each backend. Solver agreement cannot erase an existing
     frozen-inlet approximation error.
  4. At the same particle IDs and prescribed release state, inspect lab-energy
     differences in eV (mean/RMS/spread), complete event membership and the
     beam-width/covariance corrections after subtracting analytic ballistic
     focusing. Verify that the whole event-interpolation step is outside the
     prescribed magnetic fringe. Do not normalize weak energy changes by the
     large total beam energy or use total focusing as the space-charge signal.

Run only postprocessing on completed paired cases, into a NEW review directory:
  python3 blocked_review.py --cg-case /absolute/CG_CASE \
    --blocked-case /absolute/BLOCKED_CASE --output /absolute/NEW_REVIEW
Optional --replay-case LEGACY_CG_CASE requires all twelve scientific CSVs to
be byte-identical; --held-case ORDINARY_HELD_CASE adds a matched-geometry
preparation-history comparison, without declaring held an exact truth.
--no-plot omits figures. Outputs are BLOCKED_REVIEW.json, two per-ID energy
CSV files and native-field/energy plus ballistic-response PNG/PDF figures.
Structural/residual PASS is distinct from physical acceptance; the reviewer
does not invent a weak-radiation or quasi-steady tolerance.

This documentation states the procedure and implemented scope, not an unrun
physical PASS. Actual results and limitations belong in the campaign's
RESULTS/RESULT_REVIEW records under results/blocked_reference_20261028/.
Single-rank fixed-source reference agreement alone does not establish long-time
stability, actual right-particle exit, the4.5m detector, radiation separation,
MPI/snapshot compatibility, production throughput or a whole-program O(n log n)
complexity bound. Production/default/src and the approved propagation kernels
remain outside this optional research implementation.
