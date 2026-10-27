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
  aprl_held_cut_probe OUT STEPS PARTICLES PREWARM live|frozen
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
max(90,STEPS) and maximum allowed prewarm64; the width remains208. At64/90steps
the full grid is624, initial window offset248. This exactly retains the earlier
64-step fixture while allowing90steps/64prewarm with proven dynamic clearance.
Prewarm variation leaves release positions/momenta and the full/cut geometry
fixed. Do not call it a scan of production physical guards or solver parameters.

run.py prepare freezes binaries, source snapshots, production card/log,
moment sample and scripts into a NEW directory. run.py run refuses reuse,
checks hashes and physical-core binding, and executes cases serially. gnode1
mode is the authorized CPU-only1MPI<=16threads exception,180s/case/2h total.
Slurm mode reserves one exclusive128CPU node, uses1MPIx16threads per step,
600s/case/25min total; run.sbatch requests no specific node and no requeue.
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

Results: /public/home/wfollowm/APRL_Workspace/results/held_cut_20261028/
