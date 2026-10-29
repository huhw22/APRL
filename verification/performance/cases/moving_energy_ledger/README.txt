Moving-window energy-ledger acceptance (2026-10-29)

Scope: diagnose the existing mover/reference supply without modifying Cowan,
Boris, current deposition, CPML or particle release. Old binary is the accepted
pending-admission cde740f... binary. No production input or full simulation.

prepare.py --workspace /public/home/wfollowm/APRL_Workspace --root NEW_DIRECTORY
           --candidate ABSOLUTE_NEW_BINARY
does not submit jobs. It freezes the 31,675-record short input, copies scripts,
records binary/runtime/input/card SHA256 and creates three independent stages:

science: one exclusive node, 16 MPI x 8 threads. One-period undulator and free
drift, each old-on/new-on/new-off/old-off. All 257 physical steps, complete
release and moving exit, raw field/three particle planes/per-rank trajectories.
Output payloads must be bitwise identical within each physical case/layout.

performance: starts only after the same-binary science gate. Same short source
and geometry, diagnostics every 10 steps, no trajectories/particle planes;
raw field plane retained at 1/10 cadence. Three repetitions of all four variants
in alternating order for each magnet case. Profile is off. Every pair is kept;
report 40-240 step throughput medians, ranges, individual paired percentages,
and old-on/new-on incremental cost separately from complete ledger on/off cost.
These timings do not establish full-production overhead or acceleration.

mpi: starts only after science. Two exclusive nodes, same total16x8 layout
(8 ranks per node; 64 used physical cores per reserved node). Old-on/new-on
for both physical cases. Compare only within identical two-node layout.

The source card is the already exercised dense-two-node-r3 small fixture:
160x160 transverse grid; 896 reference cells,289 active longitudinal cells;
64/128 guards; original fine dz/boost; transverse CPML only; no snapshot.

run.sbatch ROOT science|performance and mpi.sbatch ROOT are submitted by the
root agent only. There is no automatic queue submission. Existing watchdog
600/60 seconds and job-local self,vader,tcp settings remain; no system edits.

Review separates three questions:
1) physical-output invariance and complete lifecycle (hard checks),
2) independent accounting identities with a source-scale floating-point bound,
3) actual energy imbalance (reported raw/boundary-only/source-corrected; no
   invented physical tolerance or automatic closure_valid override).

Frozen particles and post-exit ballistic carriers can supply current while
their accounted kinetic energy stays constrained. Their diagnosed supply work
is an approximation and is kept distinct from prescribed analytic-device work.
No-magnet prescribed work must be zero; undulator boosted-frame work must be
nonzero. Physical residual is reported separately before/after first exit.

HDF5 completeness, committed extents, all numeric finite values, MPI finalization,
binding, binary/runtime identity and scalar formula reconstruction are checked.
FAIL preserves all evidence. No automatic retries, tolerance changes, production
launch, default binary replacement or extra convergence cases are performed.
