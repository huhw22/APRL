Fixed-width boosted device-window feasibility probes — 2026-10-27
================================================================
User requested theoretical derivation and bounded experiments, not deployment.
THEORY.txt derives the fixed-lab-device strip and identifies the unresolved
exterior-history cost. Production core sources and running jobs stay unchanged.

Build independently with the accepted MPI/parallel-HDF5 toolchain:
  -DAPRL_BUILD_WINDOW_FEASIBILITY=ON -DDIRECT_EB_BOUNDS_CHECK=ON
Targets: aprl_pic_dependency_probe, aprl_window_streaming_probe,
         aprl_window_analytic_probe.
Keep consistent bounds checking across the core and consumers. No mandatory
production regression group or production configuration is changed.

Experiments
-----------
1. pic_dependency_probe: 48 samples with one real-electron step, three subcell
   positions, one/five Boris substeps and three perturbation amplitudes. Vacuum,
   zero perturbation and a farther source are negative controls. An independent
   first-order Lorentz/CIC/Ampere prediction checks the two-cell response. This
   is stencil support, not a physical signal-speed test or equilibrated beam.
2. window_streaming_probe: NxNy=6x5, window Nz48; beta0.9 and the production beta,
   durations48/96/192. Actual Cowan updates, exact reference E halo/end-face/new
   layer input. Never overwrite retained interiors. Full-reference generation
   is explicitly counted. Six-component general manufactured fields/current
   exercise the discrete operators; this is not a Gauss-compatible beam model.
3. window_analytic_probe: NxNy=8x7, window Nz48, same speeds and durations. Known
   discrete +/-z vacuum waves with two polarizations, constant longitudinal
   fields, strong amplitude1 and weak amplitude1e-6. Boundary input comes from
   analytic formulas, not a full simulation. Only the outer lateral component
   layer, E z end faces, E halo and entering layer are supplied. The entire
   6x5 interior plane is measured. This is a transversely uniform vacuum test,
   not a generic angular-spectrum, held source or space-charge validation.

All native Yee E/B space/time staggering is preserved. Plane outputs are
transformed to lab with the existing transform; B is reported multiplied by c
for common V/m units. This checks preservation, not correction of native timing.
The missing-boundary controls use local PEC and zero incoming cells; they test
that missing correct boundary information is harmful, not which individual
boundary change causes the error.

Resources and review
--------------------
run_gnode.sh BUILD FRESH_RESULT runs probes1/2. Add --analytic-only for probe3.
Requires actual gnode1, explicit APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1, no fake Slurm
identity, existing shared lock, binding1MPIx1thread and180s/probe timeouts.
The selected executable hashes and all failures remain in fresh result folders.
The build's runtime-libs may contain only the existing accepted FFTW copies.

review.py CAMPAIGN independently reads all plane CSVs, checks raw and weak
fields, FFT complex amplitudes, phase and180x180 pixel/polarization products at
every frequency. No frequency filtering, intensity subtraction or demeaning.
These are single-realization cross spectra, not an ensemble coherence claim.
It recomputes the PIC analytic response and the operation counts, and writes
SCIENCE_REVIEW.json. plot_review.py CAMPAIGN plots counted work and measured
analytic-wave errors; it does not report production acceleration.

Accepted first-stage evidence
-----------------------------
Workspace results/window_streaming_20261027/: PLAN.json, THEORY.txt,
RESULT_REVIEW.txt, SCIENCE_REVIEW.json, feasibility.png/pdf, PROVENANCE.json.
probe-r1 and analytic-r1 contain frozen sources, binding, hashes, raw JSON/CSV
and logs. All probes passed on one physical core ofgnode1. Oracle checked
12,905,052 state values with zero bit mismatches;122,040 detector scalars also
match exactly. Analytic weak-field maximum relative error1.94e-8, below1e-5;
full six-component weak FFT relative L2<=9.24e-9, cross-spectral relative
Frobenius<=2.66e-9. This supports field-window/interface research only.

Not accepted or implemented: affordable unknown exterior/held boundary supply,
self-consistent particle admission in a moving window, moving CPML, distributed
repartition, streaming checkpoints, production speedup or end-to-end O(nlogn).
The next gate is independent, constraint-consistent supply for a known moving
charge distribution plus a weak wave, followed by actual held current crossing.

Release-state startup study (user clarification, 2026-10-27)
---------------------------------------------------------
The physical target is now explicitly quasi-steady accuracy, not reproduction
of the old held transient. Prescribed release x/p must stay fixed. An initial
self-field must preserve normal space-charge expansion without imposing an
excessive field-building energy transfer on the released particles.

aprl_startup_self_field_probe uses the same accepted core, generated-gaussian
seed710271, 4096 macro particles, -1pC, rest sigma4um, beta0 or0.23. All particles
are released at t0. Same x/p/Q/m across modes and meshes; no held interval,
external device, moving boundary or production input. base48^3 at1um,
wide64^3 at1um, fine96^3 at0.5um, all15um/c physical duration. Fine also halves
dt; fixed particle count is not a particle-number convergence study.

CLI: OUT MODE BETA GRID; MODE=qs,zero,scale0.9,scale1.1,fieldonly;
GRID=base,wide,fine. fieldonly requires beta0 and keeps all sources stationary
with J=0. Zero/scaled fields violate initial Gauss for the same rho and are
negative controls, not candidate physical models. qs uses ordinary production
relativistic-Poisson at1e-10 tolerance; core physics is unchanged.

run_startup.sh BUILD FRESH_RESULT MODE BETA GRID runs on explicitly authorized
gnode1 only, after1x1 core binding, with180s per case and the shared CPU lock.
review_startup.py ROOT independently checks initial/final particle records,
stable deltaK, Gaussian force moments, constraints and observations.
plot_startup.py ROOT --require-complete writes the eight primary-case figure.

Results: results/startup_self_field_20261027/{PLAN.json,REVIEW.json,
RESULT_REVIEW.txt,startup.png/pdf}. r1 contains ten coupled/control cases and
its frozen pre-fieldonly binary/source. r2 contains three zero-current static
source isolations, separately planned after energy diagnostics failed.
All cases finish and functional checks pass; physical energy acceptance FAILS:
moving-beam deltaK changes18.3% on refinement versus the predeclared5% screen,
and standard field-energy closure fails the10% exchange-scale screen. Normal
expansion is restored; zero-field moving startup gives large unphysical loss.
Do not relabel functional checks as physical PASS or infer production errors.

Ordinary-Poisson E=-grad(phi) satisfies Gauss but generally not Cowan's modified
curl(E)=0. With particles fixed and J=0, field state and standard field energy
still change; the standard energy is not an exact Cowan invariant. Its change
is not by itself proof of real physical energy creation. Particle-grid work,
field-grid work, native and centered B diagnostics must remain separate.
This study has not implemented/accepted reference injection, a discrete steady
moving-beam solution, streaming PIC, energy-precision production acceptance or
an asymptotic speedup. The next research target is a quasi-steady reference
compatible with the actual field/deposition operators, then release and weak
wave tests, rather than a longer window or exact old-held replay.
