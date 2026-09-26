# Release readiness and production requirements

This document defines the implemented scope, unresolved limitations, and
requirements for production calculations. A successful smoke test establishes
software operability only; it does not establish physical convergence of every
radiation product.

## Capability matrix

| Area | Current status | Production note |
|---|---|---|
| Particle input | Implemented | Parallel HDF5 v4 preserves fixed-Elegant-plane arrival-time events and synchronizes them with each particle's own velocity; v3 offset records and legacy v1/v2 snapshots remain readable. Relative weights are supported. The official-library C++ SDDS converter performs direct page/ID-aware conversion without an intermediate text dump; the older bounded-memory text converter remains available. |
| Generated beam | Implemented for tests | Deterministic, uncorrelated Gaussian with total electrons and macro-particle count; it is not a beam-preparation model. |
| Relativistic transform | Implemented and audited | Boosted time zero is anchored to the reconstructed lab bunch front, the transformed bunch is centred in the numerical box, and the large virtual lab-event span is diagnostic rather than an entrance-drift constraint. SI E/B and light-front momentum transforms are checked by a lab/boost round trip. |
| Mesh geometry | Implemented and audited | YAML supplies integer `cells` and physical `cell_size`; extents use multiplication and MPI slab offsets/counts remain integer. Old length/resolution inference is rejected. |
| Maxwell/particle loop | Implemented for the no-seed path | Direct SI E/B Cowan-z or Yee update, Boris push, charge-conserving current deposition and MPI slab migration. Weak-collective particle/radiation energy agrees in sign and scale; the dense-bunch ledger is measurable but not yet converged. |
| Field boundary | Implemented for no-seed Cowan runs | Compact CFS-CPML with checked geometry and particle-carrier treatment. Seed-wave TF/SF plus Cowan/CPML is not implemented. |
| Particle boundary | Implemented | Physical histories stop at CPML entry; output-free carriers damp current and export residual charge at the outer face. |
| Stops and partial output | Implemented | Element-aware and reference-centre stops; interactive SIGINT/SIGTERM closes readable incomplete output, throughput avoids that polling. |
| Lab detector planes | Implemented | Rank zero alone writes buffered fixed-lab field and particle planes; disabled detectors allocate and communicate nothing. |
| Particle trajectories | Implemented as an optional diagnostic | Complete enough for small-particle far-field validation, but intentionally unsuitable as the primary `10^7`-particle radiation record. |
| Trajectory radiation analysis | Implemented for small tests | Angular/integrated spectra, polarization/Stokes and optional window/ensemble cross-spectral density are available. |
| Field-plane background tools | Implemented with model choices | Uniform-motion background reconstruction and matched-baseline amplitude subtraction/power comparison exist; retirement remains experimental. |
| Full field-plane radiation analysis | Implemented for downstream forward modes | Threaded FFTW analysis provides angular/integrated spectra, Stokes data, Hann-window coherent/fluctuation splitting, global transverse coherence/Gram eigenvalues, selected spatial CSD and two-frequency CSD with memory/output guards. |
| Particle/field energy closure | Boosted runtime ledger plus laboratory validation report implemented | The optional rank-zero HDF5 ledger contains boosted-frame particle kinetic energy, interior E/B energy, six inner-CPML fluxes, removed-particle energy, prescribed work and energy-spread moments. The downstream report stays entirely in the lab observer and can combine colocated entrance/exit raw signed Poynting fluxes with the matched `K=0` particle control. A three-period dense test reached 99.23% longitudinal raw closure and 91.73% in its selected forward band after detector-window convergence. A fast 50 A analytical scale tool is included; aperture, mesh, CPML and longer-device convergence remain production gates. |
| Initial particle self-field | Relativistic rigid-beam initialization implemented; convergence still required | A distributed CIC/Vay relativistic-Poisson solve initializes E and current-consistent B, enforces the interior discrete Gauss law, and releases all potential state before E/B advance. CPML starts outside a zero-potential static boundary at its inner surface, preventing Coulomb-tail/CPML startup flow. Velocity spread, macro-particle sampling and distance to that boundary still require convergence. |
| Particle subcycling | Implemented for prescribed devices | Automatic Boris substeps resolve analytical laboratory devices. Grid E/B, current deposition and detector cadence remain on the Maxwell step, which still sets radiation bandwidth. |
| Radiation-band preflight | Implemented as an optional on-axis guard | A user-declared maximum laboratory photon energy is Lorentz-transformed into the box frame. Longitudinal grid, Maxwell/current and field-plane cadence are hard-rejected only below strict Nyquist; configurable quality margins are warnings. Transverse envelope, aperture, CPML and macro-particle convergence remain explicit physics studies. |
| Output safety and provenance | Implemented and release-gated | Existing APRL and post-process outputs are refused unless the relevant card explicitly enables overwrite. A root-only manifest preserves the exact card, source/build/MPI and particle-input identity; one `run_id` and configuration digest link all APRL HDF5 products. |
| End-to-end post-processing regression | Implemented and required | Lightweight simulations drive reconstruction, field-plane analysis, trajectory radiation, matched power comparison and lab-frame energy closure; completion markers and finite key outputs are checked. All post-processing YAML mappings reject unknown keys. |
| Laser/seed injection | **Missing on the target kernel** | The generalized Cowan/CPML TF/SF injection remains the principal unimplemented source path. |
| Particle retirement | Experimental | It is not exactly charge-continuous and requires a matched zero-radiation baseline plus convergence tests. |

The remaining missing functional path is generalized laser/seed injection on
the Cowan/CPML target kernel. The no-seed solver and field-plane diagnostics
are implemented, but production conclusions still require the convergence and
model checks listed below. In particular, dense-bunch global energy closure is
now measurable rather than missing, but remains an explicit numerical
production gate; experimental retirement is not promoted by this capability
audit.

## Numerical preflight policy

The following conditions are hard errors because continuing would knowingly
violate the selected model or input contract:

- Cowan-z geometry/stability restrictions and CPML/MPI layer geometry;
- positive integer mesh counts, positive finite cell sizes, finite multiplied
  extents and at least two longitudinal cells on every MPI rank;
- initial bunch placement outside the first element interaction region and
  inside the boosted longitudinal box;
- positive finite particle weights, preserved electron charge-to-mass ratio,
  and global charge normalization to `beam.input.electrons`;
- finite Lorentz transforms with relative lab/boost round-trip momentum and
  gamma error no greater than `1e-10`;
- convergence and post-check residual of the enabled
  CIC/relativistic-Poisson initial
  Gauss-field projection, with all CIC charge confined to interior
  zero-potential-solver vertices;
- enough automatically selected Boris substeps to satisfy
  `mesh.particle_steps_per_undulator_period` across the shortest undulator for
  the fastest loaded particle, without exceeding
  `mesh.maximum_particle_substeps`;
- when a target `radiation_resolution` is enabled, strict forward on-axis
  Nyquist sampling by the z grid, Maxwell/current step and every field plane;
- detector, retirement, magnetic interaction, stop and frequency-protection
  placement rules documented by the input-card specification.
- a strict YAML key schema at every mapping level, including disabled optional
  blocks, so misspelled settings cannot silently fall back to defaults;
- a finite stop target: element-free `after-last-element` cards are rejected,
  while detector-only and element-free `reference-center-z` propagation remain
  supported.

The sampling error reports a maximum `dt` and, for Cowan-z, maximum `dz` in SI
and card units. Geometry errors similarly report a corrected grid, initial
beam reference or detector position rather than only rejecting the card.

The following are warnings because their acceptable magnitude depends on the
physics target:

- a Lorentz round-trip error above `1e-13` but below the hard threshold;
- the estimated individual-electron lab-energy roundoff scale
  `epsilon * gamma_max * m_e c^2`, with advice to compute small beam-energy
  changes from per-particle gamma using compensated or extended-precision
  accumulation instead of subtracting two rounded total beam energies;
- disabling the Gauss-consistent initial field, or using its legacy
  electrostatic approximation in a boost frame far from the bunch mean rest
  frame;
- experimental retirement use.
- an enabled runtime energy ledger whose final residual relative to initial
  boosted kinetic-plus-field energy exceeds its configured warning tolerance.
- `mesh.duration` shorter than the step-rounded startup estimate. The estimate
  is exact for the inertial reference-centre stop and ballistic-only for the
  initial particles in `after-last-element` mode.

Detector sample cadence is stored exactly and band tools reject requests above
their Nyquist limit. APRL cannot itself prove that an unspecified
future analysis band, transverse aperture, CPML reflection level or macro-
particle representation is converged; those remain required parameter scans.

A global particle-loss-versus-field-energy equality is deliberately not a
startup hard stop. It depends on detector aperture, radiation through other
boundaries, prescribed-field work, bound-field energy and matched particle
crossings. The particle-plane contract instead documents per-particle gamma
differencing followed by compensated or extended-precision accumulation, and
the standalone closure report makes the limitations explicit. The controlled
one-electron-equivalent test has the correct sign and scale: a 50--100 eV,
`+/-0.26 mrad` trajectory far field contained 48.13% of the kinetic loss, and
`+/-2 mrad` contained 54.88%. The `10^6`-electron test was dominated by a
different collective/near-field energy change between signal and `K=0`; it
cannot be reduced to kinetic loss equals one forward-plane flux. The runtime
ledger now supplies the missing stored-field, six-face flux, removed-particle
and prescribed-work terms. Its three-period driven control closed to 3.85% of
the summed exchange scale (7.996% of initial boosted kinetic-plus-field
energy), so dense-bunch mesh, macro-particle and initial-padding convergence is
still required.

## Resource report

When `runtime.resource_monitor.enabled` is true, rank zero writes flushed
`[resource]` lines suitable for a redirected Slurm log. The startup report
contains:

- current and peak resident memory measured from the operating system;
- modeled maximum-per-rank and aggregate allocated memory, including field,
  CPML, halo, particle-capacity, detector, trajectory and energy-ledger
  buffers, plus the
  transient four-slab initial Poisson solve when enabled;
- upper estimates for uncompressed field-plane and trajectory output;
- a short non-mutating step benchmark and a safety-factored wall-time upper
  estimate when the no-seed calibration is available.

Periodic reports contain progress, elapsed time, current/peak resident memory
and an upper remaining-time estimate. The final report records measured loop
seconds per step and full wall time. The model intentionally does not claim to
predict filesystem contention, scheduler interference, later migration peaks
or early physical stopping. Its two safety factors must be calibrated for the
target machine, MPI decomposition and representative input.
