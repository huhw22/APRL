# Current release audit and production gates

This document separates implemented paths from planned work. A successful
smoke test is not treated as evidence that every requested radiation product
is production-ready.

## Capability matrix

| Area | Current status | Production note |
|---|---|---|
| Particle input | Implemented | Parallel HDF5 v2 supports per-record relative weights; v1 remains equal-weight compatible. A bounded-memory text converter is provided. |
| Generated beam | Implemented for tests | Deterministic, uncorrelated Gaussian with total electrons and macro-particle count; it is not a beam-preparation model. |
| Relativistic transform | Implemented and audited | Free-drift simultaneity placement, SI E/B transforms and light-front longitudinal momentum transforms are checked by a lab/boost round trip. |
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
| Particle/field energy closure | Boosted runtime ledger plus laboratory validation report implemented | The optional rank-zero HDF5 ledger contains boosted-frame particle kinetic energy, interior E/B energy, six inner-CPML fluxes, removed-particle energy, prescribed work and energy-spread moments. The downstream report stays entirely in the lab observer, reports the `K=0` collective/numerical control, mean/spread changes, forward radiation, and an explicitly unresolved remainder. A fast 50 A analytical scale tool is included. The first dense driven control remains unconverged, so neither report alone is a production approval. |
| Initial particle self-field | Implemented with boundary caveats | A distributed CIC/Poisson projection enforces the interior discrete Gauss law and releases all potential state before E/B advance. Box padding and the near-rest-frame electrostatic assumption still require convergence. |
| Particle subcycling | Implemented for prescribed devices | Automatic Boris substeps resolve analytical laboratory devices. Grid E/B, current deposition and detector cadence remain on the Maxwell step, which still sets radiation bandwidth. |
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
- convergence and post-check residual of the enabled CIC/Poisson initial
  Gauss-field projection, with all CIC charge confined to interior
  zero-potential-solver vertices;
- enough automatically selected Boris substeps to satisfy
  `mesh.particle_steps_per_undulator_period` across the shortest undulator for
  the fastest loaded particle, without exceeding
  `mesh.maximum_particle_substeps`;
- detector, retirement, magnetic interaction, stop and frequency-protection
  placement rules documented by the input-card specification.

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
- disabling the Gauss-consistent initial field, or using its electrostatic
  approximation in a boost frame far from the bunch mean rest frame;
- experimental retirement use.
- an enabled runtime energy ledger whose final residual relative to initial
  boosted kinetic-plus-field energy exceeds its configured warning tolerance.

Detector sample cadence is stored exactly and band tools reject requests above
their Nyquist limit. The simulator cannot itself prove that an unspecified
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
