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
| Maxwell/particle loop | Implemented for the no-seed path | Direct SI E/B Cowan-z or Yee update, Boris push, charge-conserving current deposition and MPI slab migration. |
| Field boundary | Implemented for no-seed Cowan runs | Compact CFS-CPML with checked geometry and particle-carrier treatment. Seed-wave TF/SF plus Cowan/CPML is not implemented. |
| Particle boundary | Implemented | Physical histories stop at CPML entry; output-free carriers damp current and export residual charge at the outer face. |
| Stops and partial output | Implemented | Element-aware and reference-centre stops; interactive SIGINT/SIGTERM closes readable incomplete output, throughput avoids that polling. |
| Lab detector planes | Implemented | Rank zero alone writes buffered fixed-lab field and particle planes; disabled detectors allocate and communicate nothing. |
| Particle trajectories | Implemented as an optional diagnostic | Complete enough for small-particle far-field validation, but intentionally unsuitable as the primary `10^7`-particle radiation record. |
| Trajectory radiation analysis | Implemented for small tests | Angular/integrated spectra, polarization/Stokes and optional window/ensemble cross-spectral density are available. |
| Field-plane background tools | Partially implemented | Uniform-motion background reconstruction and matched-baseline band power/energy comparison exist. |
| Full field-plane radiation analysis | **Missing** | There is no durable detector-to-full-spectrum, mutual spectral density and spatial-coherence pipeline yet. Existing power comparison is not a substitute. |
| Initial particle self-field | **Missing** | The Maxwell state starts without a Gauss-consistent bunch self-field. Startup transients may contaminate quantitative production results. |
| Particle subcycling | **Missing** | Particle and field steps are identical. The new undulator sampling guard can force a fine field grid but cannot reduce that cost. |
| Laser/seed injection | **Missing on the target kernel** | The generalized Cowan/CPML TF/SF injection remains the principal unimplemented source path. |
| Particle retirement | Experimental | It is not exactly charge-continuous and requires a matched zero-radiation baseline plus convergence tests. |

The release therefore is **not** complete apart from laser injection. It is a
coherent no-seed development solver with working diagnostics, but the three
other bold missing areas are production gates separate from laser injection.

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
- at least `mesh.particle_steps_per_undulator_period` actual steps across the
  shortest undulator period, computed from the fastest loaded particle;
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
- absence of a Gauss-consistent initial particle self-field;
- absence of particle subcycling and experimental retirement use.

Detector sample cadence is stored exactly and band tools reject requests above

A global particle-loss-versus-field-energy equality is deliberately not a
startup hard stop. It depends on detector aperture, radiation through other
boundaries, prescribed-field work and matched particle crossings. The particle
plane contract instead documents per-particle gamma differencing followed by
compensated or extended-precision accumulation.
their Nyquist limit. The simulator cannot itself prove that an unspecified
future analysis band, transverse aperture, CPML reflection level or macro-
particle representation is converged; those remain required parameter scans.

## Resource report

When `runtime.resource_monitor.enabled` is true, rank zero writes flushed
`[resource]` lines suitable for a redirected Slurm log. The startup report
contains:

- current and peak resident memory measured from the operating system;
- modeled maximum-per-rank and aggregate allocated memory, including field,
  CPML, halo, particle-capacity, detector and trajectory buffers;
- upper estimates for uncompressed field-plane and trajectory output;
- a short non-mutating step benchmark and a safety-factored wall-time upper
  estimate when the no-seed calibration is available.

Periodic reports contain progress, elapsed time, current/peak resident memory
and an upper remaining-time estimate. The final report records measured loop
seconds per step and full wall time. The model intentionally does not claim to
predict filesystem contention, scheduler interference, later migration peaks
or early physical stopping. Its two safety factors must be calibrated for the
target machine, MPI decomposition and representative input.
