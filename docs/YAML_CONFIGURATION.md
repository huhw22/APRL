# YAML input-card specification

APRL accepts one ordinary YAML document. Every mapping is checked
against the documented schema before physical parsing: an unknown key is a
hard error that reports its YAML line, mapping path, and the recognized keys.
This catches misspellings even inside a disabled optional block. The five
post-processing cards follow the same strict rule. Unknown physical defaults
are deliberately avoided; required quantities are likewise reported with
their YAML line when missing or invalid.

## Units

```yaml
units:
  length: micrometer
  time: picosecond
```

Length accepts `m`, `mm`, `um`/`micrometer`, and `nm`. Time accepts `s`, `ms`,
`us`, `ns`, `ps`, `fs`, and `as`. These units apply to numeric length and time
values in the YAML document. The core converts them once and then uses SI.
Particle HDF5 positions are always metres and do not inherit the YAML unit.

## Output safety and run identity

```yaml
output:
  overwrite: false
  manifest: output/run-manifest.yaml
```

This top-level block controls every APRL output. Existing files are a
hard error by default; `overwrite: true` is required for intentional
replacement. The root-only manifest preserves the exact card and records a
unique run identity, source/build information, MPI environment and particle
input identity. The same identity is embedded in all output HDF5 groups. See
[RUN_PROVENANCE.md](RUN_PROVENANCE.md) for the complete contract.

## Mesh and boost

```yaml
mesh:
  field_solver: cowan-z
  cells: [40, 40, 60]
  cell_size: [1.0, 1.0, 0.2]
  center: [0.0, 0.0, 0.0]
  duration: 0.01
  boost_gamma: 2.0
  particle_steps_per_undulator_period: 32
  maximum_particle_substeps: 4096
```

`field_solver` is required. Use `cowan-z` for the z-priority controlled-
dispersion stencil or `yee` for the original regression path. Both operate on
the same staggered E/B lattice; neither allocates A/phi. The Cowan path applies
the transverse smoothing directly while evaluating Faraday's law, rather than
storing three extra smoothed-field volumes.

For `cowan-z`, z must have the smallest cell size: `dx >= dz` and
`dy >= dz`. The program checks this before allocating fields or reading the
particle file, exits on violation, and reports the minimum valid transverse
`cell_size` derived from the configured `dz`. Its time step is
`dt = dz/c`, which makes resolved vacuum propagation exactly dispersion-free
on the z axis. Startup logging reports the squared aspect ratios, all Cowan
stencil coefficients, and transverse-axis phase/group velocity ratios at 16
cells per wavelength. Those transverse figures are diagnostics, not a claim
that an arbitrary oblique mode is dispersion-free.

The present TF/SF seed-wave boundary correction is still the Yee form, so a
configuration that combines `cowan-z` with a nonempty `incident_waves` list is
rejected before particle input. Use `yee` for seed-injection regression until
the generalized Cowan TF/SF correction is implemented.

`cells`, `cell_size`, and `center` describe the boosted computational box.
`cells` contains the authoritative integer counts and each direction needs at
least three cells; z additionally needs at least two cells per MPI rank.
`cell_size` is converted to SI once. The full physical extent is then derived
only by multiplication:

```text
extent_i = cells_i * cell_size_i.
```

The obsolete `lengths` and `resolution` keys are rejected rather than mixed
with the new contract. MPI z slabs are formed from the integer z count using
integer quotient/remainder offsets; ranks may own counts differing by one, but
no floating-point division determines a slab boundary. Startup prints the
counts, cell sizes, derived extents and per-rank z-count range.

`duration` is boosted-frame time and remains a hard upper guard, not a requested
physical propagation distance. After particles are transformed, startup logs a
`[duration-estimate]` line. For `reference-center-z` it uses the exact inertial
boost-reference worldline used by the stop predicate. For
`after-last-element` it extrapolates every initial particle ballistically to the
last interaction exit and takes the slowest crossing. The latter is an estimate
because magnetic forces can alter the orbit. If the step-rounded estimate
exceeds the configured duration, startup prints a value in the YAML time unit;
the run is not enlarged automatically.

`particle_steps_per_undulator_period` is the
minimum requested Boris sampling of the shortest configured undulator. After
the real particles are read and boosted, the program derives the largest
laboratory z advance made by any particle in one Maxwell step and chooses the
smallest integer number of particle substeps that meets the request.
`maximum_particle_substeps` is a positive safety/cost ceiling (default 4096).
Exceeding it is a hard preflight error with both the required ceiling and a
field-step/z-cell refinement suggestion.

The substeps resample analytical laboratory devices along the evolving
particle orbit. The staggered grid E/B sample is held fixed during the
enclosing Maxwell step, and charge-conserving current deposition and detector
output use the start-to-final-position field-step chord. Thus this setting
improves orbit integration through prescribed devices, but it does not raise
the Maxwell Nyquist frequency, resolve intra-step radiation current, or permit
a coarser radiation grid. The value `1` in the committed examples is only for
smoke testing. See `docs/PARTICLE_SUBCYCLING.md` for the numerical contract.

`boost_gamma` selects one constant inertial computational frame for the whole
run. It is editable per card but is not changed when particles enter or leave
an element: a time-dependent boost would be a non-inertial coordinate system
and would invalidate the shared Maxwell grid and laboratory element/detector
worldlines. For a planar undulator, static magnetic deflection preserves total
laboratory gamma while reducing the mean longitudinal velocity. Run
`undulator_resonance INPUT.yaml [LAB_GAMMA]` to obtain the standard estimate
`gamma_lab/sqrt(1+K^2/2)`, compare it with the configured boost, and inspect the
predicted box-z drift across the characteristic undulator.

## Optional radiation-resolution preflight

```yaml
radiation_resolution:
  enabled: true
  maximum_photon_energy_eV: 124.0
  warning_grid_points_per_wavelength: 8.0
  warning_maxwell_samples_per_cycle: 8.0
  warning_detector_samples_per_cycle: 8.0
```

This block declares the highest laboratory photon energy that the current run
is intended to interpret. It is independent of how that EUV frequency arose:
the program does not assume an FEL resonance formula, a bunching harmonic, or
a particular magnetic element. The check models a paraxial forward (`+z`)
laboratory mode, transforms its frequency and wavelength into the boosted
frame, and reports:

- longitudinal grid points per boosted wavelength;
- Maxwell/current-deposition samples per boosted optical cycle;
- for every field plane, the conservative realized laboratory sample gap and
  samples per target cycle after quantization to complete Maxwell steps.

Only an unavoidable Nyquist failure is fatal: the grid and Maxwell update must
both exceed two samples per boosted wavelength/cycle, and an enabled field
plane must exceed two samples per laboratory target cycle. Error messages give
the limiting `dz`, `dt`, or detector `rhythm`. The three configurable warning
levels default to eight; falling below them only emits a warning and never
changes the algorithm. Cowan-z still has exact axial vacuum phase velocity for
resolved modes, so the appropriate production margin must be established by a
refinement scan rather than by a universal hard points-per-wavelength rule.

When the block is absent or `enabled: false`, no target frequency is inferred
from the undulator or bunch and no target-band restriction is applied; startup
prints that the check is disabled. No higher-order E/B sampling, extra field
state, communication, or time-loop work is introduced.

The preflight intentionally does not impose `dx/dy < wavelength` for a
paraxial on-axis carrier. Transverse envelope resolution, accepted angle and
aperture, macro-particle noise, CPML reflection, pulse-window length and
off-axis numerical dispersion remain problem-specific convergence checks.

## Field boundary

The boundary policy is explicit and required. Production-oriented no-seed
tests use compact unsplit CFS-CPML:

```yaml
boundary:
  type: cpml
  cells: [8, 8, 8]
  polynomial_order: 3.0
  target_reflection: 1.0e-8
  kappa_max: 8.0
  alpha_fraction: 0.0
```

`cells` is the thickness per face in x, y and z cells. A zero disables CPML
on both faces of that axis without allocating histories; a nonzero value must
be at least two. Opposite layers must leave a non-PML interior. The current z
slab decomposition additionally requires each complete z layer to fit on its
endpoint MPI rank; invalid rank counts fail before particle input.

`target_reflection` is used to grade the maximum conductivity and must be in
`(0,1)`. It is not an assertion that every discrete mode will achieve exactly
that reflection. `kappa_max>=1`, positive `polynomial_order`, and nonnegative
`alpha_fraction` control the CFS profile. The example keeps
`alpha_fraction: 0.0`, which gave the lower reflection for the near-axis
forward pulse relevant to the current FEL target. A nonzero shift is retained
as an explicit experimental option for later low-frequency or evanescent-mode
studies and must be revalidated for that spectrum. Startup logging reports the
total and maximum-per-rank auxiliary memory actually allocated.

For a boundary regression without absorption use:

```yaml
boundary:
  type: pec
```

PEC allocates no boundary histories. CPML currently rejects all nonempty
`incident_waves` lists so the no-seed absorbing boundary can be validated
before the TF/SF surface is redesigned. The numerical construction and Cowan
stability conditions are detailed in
[MAXWELL_COWAN_CPML.md](MAXWELL_COWAN_CPML.md).

## Initial particle self-field

```yaml
initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-10
  maximum_iterations: 10000
```

The block is optional and defaults to the values above. When enabled, the
APRL deposits the boosted bunch with the production CIC shape, sums
shared MPI vertex planes, and solves the discrete relativistic-Poisson
equation before the first physical step. `relativistic-poisson` is the default:
it derives the represented-mass-weighted mean axial velocity in the simulation
frame, uses the Vay rigid-beam elliptic operator, and initializes both Yee E
and the associated B field. `electrostatic-poisson` retains the old B=0
construction only for controlled regression. The resulting Yee-edge E field
is checked against the same charge deposition. Non-convergence and a
post-check residual above the configured tolerance are hard errors;
`maximum_iterations` must be positive and `relative_tolerance` must lie
strictly between zero and one.

With CPML, the zero-potential surface for the static solve is the **inner CPML
surface**, and E/B inside CPML start at zero. This prevents a non-radiative
Coulomb tail from charging zeroed CPML memory variables and then flowing back
into the physical box. Radiation still sees the normal CPML update. The run
stops if any CIC charge weight reaches this static boundary. Move every
particle at least one complete cell inside the CPML entrance or enlarge the
physical mesh padding rather than silently dropping boundary charge.

The scalar potential and four distributed CG work slabs exist only during
initialization and are included in the peak-memory estimate. They are released
before time advance, so the propagation state remains E/B-only. The static
boundary is still a finite-domain approximation: production work must converge
the distance from the bunch to the CPML entrance. The relativistic construction
assumes a rigid common axial velocity; velocity spread, envelope mismatch and
macroparticle noise remain physical/numerical convergence questions.
Disabling the block is intended for controlled regressions and emits a startup
warning. See [INITIAL_SELF_FIELD.md](INITIAL_SELF_FIELD.md) for the discrete
MPI ownership and physical limitations.

## Beam reference and input

Every particle position is relative to a laboratory-frame beam reference:

```yaml
beam:
  reference:
    input_plane_z: -2.0
    initial_center_z: -1.8
```

All particle planes, magnetic elements and detector planes use the same
absolute laboratory z coordinate. No element is forced to lie at zero.
For production HDF5-v3/v4 input, `input_plane_z` is the Elegant observation-plane
coordinate. `initial_center_z` is the coordinate assigned to the reconstructed
common-time bunch reference after straight-line forward projection. In v3,
positive particle longitudinal offset points downstream. In v4, the file keeps
the signed fixed-plane arrival-time offset and startup reconstructs the
longitudinal position with each particle's own velocity. Versions 1/2 and the
generated Gaussian path are already common-time snapshots and ignore
`input_plane_z`.

For HDF5-v3/v4 records the program requires every reconstructed particle position
to be at or downstream of `input_plane_z`. It then requires the snapshot bunch
front, plus one lab-equivalent boosted z cell, to remain before the first
magnetic-element **interaction** entrance. These bounds define an admissible
interval for `initial_center_z`. If it is empty, the error reports both limits
and recommends moving the Elegant plane, changing the magnetic fringe/action
region, reducing the bunch span, or refining `dz`.

The subsequent Lorentz transformation anchors boosted time zero to the
downstream bunch-front event. Other particles are synchronized relative to
that event and the transformed bunch is translated to the centre of the
boosted numerical box. A potentially metre-scale span between the corresponding
laboratory events is logged only as a relativity-of-simultaneity diagnostic; it
is not interpreted as a required physical drift and does not move
`initial_center_z` upstream.

Production input uses HDF5:

```yaml
  input:
    type: hdf5
    file: ../examples/particles/example_particles.h5
    electrons: 4.0
    position_offset: [0.0, 0.0, 0.0]
```

Relative `file` paths are resolved from the YAML file directory. `electrons`
is always the total physical electron count represented by all records.
HDF5 v2/v3/v4 records carry a positive relative `macro_weight`; APRL
normalizes their sum to `electrons` and scales each macro charge and mass
together. Legacy v1 files remain readable and imply equal weights.
`position_offset` is optional and uses the YAML length unit. Its x/y components
shift the fixed-plane coordinates; for v4 its z component shifts the relative
snapshot after arrival-time synchronization. The binary schema is defined in
[PARTICLE_INPUT_HDF5.md](PARTICLE_INPUT_HDF5.md).

For small integration tests only, a deterministic Gaussian can be generated
without an input file:

```yaml
  input:
    type: generated-gaussian
    electrons: 8.0
    macroparticles: 8
    gamma: 4.0
    direction: [0.0, 0.0, 1.0]
    center: [0.0, 0.0, 0.025]
    sigma_position: [0.02, 0.02, 0.01]
    sigma_proper_velocity: [0.0, 0.0, 0.0]
    random_seed: 17
```

`proper_velocity` means the dimensionless vector gamma*v/c. Generation is
counter-based and therefore gives the same global particles for different MPI
rank counts. This path is intentionally limited to one uncorrelated Gaussian
and is not a beam-preparation model.

## Sources

Incident waves are injected through a closed TF/SF surface and then advanced
by Maxwell's equations. Available profiles are `plane`, `truncated-plane`,
`gaussian`, `super-gaussian`, and their `standing-*` forms. Amplitude is either
`peak_electric_field_V_per_m` or `normalized_amplitude`.

```yaml
sources:
  incident_waves:
    - profile: gaussian
      position: [0.0, 0.0, -4.0]
      direction: [0.0, 0.0, 1.0]
      polarization: [1.0, 0.0, 0.0]
      normalized_amplitude: 0.001
      wavelength: 0.8
      radius: [5.0, 5.0]
      order: [0, 0]
      envelope:
        type: gaussian
        center_time: 0.0
        duration: 0.02
        carrier_phase_rad: 0.0
```

Envelope types are `neumann`, `gaussian`, `secant`, `flat-top`, and
`inverse-gaussian`. `rising_cycles` and two-value `inverse_gaussian_sigma` are
available where relevant.

Magnetic devices remain prescribed laboratory-frame fields, transformed only
at particle events; they are not duplicated on the Maxwell grid. They are
optional. With no magnetic device, initialization keeps the same user-defined
beam-centre placement and Lorentz synchronization but omits only the
first-magnetic-interaction clearance test. This supports detector-only and
element-free free-propagation cards without inventing a dummy magnet.

```yaml
  magnetic_elements:
    - type: planar-undulator
      characteristic: true
      strength_parameter: 0.1
      period: 10.0
      periods: 1
      entrance_z: 0.0
      polarization_angle_rad: 0.0
      gaussian_fringe: true
      fringe_relative_cutoff: 1.0e-9
```

`uniform-dipole` instead uses `length` and `field_T`. When magnets exist, the
first entrance is the minimum physical `entrance_z` across them, independent of
YAML ordering. It shares the laboratory coordinate system with
`beam.reference.initial_center_z`; it is not forced to zero.

A planar undulator keeps the legacy divergence-compatible analytical end-field
idea, but fixes its ambiguous infinite support. The Gaussian is multiplied by
a quintic compact-support taper; the transverse component is derived from the
longitudinal envelope so `div(B)=0` remains satisfied through both the entrance
and exit. `fringe_relative_cutoff` defines the raw Gaussian value used to place
the compact edge and must lie in `(0,1)`. With `gaussian_fringe: false`, the
interaction and physical ranges coincide.

`characteristic: true` is optional simulation metadata used by
`undulator_resonance`. When a card has multiple planar undulators, mark
exactly one to select the element used for the preflight resonance and
constant-boost recommendation.

## Runtime strategy

The runtime strategy is global and independent of trajectory or detector
selection. Use the small-server test path when a run must be stoppable at any
time:

```yaml
runtime:
  mode: interactive
  stop_check_interval_steps: 16
```

`interactive` installs minimal SIGINT/SIGTERM handlers. Every configured
number of completed field steps, ranks coordinate one stop flag. A requested
stop never interrupts a Maxwell or particle update: the program leaves the
loop at a complete step, commits and closes every enabled trajectory and
detector file, and marks the files incomplete but readable. Interactive
trajectory output also makes committed prefixes durable every
`flush_every_samples`. This remains available with trajectory output disabled.

For scheduled supercomputer production use:

```yaml
runtime:
  mode: throughput
  resource_monitor:
    enabled: true
    progress_interval_steps: 1000
    calibration_steps: 1
    memory_safety_factor: 1.25
    time_safety_factor: 1.25
```

`throughput` is the default. It installs no signal handlers, performs no
stop-signal MPI polling, and disables periodic trajectory durability flushes;
normal physical-stop shutdown still closes every output. This avoids the
testing path's synchronization and filesystem costs. `local-test` and `hpc`
are accepted aliases for `interactive` and `throughput`, respectively.


`resource_monitor` is independent of the stop policy. When enabled, startup
models the allocated field, CPML, particle, detector, trajectory and MPI halo
buffers and reports maximum per-rank and aggregate memory. With no incident
wave it also benchmarks `calibration_steps` zero-field Maxwell updates and a
bounded sample of particle push/deposition work without advancing the physical
state. The estimate applies the configured factors, which must be at least one.
The current seed-wave combinations are rejected earlier, so skipping the
microbenchmark for a nonempty incident-wave list is an explicit future-facing
safeguard rather than a hidden fallback.

`progress_interval_steps: 0` suppresses periodic records but retains the
startup estimate and final measurement. Otherwise rank zero writes a flushed,
single-line `[resource]` record at that interval. The final record includes
current and peak resident memory plus measured loop and full wall time. This
format is intentionally suitable for stdout/stderr capture by `sbatch`; it
does not require an interactive terminal or a separate monitoring process.
The time estimate excludes filesystem contention, workload changes caused by
particle migration, and earlier physical stopping, so production allocations
should use safety factors calibrated on the target machine and rank layout.

For compatibility with input cards from the preceding commits,
`trajectory.mode` is still accepted as a global runtime-mode alias when
`runtime.mode` is absent. New cards should use `runtime.mode`; specifying both
with different values is an error.

## Runtime energy ledger

The optional global particle/field audit is independent of detector and
trajectory output:

```yaml
energy_ledger:
  enabled: true
  directory: output/example
  filename: energy-ledger.h5
  sample_interval_steps: 10
  buffer_records: 64
  compression: 0
  warning_relative_tolerance: 0.01
```

It records boosted-frame particle kinetic energy, physical-interior E/B
energy, all six inner-CPML Poynting fluxes, removed-particle kinetic energy and
prescribed-device work. Complete volume and particle statistics are evaluated
only at `sample_interval_steps`; boundary power and prescribed work are
accumulated every Maxwell step. Only MPI rank zero writes the small HDF5 file.
Omitting or disabling the block removes all of this work and storage.

The same records contain laboratory gamma mean/rms, a linear longitudinal
chirp and the rms remaining after that chirp is removed. These momentum
statistics share an equal-box-time slice and must not be added to the
boosted-frame energy terms. Use fixed laboratory particle planes for rigorous
accelerator entrance/exit or slice diagnostics. See
[ENERGY_LEDGER.md](ENERGY_LEDGER.md) for the equation, schema, source-validity
flags and current gamma-1000 control result.

## Stop strategy

Exactly one physical stop strategy is selected in addition to `mesh.duration`,
which remains a maximum-runtime guard. Exhausting that guard before the chosen
physical condition is an error and leaves trajectory output marked incomplete.

Stop after every valid particle has passed the last element's interaction
exit:

```yaml
stop:
  mode: after-last-element
```

Here an element may be a magnetic device, a field plane, or a particle plane.
A detector can therefore be the first and only element in a pure-drift run.
If there are no enabled elements at all, `after-last-element` is rejected while
loading the card because it has no finite physical target; use
`reference-center-z` instead. This prevents an element-free job from merely
running until `mesh.duration` and failing late.

With CPML enabled, its inner surface is the physical particle boundary. A
particle trajectory and detector participation end exactly at the first CPML
entrance. The rest of that step continues as an output-free ballistic carrier
whose deposited current is damped by the CPML conductivity profile. Residual
carrier charge is exported through a virtual outward current at the outer box.
The physical particle is immediately excluded from the all-particles predicate;
carriers do not postpone this stop mode. If no valid particles remain, this
mode stops and reports that fact instead of claiming that the bunch crossed
the downstream boundary. Use `reference-center-z` when field ring-down after
the physical particles leave must be retained.

This policy has no YAML switch. A face without CPML uses the direct
charge-conserving outer boundary instead. Production runs allocate no
face-sized particle-flux arrays and do no boundary I/O; shutdown reports
separate CPML-entry, direct-exit, and residual-carrier totals. See
[PARTICLE_OPEN_BOUNDARY.md](PARTICLE_OPEN_BOUNDARY.md) for the damping,
continuity, and MPI details.

Alternatively, stop when the boost reference centre reaches a fixed lab z:

```yaml
stop:
  mode: reference-center-z
  z: 2.5
```

The reference centre starts at `beam.reference.initial_center_z` and follows
the boost-frame origin worldline. The target must lie strictly downstream of
that initial centre and every current element interaction region, allowing
deliberate observation after the last device or a completely element-free
free-propagation run.

Stopping operates on generic beamline-element extents. Magnetic devices and
both laboratory detector-plane types participate in the same first/last
boundary logic. The stop boundary of a field detector is its sampling plane,
not the entrance of its left diagnostic region. Particle planes remain exempt
from element-overlap rules, so they may also precede every magnetic device.

## Experimental particle retirement

The optional retirement layer is an alternative to the field detector's
ballistic charged-particle background reconstruction:

```yaml
particle_retirement:
  enabled: true
  entrance_z: 5.0
  length: 2.0
```

A field detector can bind this global layer to a manual spectral guard:

```yaml
      retirement_frequency_protection:
        enabled: true
        minimum_photon_energy_eV: 50.0
        cycles: 2
```

Both positions use the configured length unit and are fixed laboratory-frame
coordinates. At `entrance_z`, a physical particle stops participating in the
push, ordinary trajectory cadence, and detector crossings. A compact
ballistic carrier continues along its instantaneous velocity and deposits
current with the C2 profile

```text
w(s) = 1 - 10 s^3 + 15 s^4 - 6 s^5,   s in [0,1].
```

It is removed at `entrance_z + length`. This acts only on particle current;
the E/B solver and fields already present on the grid continue unchanged.
There is no right-side diagnostic region.

The method is deliberately experimental. Tapering current reduces the sharp
high-frequency impulse of immediate deletion, but removal of net charge inside
the Maxwell domain is not an exactly charge-continuous operation. Therefore a
matched zero-radiation baseline and the power tests in
[FIELD_POWER_COMPARISON.md](FIELD_POWER_COMPARISON.md) are mandatory. The
baseline must be subtracted at field amplitude, never by subtracting scalar
powers.

The guard never chooses or changes `particle_retirement.length`. After the
particles are read, it uses their maximum laboratory gamma and checks

```text
beta_max = sqrt(1 - 1/gamma_max^2)
delay    = length/c * (1/beta_max - 1)
observed_cycles = minimum_photon_energy_eV * delay / h
observed_cycles >= cycles
```

The photon energy is the lower edge of the radiation band that must be
protected, not necessarily the undulator resonance centre. The same detector
also reserves the conservative left causal interval

```text
extent_x        = mesh.cells.x * mesh.cell_size.x
extent_y        = mesh.cells.y * mesh.cell_size.y
rho_guard       = sqrt(extent_x^2 + extent_y^2)
causal_distance = gamma_max * rho_guard
causal_z         = detector_z - causal_distance
```

The retirement exit plus one lab-equivalent z cell must lie left of
`causal_z`. An invalid length or placement exits before field allocation and
prints a recommended retirement length or detector z. This protected detector
must be the final beamline element; no magnetic, field-detector, or
particle-detector element may lie downstream. Field-detector regions may still
overlap, particle planes remain absent from region-exclusion checks, and the
protected interval excludes magnetic interaction regions in the same way as
the ballistic-reference interval.

Only one field detector can bind the global retirement layer. Its
`particle_background.enabled` must be false. The ballistic-reference setting
remains the independent alternative: it changes no particle or Maxwell state
and therefore has no restriction on elements downstream of its sampling plane.

Initialization enforces all of the following:

- the entrance is at least one lab-equivalent z cell beyond the final magnetic
  interaction (including its configured fringe support);
- every field plane is at least one such cell beyond the taper exit;
- field-plane `particle_background` is disabled because the two removal routes
  are alternatives;
- when frequency protection is enabled, the manual taper length spans the
  requested number of cycles at the minimum protected energy;
- the taper exit is upstream of the detector's derived causal interval and the
  protected detector is the final beamline element;
- at least one field plane exists;
- `stop.mode` is `reference-center-z`, and its z lies beyond the final field
  plane, so retirement does not terminate the run before the field is sampled.

Particle planes remain non-exclusive and may be placed at the retirement
entrance for small diagnostic runs. Other field planes may overlap one another
but must all be downstream of the retirement exit. When retirement is disabled,
it allocates no carrier vector and adds no active carrier work.

## Laboratory detector planes

The only detector mechanisms in this program are fixed laboratory-frame
field planes and particle-crossing planes. The predecessor's boosted-frame
tracking and moving diagnostic surfaces are deliberately not accepted.

```yaml
detectors:
  enabled: true
  directory: output/example/detectors
  field_planes:
    - name: exit-fields
      z: 2.5
      rhythm: 0.002
      buffer_samples: 2
      compression: 0
      particle_background:
        enabled: true
        buffer_records: 16384
        compression: 0
        validation:
          enabled: false
          maximum_particles: 100000
  particle_planes:
    - name: exit-particles
      z: 2.5
      buffer_records: 4096
      compression: 0
```

The two charged-particle-field treatments are detector settings:
`particle_background` records a virtual straight-line reference for offline
analytic subtraction, whereas `retirement_frequency_protection` audits a
manually configured global retirement layer. They cannot be enabled together
on one field plane or in one retirement run. `particle_background` is disabled
by default and must be enabled explicitly. Omitting it creates no companion
file, crossing events, buffers, communication, or detector reference region.

A particle plane is a zero-length beamline element at lab `z`. A field plane's
physical and stop location is also `z`, but when `particle_background.enabled`
is true it owns a diagnostic-only interaction interval immediately to the
left. Either plane can therefore be the last element used by
`after-last-element`, and a `reference-center-z` stop must lie downstream of
its sampling plane. Detector names are unique and may contain letters, digits,
`.`, `-`, and `_`.

After laboratory particles are read, the field detector derives

```text
extent_x          = mesh.cells.x * mesh.cell_size.x
extent_y          = mesh.cells.y * mesh.cell_size.y
rho_guard         = sqrt(extent_x^2 + extent_y^2)
gamma_guard       = maximum initial laboratory particle gamma
reference_length  = gamma_guard * rho_guard
reference_z       = detector_z - reference_length
```

The multiplied mesh extents are full transverse widths, so this is a deliberately
conservative full-diagonal rule. When a physical particle crosses
`reference_z` downstream, the detector records one laboratory position,
proper velocity, charge, mass, and weight. This record defines a virtual
straight line for charged-particle background reconstruction. It never changes
the particle push, current deposition, MPI migration, or Maxwell fields. The
main field file remains the raw total Maxwell field; the companion reference
file is intended for later convolution/subtraction on the analysis machine.

`particle_background.validation` is a test-only two-plane diagnostic. When
enabled, the real particle is sampled once more as it crosses the field plane.
Rank zero pairs that state with the left-boundary reference by particle ID,
stores the actual and straight-line-predicted exit states, and reports
transverse-position, arrival-time, proper-velocity, and direction errors.
`maximum_particles` is mandatory protection against accidentally allocating a
large rank-zero pairing table: startup exits if the global input macroparticle
count exceeds the configured value. Validation requires
`particle_background.enabled: true`.

When validation is false, no pairing table or validation dataset exists and no
extra field-plane crossing event is generated.

Initialization applies the following independent rules:

- magnetic interaction regions may not overlap other magnetic regions;
- a magnetic interaction region may not overlap a field detector's left
  reference region;
- field detector reference regions may overlap one another;
- particle detector planes never participate in exclusion checks.

An invalid magnetic/field-detector placement exits before field allocation and
prints the minimum recommended downstream detector `z`. The left-only region
and its Lorentz-transformed moving planes are detailed in
[FIELD_DETECTOR_REFERENCE.md](FIELD_DETECTOR_REFERENCE.md).

A field plane stores laboratory E and B over all x-y cell centres. `rhythm` is
a laboratory-time minimum interval; actual sample times are stored because
they are quantized to completed Maxwell steps. `buffer_samples` controls the
rank-zero HDF5 batch and should normally stay small because one sample is a
full x-y plane. The stored Yee B value remains half a Maxwell step staggered
from E; that fact is recorded in file metadata rather than hidden by a second
full-domain field copy.

A particle plane stores a record only when a particle segment actually
crosses its fixed lab z, including downstream/upstream direction. It does not
change the field reference records, normal trajectory cadence, or trajectory
records. `buffer_records` controls only the separate particle-detector writer.

Omitting `detectors`, using empty plane lists, or setting `enabled: false`
constructs no detector manager, allocates no detector buffers, opens no
detector files, and enters no detector MPI communication. When enabled, only
the rank owning a field plane samples it; particle crossings are sent in
bounded batches. MPI rank zero alone opens and appends detector HDF5 files, so
ranks never compete for a shared detector file. Detector output uses
throughput buffering and is finalized after trajectory files, keeping the two
I/O paths independent. See [DETECTOR_OUTPUT_HDF5.md](DETECTOR_OUTPUT_HDF5.md).

## Trajectory output

```yaml
trajectory:
  enabled: false
  directory: output/example
  basename: particles
  rhythm: 0.002
  buffer_records: 64
  flush_every_samples: 1
  compression: 0
```

Trajectory output is intended for small-particle validation rather than the
production radiation path. When disabled, APRL opens no trajectory
file, reserves no record buffer, performs no periodic trajectory sampling,
and skips terminal-record conversion. When enabled, each MPI rank writes one
HDF5 file. `rhythm` uses the YAML time unit. In global
interactive runtime mode, `flush_every_samples` periodically commits a
readable prefix. Global throughput mode buffers normally until a full batch or
clean close. `compression` is 0-9; zero minimizes CPU cost. When `enabled` is
false, `rhythm` may be omitted. Format version 2 additionally stores exact
CPML-entry or direct-domain-exit terminal events between cadence samples.
Numerical CPML carriers are never written. See
[TRAJECTORY_OUTPUT_HDF5.md](TRAJECTORY_OUTPUT_HDF5.md).
