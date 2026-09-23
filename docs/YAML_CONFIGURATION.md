# YAML input-card specification

The simulator accepts one ordinary YAML document. Unknown physical defaults
are deliberately avoided: required quantities are reported with their YAML
line when missing or invalid.

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

## Mesh and boost

```yaml
mesh:
  field_solver: cowan-z
  lengths: [40.0, 40.0, 12.0]
  resolution: [1.0, 1.0, 0.2]
  center: [0.0, 0.0, 0.0]
  duration: 0.01
  boost_gamma: 2.0
  particle_steps_per_undulator_period: 32
```

`field_solver` is required. Use `cowan-z` for the z-priority controlled-
dispersion stencil or `yee` for the original regression path. Both operate on
the same staggered E/B lattice; neither allocates A/phi. The Cowan path applies
the transverse smoothing directly while evaluating Faraday's law, rather than
storing three extra smoothed-field volumes.

For `cowan-z`, z must have the smallest mesh spacing: `dx >= dz` and
`dy >= dz`. The program checks this before allocating fields or reading the
particle file, exits on violation, and reports a valid transverse cell-count
and spacing suggestion derived from the configured `dz`. Its time step is
`dt = dz/c`, which makes resolved vacuum propagation exactly dispersion-free
on the z axis. Startup logging reports the squared aspect ratios, all Cowan
stencil coefficients, and transverse-axis phase/group velocity ratios at 16
cells per wavelength. Those transverse figures are diagnostics, not a claim
that an arbitrary oblique mode is dispersion-free.

The present TF/SF seed-wave boundary correction is still the Yee form, so a
configuration that combines `cowan-z` with a nonempty `incident_waves` list is
rejected before particle input. Use `yee` for seed-injection regression until
the generalized Cowan TF/SF correction is implemented.

`lengths`, `resolution`, and `center` describe the boosted computational box.
Each length must be an integral number of cells, every dimension needs at
least three cells, and z needs at least two cells per MPI rank. `duration` is
boosted-frame time. The particle-step setting is reserved for the future
subcycling implementation.

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

## Beam reference and input

Every particle position is relative to a laboratory-frame beam reference:

```yaml
beam:
  reference:
    initial_center_z: -1.8
```

The leftmost physical entrance of the first beamline element defines lab z=0.
`initial_center_z` is the lab coordinate of the particle file's relative z=0
reference at the input snapshot, so it is normally negative. It is also the
initial point of the boost reference-centre worldline used by the fixed-z stop
mode. Positive particle-relative z points downstream.

The simulator performs the free-drift Lorentz simultaneity transform to
boosted time zero and checks the transformed bunch front against the first
element's **interaction** entrance, not merely its physical entrance. If the
front touches that region, the run stops before field allocation and reports a
recommended maximum (more negative) `initial_center_z`. The recommendation
includes one lab-equivalent longitudinal cell as a safety margin.

Production input uses HDF5:

```yaml
  input:
    type: hdf5
    file: ../examples/particles/example_particles.h5
    electrons: 4.0
    position_offset: [0.0, 0.0, 0.0]
```

Relative `file` paths are resolved from the YAML file directory. `electrons`
is the total physical electron count represented by all records; equal macro
charge and mass are assigned to each. `position_offset` is optional, uses the
YAML length unit, and is added before reference placement. The binary schema is
defined in [PARTICLE_INPUT_HDF5.md](PARTICLE_INPUT_HDF5.md).

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
at particle events; they are not duplicated on the Maxwell grid. At least one
element is required because beam placement is defined from the first entrance.

```yaml
  magnetic_elements:
    - type: planar-undulator
      strength_parameter: 0.1
      period: 10.0
      periods: 1
      entrance_z: 0.0
      polarization_angle_rad: 0.0
      gaussian_fringe: true
      fringe_relative_cutoff: 1.0e-9
```

`uniform-dipole` instead uses `length` and `field_T`. The first entrance is the
minimum physical `entrance_z` across all elements, independent of YAML
ordering, and must be zero.

A planar undulator keeps the legacy divergence-compatible analytical end-field
idea, but fixes its ambiguous infinite support. The Gaussian is multiplied by
a quintic compact-support taper; the transverse component is derived from the
longitudinal envelope so `div(B)=0` remains satisfied through both the entrance
and exit. `fringe_relative_cutoff` defines the raw Gaussian value used to place
the compact edge and must lie in `(0,1)`. With `gaussian_fringe: false`, the
interaction and physical ranges coincide.

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
```

`throughput` is the default. It installs no signal handlers, performs no
stop-signal MPI polling, and disables periodic trajectory durability flushes;
normal physical-stop shutdown still closes every output. This avoids the
testing path's synchronization and filesystem costs. `local-test` and `hpc`
are accepted aliases for `interactive` and `throughput`, respectively.

For compatibility with input cards from the preceding commits,
`trajectory.mode` is still accepted as a global runtime-mode alias when
`runtime.mode` is absent. New cards should use `runtime.mode`; specifying both
with different values is an error.

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

Particles that cross a transverse or longitudinal computational boundary are
clipped to that boundary for their final current segment, removed from the
valid set, and excluded from this all-particles predicate. This is an explicit
domain loss, not the old post-undulator soft deletion. If no valid particles
remain, this mode stops and reports that fact instead of claiming that the
bunch crossed the downstream boundary.

Alternatively, stop when the boost reference centre reaches a fixed lab z:

```yaml
stop:
  mode: reference-center-z
  z: 2.5
```

The reference centre starts at `beam.reference.initial_center_z` and follows
the boost-frame origin worldline. The target must lie strictly beyond every
current element interaction region, allowing deliberate observation after the
last device.

Stopping operates on generic beamline-element extents. Magnetic devices and
both laboratory detector-plane types participate in the same first/last
boundary logic.

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
  particle_planes:
    - name: exit-particles
      z: 2.5
      buffer_records: 4096
      compression: 0
```

Every plane is a zero-length beamline element at lab `z`. It can therefore be
the last element used by `after-last-element`, and a `reference-center-z` stop
must lie downstream of it. Detector names are unique and may contain letters,
digits, `.`, `-`, and `_`.

A field plane stores laboratory E and B over all x-y cell centres. `rhythm` is
a laboratory-time minimum interval; actual sample times are stored because
they are quantized to completed Maxwell steps. `buffer_samples` controls the
rank-zero HDF5 batch and should normally stay small because one sample is a
full x-y plane. The stored Yee B value remains half a Maxwell step staggered
from E; that fact is recorded in file metadata rather than hidden by a second
full-domain field copy.

A particle plane stores a record only when a particle segment actually
crosses its fixed lab z, including downstream/upstream direction. It does not
change the normal trajectory cadence or records. `buffer_records` controls
only the separate detector writer.

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
  enabled: true
  directory: output/example
  basename: particles
  rhythm: 0.002
  buffer_records: 64
  flush_every_samples: 1
  compression: 0
```

Each MPI rank writes one HDF5 file. `rhythm` uses the YAML time unit. In global
interactive runtime mode, `flush_every_samples` periodically commits a
readable prefix. Global throughput mode buffers normally until a full batch or
clean close. `compression` is 0-9; zero minimizes CPU cost. When `enabled` is
false, `rhythm` may be omitted.
