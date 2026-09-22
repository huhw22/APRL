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
  lengths: [40.0, 40.0, 12.0]
  resolution: [1.0, 1.0, 0.2]
  center: [0.0, 0.0, 0.0]
  duration: 0.01
  boost_gamma: 2.0
  particle_steps_per_undulator_period: 32
```

`lengths`, `resolution`, and `center` describe the boosted computational box.
Each length must be an integral number of cells, every dimension needs at
least three cells, and z needs at least one cell per MPI rank. `duration` is
boosted-frame time. The particle-step setting is reserved for the future
subcycling implementation.

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
  mode: throughput
  buffer_records: 64
  flush_every_samples: 1
  compression: 0
```

Each MPI rank writes one HDF5 file. `rhythm` uses the YAML time unit.
`interactive` mode periodically commits a readable prefix and coordinates
SIGINT/SIGTERM at a complete field step. `throughput` avoids those durability
flushes for scheduled HPC runs. `compression` is 0-9; zero minimizes CPU cost.
When `enabled` is false, `rhythm` may be omitted.
