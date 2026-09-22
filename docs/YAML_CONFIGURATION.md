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
    distance_to_first_magnet: 0.45
```

If the first magnetic entrance is at `z_entry`, particle relative z=0 is
placed at `z_entry - distance_to_first_magnet`. Positive particle z points
towards/into the first magnet. The simulator then performs the free-drift
Lorentz simultaneity transform to boosted time zero and checks the transformed
bunch front in the laboratory frame. If it reaches or passes the entrance,
the run stops before field allocation and prints a recommended minimum
distance. The recommendation includes one boosted longitudinal cell expressed
in laboratory length as a safety margin.

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
```

`uniform-dipole` instead uses `length` and `field_T`. The first entrance is the
minimum `entrance_z` across all elements, independent of YAML ordering.

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
