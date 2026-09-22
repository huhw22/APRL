# YAML configuration

The program accepts one YAML document. Underscore-separated keys are used
consistently; the parser reports missing keys and YAML line numbers instead of
silently inventing physical parameters.

## Units

`units.length` accepts `m`, `mm`, `um`/`micrometer`, and `nm`.
`units.time` accepts `s`, `ms`, `us`, `ns`, `ps`, `fs`, and `as`. All numeric
length and time values are converted once at input. The simulation core stores
only metres and seconds. Magnetic field values marked `_T` are tesla and
electric field values marked `_V_per_m` are V/m.

## Mesh and boosted frame

```yaml
mesh:
  lengths: [40.0, 40.0, 12.0]
  resolution: [1.0, 1.0, 0.2]
  center: [0.0, 0.0, 0.0]
  duration: 0.01
  boost_gamma: 2.0
  particle_steps_per_undulator_period: 32
```

The mesh values describe the computational box directly. There is no hidden
longitudinal rescaling. Cell counts must be integral, every axis needs at least
three cells, and z needs at least one cell per MPI rank. `duration` is boosted-
frame time. The particle-steps setting is retained in the configuration but is
not enforced until particle subcycling is implemented.

## Beam

Placement is mandatory. `absolute-lab` preserves input z coordinates.
`head-to-first-element` translates the complete MPI-distributed bunch so that
the global maximum particle z is exactly `distance` before the first magnetic
element entrance.

The `distributions` list supports `ellipsoid`, `file`, `manual`, and
`3d-crystal`. Every distribution specifies electron count, macro-particle
count, gamma, direction, and one or more laboratory-frame positions. Proper
velocity is stored as `gamma*v/c`. File input contains six whitespace-separated
columns per record: x, y, z, ux, uy, uz; positions use `units.length`.

## Sources

`incident_waves` are injected through a closed TF/SF surface and subsequently
advanced by Maxwell's equations. Profiles are `plane`, `truncated-plane`,
`gaussian`, `super-gaussian`, and their `standing-*` forms. An amplitude can be
given either as `normalized_amplitude` or as
`peak_electric_field_V_per_m`.

Example wave:

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

`magnetic_elements` are laboratory-frame device fields used only by the
particle pusher. Supported types are `planar-undulator` and `uniform-dipole`.
This keeps incident Maxwell fields and externally maintained device fields
non-overlapping.

## Trajectories

Each MPI rank owns one HDF5 file. `interactive` mode periodically commits a
readable prefix and coordinates SIGINT/SIGTERM at a complete field step.
`throughput` avoids that synchronization and is intended for scheduled HPC
runs. Stored records are in the laboratory frame and retain stable particle
IDs, source IDs, event time, position, proper velocity, charge, and weight.
