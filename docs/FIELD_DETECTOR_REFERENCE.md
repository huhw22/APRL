# Field-detector ballistic reference region

## Purpose and scope

A laboratory field plane stores the total propagated Maxwell field. Near a
relativistic electron bunch that field contains both the desired radiation and
the bunch's bound charged-particle field. The detector therefore supports a
compact reference from which a downstream tool can reconstruct the field of a
uniform-velocity continuation and assess or subtract that contribution.

This mechanism is diagnostic only. The physical particle remains in the
self-consistent pusher, deposits its normal current, migrates between MPI
ranks, and responds to fields exactly as before. No particle is frozen or
deleted, and the Maxwell kernel is not modified.

## Left-only geometry

For a field plane at laboratory position `z_D`, startup derives

```text
rho_guard = sqrt(size_x^2 + size_y^2)
gamma_guard = maximum initial laboratory particle gamma
L_reference = gamma_guard * rho_guard
z_F = z_D - L_reference
```

The full transverse mesh diagonal is intentionally conservative. A particle
crossing `z_F` downstream supplies one laboratory-frame state for a virtual
uniform-velocity worldline. The interval `[z_F, z_D]` is called the detector's
reference region, but it is not a new simulation domain and allocates no
volume field arrays.

Only a left interval is needed. The reference is used to predict the charged
particle contribution at the event where the field plane is sampled. Particle
motion after that event cannot change the already sampled event, so a matching
right-side buffer would add cost without defining additional information for
this detector. This remains an approximation to be validated for the chosen
beam, detector aperture, and drift length; it is not an exact mathematical
radiation/bound-field projector.

Both `z_F` and `z_D` are fixed laboratory planes. In the boosted computational
frame they are moving spacetime surfaces. The detector locates crossings by
transforming the endpoints of each completed particle segment to the lab and
interpolating the crossing event there. This avoids treating the two planes as
simultaneous fixed-z surfaces in the boosted frame.

## Placement rules

- Magnetic interaction regions may not overlap one another.
- A magnetic interaction region may not overlap `[z_F, z_D]`.
- Reference regions belonging to different field planes may overlap.
- Particle detector planes do not participate in these exclusion checks.
- Stopping at `after-last-element` uses `z_D`, not `z_F`, for the field plane.

An overlap is rejected during initialization. The error reports a minimum
recommended field-plane position based on the last magnetic interaction exit,
the derived reference length, and one boosted longitudinal-cell margin.

## YAML configuration

```yaml
detectors:
  enabled: true
  directory: output/detectors
  field_planes:
    - name: radiation
      z: 0.4
      rhythm: 0.0002
      buffer_samples: 2
      compression: 0
      particle_background:
        enabled: true
        buffer_records: 16384
        compression: 0
        validation:
          enabled: false
          maximum_particles: 100000
```

The reference file is `<name>-ballistic-reference.h5`. Disabling
`particle_background` creates neither the reference file nor its crossing
events. The total field file `<name>.h5` is unchanged.

## Two-plane validation for small tests

The optional validation samples the same physical particle once at `z_F` and
once at `z_D`, pairs the records by particle ID on rank zero, and compares the
second state with straight-line propagation from the first. For entry proper
velocity `u`,

```text
gamma = sqrt(1 + |u|^2)
v = c u / gamma
Delta t = (z_D - z_F) / v_z
r_pred = r_F + v Delta t
```

The file preserves one validation record per matched particle and summary
attributes for:

- transverse position error;
- signed and absolute arrival-time error;
- relative proper-velocity change;
- change of velocity direction.

This is deliberately not a sampled trajectory. Its communication and storage
scale as two crossing events per particle, not as particles multiplied by
time steps. Nevertheless, the rank-zero pairing table scales with particle
count, so validation is off by default and startup enforces
`maximum_particles` (the macroparticle count, not the represented electron
weight). When disabled, the table, validation dataset, and second crossing
event do not exist.

The comparison measures whether a force-free continuation is accurate over
the chosen reference distance. It does not by itself prove that a later field
subtraction is exact. A useful validation sequence is: check the two-plane
orbit errors, then use full trajectories only for a very small bunch to compare
the existing trajectory-to-far-field result with the expected analytical
spectrum.

## Background-subtraction boundary

The main field file currently remains the raw total Maxwell field. The detector
API also has a nullable laboratory background sampler for future analytical
seed-laser subtraction; it is currently empty and recorded as
`external_background_subtracted=none`. Prescribed undulator fields are not
written into the radiation-oriented field plane and therefore need no such
subtraction.
