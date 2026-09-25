# Particle HDF5 input specification

## Purpose and ownership

This is the production particle-input contract. It is compact, directly
sliceable by MPI ranks, and independent of YAML display units. Version 3
describes records crossing one fixed Elegant observation plane. The YAML card
places both that plane and the reconstructed bunch reference in one shared,
user-defined laboratory coordinate system. No magnet or detector is required
to define global `z=0`.

## Required layout

```text
/particles                         group
  @format_version = 3              signed integer attribute
  /records                         1-D compound dataset, N > 0
```

Every `/particles/records` item contains:

| Field | HDF5 type | Meaning |
|---|---:|---|
| `position_m[3]` | little-endian float64 | x,y at the fixed input plane; component 2 is the signed longitudinal offset from the reconstructed bunch reference, all in metres |
| `proper_velocity[3]` | little-endian float64 | dimensionless gamma*v/c |
| `source_id` | little-endian uint64 | stable source record identifier |
| `macro_weight` | little-endian float64 | positive relative macro-particle weight |

All floating-point values must be finite and `macro_weight` must be positive.
`source_id` need not be contiguous, but should be stable across regenerated
inputs so that later trajectory analysis can join records. The provided
converter assigns one-based input-row IDs.

For a beam travelling in `+z`, the longitudinal offset is normally constructed
from the Elegant arrival-time difference as `zeta_i=-v_reference*(t_i-t_ref)`,
with the exact sign convention checked so that earlier particles have positive
`zeta`. The simulator forms the target common-time snapshot position

```text
z_target_i = beam.reference.initial_center_z + zeta_i
```

and advances the plane coordinates along the supplied straight-line slope:

```text
x_target_i = x_plane_i + (ux_i/uz_i)*(z_target_i-input_plane_z)
y_target_i = y_plane_i + (uy_i/uz_i)*(z_target_i-input_plane_z).
```

Every particle must have positive `uz`, and every target distance must be
nonnegative. A violation stops before field allocation and reports the minimum
valid centre position. This forward projection is the only Elegant-to-snapshot
transport performed by the simulator.

Versions 1 and 2 remain readable as legacy common-time relative snapshots.
Version 1 has no `macro_weight` field and is interpreted as unit weight for
every record. They do not use `input_plane_z`. New files should use version 3.

Files produced by the converter also carry descriptive attributes
`coordinate_frame=fixed-lab-plane-plus-longitudinal-offset`, `position_unit=m`,
`proper_velocity_unit=gamma*v/c`, `source_id_definition`, and
`macro_weight_definition`. Readers rely on the version and compound field
names/types; the descriptive strings are for inspection.

## MPI read behavior

For N records and P ranks, the simulator assigns one contiguous, nearly equal
hyperslab to each rank and performs a collective MPI-IO dataset read. It does
not make rank zero read and redistribute the whole file. A multi-rank run
therefore requires a parallel-HDF5 build; a serial HDF5 build is accepted only
for a one-rank run and otherwise fails explicitly.

The YAML `beam.input.electrons` value is the total number of physical electrons
represented by the complete dataset. For versions 2 and 3, record `i` represents

```text
electrons_i = beam.input.electrons * macro_weight_i / sum(macro_weight).
```

The simulator evaluates the global weight sum collectively. Macro-particle
charge and mass scale together, leaving the physical charge-to-mass ratio
unchanged. Startup independently sums the represented electron count with
extended precision and exits if normalization or charge-to-mass invariance is
lost. Output `charge_C` is the authoritative represented charge; output
`weight` preserves the relative input weight. `position_offset`, if present,
is applied to the plane x/y and longitudinal offset before forward projection.

## Text conversion input

The compiled `particle_text_to_hdf5` utility accepts:

```text
x_plane  y_plane  zeta  ux  uy  uz  [macro_weight]
```

Position columns use the command-line `--length-unit`; proper velocity is
always gamma*v/c. The optional seventh column is the positive relative weight;
when omitted it defaults to one. The converter validates every value, rejects
extra columns, ignores blank lines and trailing `#` comments, and writes in
bounded chunks. See [BUILD_AND_RUN.md](BUILD_AND_RUN.md) for its command line
and `examples/particles/example_particles.txt` for a minimal source file.
