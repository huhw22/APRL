# Particle HDF5 input specification

## Purpose and ownership

This is the production particle-input contract. It is compact, directly
sliceable by MPI ranks, and independent of YAML display units. A file describes
one laboratory-frame bunch snapshot relative to a user-defined bunch-centre
reference. The YAML card supplies the total physical electron count and gives
that reference an initial lab z coordinate; the first physical beamline
entrance defines global z=0.

## Required layout

```text
/particles                         group
  @format_version = 1              signed integer attribute
  /records                         1-D compound dataset, N > 0
```

Every `/particles/records` item contains:

| Field | HDF5 type | Meaning |
|---|---:|---|
| `position_m[3]` | little-endian float64 | relative laboratory x,y,z in metres |
| `proper_velocity[3]` | little-endian float64 | dimensionless gamma*v/c |
| `source_id` | little-endian uint64 | stable source record identifier |

All numeric values must be finite. `source_id` need not be contiguous, but
should be stable across regenerated inputs so that later trajectory analysis
can join records. The provided converter assigns one-based input-row IDs.

Files produced by the converter also carry descriptive attributes
`coordinate_frame=relative-lab`, `position_unit=m`,
`proper_velocity_unit=gamma*v/c`, and `source_id_definition`. Readers rely on
the version and compound field names/types; the descriptive strings are for
inspection.

## MPI read behavior

For N records and P ranks, the simulator assigns one contiguous, nearly equal
hyperslab to each rank and performs a collective MPI-IO dataset read. It does
not make rank zero read and redistribute the whole file. A multi-rank run
therefore requires a parallel-HDF5 build; a serial HDF5 build is accepted only
for a one-rank run and otherwise fails explicitly.

The YAML `beam.input.electrons` value is divided uniformly over N records.
Macro-particle charge and mass therefore scale together, leaving the physical
charge-to-mass ratio unchanged. `position_offset`, if present, is applied after
reading and before placement relative to the first magnet.

## Text conversion input

The compiled `particle_text_to_hdf5` utility accepts:

```text
x  y  z  ux  uy  uz
```

Position columns use the command-line `--length-unit`; proper velocity is
always gamma*v/c. The converter validates all six values, rejects extra
columns, ignores blank lines and trailing `#` comments, and writes in bounded
chunks. See [BUILD_AND_RUN.md](BUILD_AND_RUN.md) for its command line and
`examples/particles/example_particles.txt` for a minimal source file.
