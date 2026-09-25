# Particle HDF5 input specification

## Purpose and ownership

This is the production particle-input contract. It is compact, directly
sliceable by MPI ranks, and independent of YAML display units. Version 4
describes particle events crossing one fixed Elegant observation plane and
preserves their arrival-time coordinate. The YAML card places both that plane
and the reconstructed bunch reference in one shared, user-defined laboratory
coordinate system. No magnet or detector is required to define global `z=0`.

## Required layout

```text
/particles                         group
  @format_version = 4              signed integer attribute
  /records                         1-D compound dataset, N > 0
```

Every `/particles/records` item contains:

| Field | HDF5 type | Meaning |
|---|---:|---|
| `plane_position_m[2]` | little-endian float64 | x,y at the fixed input plane, in metres |
| `arrival_time_offset_s` | little-endian float64 | plane-crossing time minus the file's weighted reference arrival time, in seconds |
| `proper_velocity[3]` | little-endian float64 | dimensionless gamma*v/c |
| `source_id` | little-endian uint64 | stable source record identifier |
| `macro_weight` | little-endian float64 | positive relative macro-particle weight |

All floating-point values must be finite and `macro_weight` must be positive.
`source_id` need not be contiguous, but should be stable across regenerated
inputs so that later trajectory analysis can join records. The provided
text converter assigns one-based input-row IDs; the SDDS converter preserves
`particleID` when present and otherwise uses the selected-page row number.

The converter obtains `proper_velocity` from Elegant's standard
`(xp,yp,p)` columns without a paraxial approximation:

```text
uz = p/sqrt(1+xp^2+yp^2)
ux = xp*uz
uy = yp*uz.
```

Here `p=beta*gamma` is the momentum magnitude. At startup the simulator uses
each record's own velocity and crossing time to construct one common-lab-time
snapshot. The common time is chosen so that the macro-weighted longitudinal
centroid lies at `beam.reference.initial_center_z`; correlations between
arrival time and longitudinal velocity are retained. In symbols,

```text
T = [(initial_center_z-input_plane_z)/c + <beta_z*tau>_w] / <beta_z>_w
r_i(T) = r_i(tau_i) + beta_i*c*(T-tau_i),
```

where `tau=arrival_time_offset_s`. Every particle must have positive `uz`, and
`T-tau_i` must be nonnegative. A violation stops before field allocation and
reports the minimum valid centre position. This ballistic synchronization is
the only Elegant-to-snapshot transport performed by the simulator; it does not
retrack the upstream Elegant lattice.

Version 3 remains readable for the former
`position_m=(x_plane,y_plane,zeta)` representation. It uses the reference-speed
approximation `zeta=-v_reference*(t-t_ref)` and advances x/y with `ux/uz` and
`uy/uz`. Versions 1 and 2 are legacy common-time relative snapshots; version 1
has no `macro_weight` field and is interpreted as unit weight for every record.
Versions 1/2 do not use `input_plane_z`. New Elegant conversions should use
version 4.

Files produced by the SDDS converter also carry the selected SDDS page,
protocol version, source binary/ASCII mode, absolute reference arrival time,
and descriptions of the source-ID and weight choices. Readers rely on the
version and compound field names/types; the descriptive attributes are for
inspection.

## MPI read behavior

For N records and P ranks, the simulator assigns one contiguous, nearly equal
hyperslab to each rank and performs a collective MPI-IO dataset read. It does
not make rank zero read and redistribute the whole file. A multi-rank run
therefore requires a parallel-HDF5 build; a serial HDF5 build is accepted only
for a one-rank run and otherwise fails explicitly.

The YAML `beam.input.electrons` value is the total number of physical electrons
represented by the complete dataset. For versions 2 through 4, record `i` represents

```text
electrons_i = beam.input.electrons * macro_weight_i / sum(macro_weight).
```

The simulator evaluates the global weight sum collectively. Macro-particle
charge and mass scale together, leaving the physical charge-to-mass ratio
unchanged. Startup independently sums the represented electron count with
extended precision and exits if normalization or charge-to-mass invariance is
lost. Output `charge_C` is the authoritative represented charge; output
`weight` preserves the relative input weight. For version 4, transverse
`position_offset` is applied at the Elegant plane and its z component is
applied to the reconstructed relative snapshot after synchronization. This
keeps seconds and metres from being mixed in the file reader.

## Direct Elegant SDDS conversion

The `elegant_sdds_to_hdf5` C++ utility links the official SDDS C library and
reads Elegant's native SDDS file directly. It does not create or parse an
intermediate text dump. The required columns are
`x[m], xp, y[m], yp, t[s], p[m$be$nc]`; the standard optional `particleID`
column is preserved exactly, including SDDS5 `long64`/`ulong64` data. If it is
absent, IDs are the one-based rows of the selected page. Elegant normally uses
equal-charge macroparticles, so weights default to one; a positive numeric
per-row weight column may be selected explicitly.

SDDS files can contain several pages (for example, passes, steps, or bunches).
The converter deliberately writes exactly one page; page 1 is the default and
`--page` makes another choice explicit. Pages are never concatenated silently.
`SDDS5` denotes protocol version 5, not the storage mode: both
`&data mode=binary` and `mode=ascii` are valid. The official library handles
protocol version, row/column-major order, endianness, compression, and page
layout.

```bash
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --page 1 \
  --output particles.h5

# Only when the source actually has unequal per-row weights:
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --output particles.h5 \
  --weight-column macroparticleWeight
```

The Elegant page-level `Charge` parameter is not a per-particle weight. Keep
the authoritative total physical population in YAML
`beam.input.electrons`; use `--weight-column` only for a true per-row relative
weight column.

## Text conversion input

The legacy `particle_text_to_hdf5` utility still accepts:

```text
x_plane  y_plane  zeta  ux  uy  uz  [macro_weight]
```

Position columns use the command-line `--length-unit`; proper velocity is
always gamma*v/c. The optional seventh column is the positive relative weight;
when omitted it defaults to one. It writes format version 3 for backward
compatibility. The converter validates every value, rejects
extra columns, ignores blank lines and trailing `#` comments, and writes in
bounded chunks. See [BUILD_AND_RUN.md](BUILD_AND_RUN.md) for its command line
and `examples/particles/example_particles.txt` for a minimal source file.
