# Laboratory detector HDF5 output

Detector files are independent of the per-rank trajectory files. They use
ordinary serial HDF5 and are opened only by MPI rank zero. No parallel-HDF5
file handle is shared between ranks.

All numeric values use SI and all events are transformed to the laboratory
frame before writing. A clean configured stop sets `complete` to 1. A clean
interrupt or maximum-duration exhaustion leaves it 0. Readers must use only
the prefix named by the corresponding `committed_*` scalar.

## Field plane

Each configured field plane creates `<name>.h5` with group
`/field_plane`:

- `time_s`: float64, shape `[samples]`, actual laboratory sample times;
- `fields`: compound records, shape `[samples, ny, nx]`, with
  `electric_V_per_m[3]` and `magnetic_T[3]`;
- `committed_samples`: uint64 scalar giving the readable sample prefix;
- `complete`: uint8 scalar.

The group attributes include `plane_z_m`, `x_first_m`, `y_first_m`, `dx_m`,
`dy_m`, `nx`, and `ny`. x and y are cell-centre coordinates. E and B are
spatially reconciled from the Yee lattice without allocating a collocated
three-dimensional copy. B is sampled at its stored leapfrog half time; the
`magnetic_time_stagger` attribute states this explicitly.

The plane contains the Maxwell-grid field, including self-consistent fields
and injected Maxwell seed waves. Prescribed magnetic-device fields are not
duplicated into this radiation-oriented output; they remain available from
the input card if a downstream analysis explicitly needs them.

Poynting flux is intentionally not duplicated in the file: downstream tools
can compute `S = E cross B / mu0` from the saved laboratory fields.

## Particle plane

Each configured particle plane creates `<name>.h5` with group
`/particle_plane`:

- `records`: an extensible compound dataset;
- `committed_records`: uint64 scalar giving the readable record prefix;
- `complete`: uint8 scalar.

Each record contains:

- `particle_id`, `source_id` (uint64);
- `time_s` (float64);
- `position_m[3]` (float64, with z exactly equal to the plane position);
- `proper_velocity[3]` (float64, dimensionless `gamma*v/c`);
- `charge_C`, `mass_kg`, `weight` (float64);
- `crossing_direction` (int32: `+1` downstream, `-1` upstream).

The crossing event is linearly located inside the completed particle step in
laboratory spacetime. A particle can legitimately appear more than once if it
later reverses and crosses the plane again.

## I/O and trajectory isolation

Field ownership changes as a fixed lab plane moves through boosted-frame MPI
z slabs. The current owner sends one transient x-y sample to rank zero.
Particle events are aggregated only on steps with crossings and transferred
in bounded point-to-point batches. Rank zero serializes all detector writes;
other ranks never open these files.

Detector writers do not call, flush, or modify the trajectory writer. Their
datasets, buffers, completion markers, and filenames are separate. During
normal shutdown trajectory files are closed first, then detector files, so
the two output paths do not race each other.
