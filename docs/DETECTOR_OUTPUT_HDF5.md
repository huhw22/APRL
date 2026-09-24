# Laboratory detector HDF5 output

Detector files are independent of the per-rank trajectory files. They use
ordinary serial HDF5 and are opened only by MPI rank zero. No parallel-HDF5
file handle is shared between ranks.

All numeric values use SI and all events are transformed to the laboratory
frame before writing. A clean configured stop sets `complete` to 1. A clean
interrupt in interactive runtime mode or maximum-duration exhaustion leaves
it 0. Readers must use only the prefix named by the corresponding
`committed_*` scalar.

## Field plane

Each configured field plane creates `<name>.h5` with group
`/field_plane`:

- `time_s`: float64, shape `[samples]`, actual laboratory sample times;
- `fields`: compound records, shape `[samples, ny, nx]`, with
  `electric_V_per_m[3]` and `magnetic_T[3]`;
- `committed_samples`: uint64 scalar giving the readable sample prefix;
- `complete`: uint8 scalar.

The group attributes include `plane_z_m`, `x_first_m`, `y_first_m`, `dx_m`,
`dy_m`, `nx`, `ny`, `reference_entrance_z_m`, `reference_distance_m`,
`reference_rho_guard_m`, and `reference_gamma_guard`. x and y are cell-centre
coordinates. E and B are
spatially reconciled from the Yee lattice without allocating a collocated
three-dimensional copy. B is sampled at its stored leapfrog half time; the
`magnetic_time_stagger` attribute states this explicitly.

The plane contains the Maxwell-grid field, including self-consistent fields
and injected Maxwell seed waves. Prescribed magnetic-device fields are not
duplicated into this radiation-oriented output; they remain available from
the input card if a downstream analysis explicitly needs them.

The stored field remains the raw total Maxwell field. The attribute
`particle_background_status` states whether a companion ballistic-reference
file was produced, and `particle_background_validation` states whether its
two-plane comparison is active. `external_background_subtracted` is currently `none`; the
detector API already accepts a detector-only laboratory background sampler so
a future analytical injected-laser field can be removed without changing the
propagated grid.

`particle_current_policy` records whether the run used ordinary physical
particles or the experimental C2-quintic retirement current. Retirement files
also carry `particle_retirement_entrance_z_m` and
`particle_retirement_exit_z_m`. These attributes are audit metadata: the field
file remains raw, and no baseline is silently subtracted while writing it.


When the field plane enables the detector-bound frequency guard,
`retirement_frequency_protection` records that initialization passed. The
file also stores:

- `retirement_protected_minimum_photon_energy_eV`;
- `retirement_protected_minimum_cycles`;
- `retirement_required_length_m`;
- `retirement_observed_cycles_at_minimum_energy`;
- `retirement_causal_guard_distance_m`;
- `retirement_causal_guard_entrance_z_m`.

These make the manual length choice and actual maximum-particle-gamma check
auditable without reopening the input particle file.

Poynting flux is intentionally not duplicated in the file: downstream tools
can compute `S = E cross B / mu0` from the saved laboratory fields.
The committed `field_reconstruction` post-processor performs this calculation
for both the raw and particle-background-subtracted fields and stores their
forward power and time-integrated energy. See
[FIELD_RECONSTRUCTION.md](FIELD_RECONSTRUCTION.md).

## Ballistic particle-background reference

For a field plane named `<name>`, enabling `particle_background` creates
`<name>-ballistic-reference.h5` with group `/ballistic_reference`:

- `records`: one compound record per downstream crossing of the automatically
  derived left reference boundary;
- `committed_records`: uint64 readable-prefix marker;
- `complete`: uint8 completion marker.

The record members are `particle_id`, `source_id`, `time_s`, `position_m[3]`,
`proper_velocity[3]`, `charge_C`, `mass_kg`, `weight`, and
`crossing_direction`. The format deliberately matches the particle-plane
crossing payload, but its role is different: it initializes a virtual uniform
velocity worldline and does not freeze or replace the physical particle.

Attributes include `plane_z_m` for the reference entrance,
`field_detector_z_m`, the guard geometry, `role`, and explicit
`affects_particle_push=false` and `affects_maxwell_current=false` contracts.
Only rank zero opens this file. Reference crossings share the detector's
bounded MPI crossing batches and are written once per particle rather than at
trajectory cadence.

### Optional two-plane validation

When `particle_background.validation.enabled` is true, the same group also
contains:

- `validation_records`: one compound record for every reference/field-plane
  pair matched by `particle_id`;
- `committed_validation_records`: uint64 readable-prefix marker.

Each validation record contains the actual field-plane time, position, and
proper velocity; the straight-line-predicted time and position; transverse
position error, signed arrival-time error, relative proper-velocity change,
and direction error. The entry state remains in `records`, so the two datasets
can be joined by `particle_id` without duplicating it.

On a clean or coordinated interrupted close, group attributes report matched,
unmatched-entry, unmatched-exit, duplicate-entry, and invalid-prediction
counts, together with RMS and maximum error measures. Partial runs can have
unmatched entries simply because those particles have not reached the field
plane; readers must still honor both committed-prefix scalars. The validation
path is intended for small test bunches and is guarded by the configured
`maximum_particles` limit.

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
