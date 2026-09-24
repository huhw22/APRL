# Laboratory trajectory HDF5 output

## File layout

Trajectory output uses one HDF5 file per MPI rank:

```text
<directory>/<basename>-rank-00000.h5
<directory>/<basename>-rank-00001.h5
...
```

There is no shared-file writer and therefore no MPI contention on a trajectory
file. Records are buffered independently by each rank. Interactive runs may
periodically commit a readable prefix; normal completion marks the file
complete, while a coordinated signal stop leaves the committed prefix readable
and marks it incomplete.

The `/trajectory` group attribute `format_version` is `2`. Its
one-dimensional compound dataset `/trajectory/records` contains:

| field | meaning |
| --- | --- |
| `particle_id` | globally stable macroparticle identifier |
| `source_id` | configured particle-source identifier |
| `time_s` | laboratory time in seconds |
| `position_m[3]` | laboratory x, y, z in metres |
| `proper_velocity[3]` | dimensionless laboratory proper velocity `gamma*v/c` |
| `charge_C` | signed macroparticle charge in coulombs |
| `weight` | relative input macro-particle weight (diagnostic metadata) |
| `event_type` | record type listed below |
| `boundary_face` | boundary code listed below, or `-1` |

## Event types

`event_type` has the following stable integer values:

| value | name | meaning |
| ---: | --- | --- |
| 0 | sample | normal cadence sample; `boundary_face = -1` |
| 1 | CPML entry | exact physical-trajectory endpoint at an inner CPML surface |
| 2 | domain exit | exact endpoint at an outer face without an intervening CPML surface |
| 3 | particle-retirement entry | exact endpoint at the experimental retirement entrance; `boundary_face = -1` |

Boundary-face codes are `0=x-`, `1=x+`, `2=y-`, `3=y+`, `4=z-`, and `5=z+`.
The mappings are also stored as HDF5 attributes so a reader need not rely only
on this document.

CPML and retirement numerical current carriers never produce trajectory records. Consequently
the trajectory dataset is suitable for radiation reconstruction without
including paths through the non-physical absorbing layer. An event can occur
between normal sample times and should be retained when reconstructing the end
of a particle history.

## Sampling cadence and radiation time windows

`trajectory.rhythm` is converted from the input card's time unit to seconds
and applied as a cadence threshold in the boosted simulation frame. A sample
is taken at the first completed E/B step at or beyond that threshold; the
particle state is not interpolated back to an exact cadence time. Each sampled
state is then Lorentz-transformed and stored as its own laboratory event.

Consequently, this dataset is an event table rather than a rectangular
`[common_time, particle]` array:

- particles recorded during one boosted-frame step can have different
  `time_s` after the transformation because laboratory time depends on the
  particle's longitudinal position;
- a particle that migrates between MPI slabs continues in another rank file;
  its `particle_id` does not change;
- CPML-entry, retirement-entry, and domain-exit records can occur between periodic samples;
- records are appended in rank-local write order, not globally sorted by
  particle or laboratory time.

The radiation reader therefore reads only the committed prefix of every rank
file, redistributes records by `particle_id`, and sorts each reconstructed
history by `time_s`. It derives segment velocities from consecutive laboratory
positions and times; `proper_velocity` is retained as a diagnostic rather than
used to replace those chords.

For a far-field direction `n`, the relevant window coordinate is not
`time_s` alone but the reduced observer time

```text
u = time_s - n dot position_m/c.
```

For forward relativistic radiation, `time_s` and `z/c` nearly cancel, so the
useful `u` interval can be much shorter than either the laboratory simulation
duration or `trajectory.rhythm`. The far-field tool reports the actual
central-axis range of internal-knot `u` values. Production window choices
should be based on that reported range and on convergence under a finer
trajectory cadence, not inferred from the number of HDF5 rows alone.

## Disabled-output cost

When `trajectory.enabled` is false, no trajectory writer, record buffer, or
file is created. Boundary accounting and carrier evolution remain independent
of this output setting.
