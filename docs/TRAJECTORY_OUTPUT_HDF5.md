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
| `weight` | statistical macroparticle weight |
| `event_type` | record type listed below |
| `boundary_face` | boundary code listed below, or `-1` |

## Event types

`event_type` has the following stable integer values:

| value | name | meaning |
| ---: | --- | --- |
| 0 | sample | normal cadence sample; `boundary_face = -1` |
| 1 | CPML entry | exact physical-trajectory endpoint at an inner CPML surface |
| 2 | domain exit | exact endpoint at an outer face without an intervening CPML surface |

Boundary-face codes are `0=x-`, `1=x+`, `2=y-`, `3=y+`, `4=z-`, and `5=z+`.
The mappings are also stored as HDF5 attributes so a reader need not rely only
on this document.

CPML numerical current carriers never produce trajectory records. Consequently
the trajectory dataset is suitable for radiation reconstruction without
including paths through the non-physical absorbing layer. An event can occur
between normal sample times and should be retained when reconstructing the end
of a particle history.

## Disabled-output cost

When `trajectory.enabled` is false, no trajectory writer, record buffer, or
file is created. Boundary accounting and carrier evolution remain independent
of this output setting.
