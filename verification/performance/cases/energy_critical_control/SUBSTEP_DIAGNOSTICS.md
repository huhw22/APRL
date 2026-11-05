# Converged particle-substep plane diagnostic

Enable `detectors.substep_particle_planes: true` to add an independent
`<particle-plane-name>-substeps.h5` file for every existing particle plane.
The legacy `<name>.h5` file and its whole-Maxwell-segment crossing algorithm
are retained. The option defaults to false; no extra writer, event allocation
or detector collective is entered when disabled.

The added `/particle_substep_plane/records` stores `source_id`, `particle_id`,
`time_s` (lab), `time_box_s`, `position_m`, `proper_velocity` (lab),
`proper_velocity_box`, `gamma_lab`, `gamma_box`, `kinetic_lab_J`,
`kinetic_box_J`, `charge_C`, `mass_kg`, `weight`, and `crossing_direction`.
As in the existing files, macro mass and charge already include the macro
weight. `(gamma, proper_velocity)` is a dimensionless four-momentum: multiply
by macro mass times c squared to obtain `(energy, c*P)`.

Each event uses one actual converged particle substep's endpoint worldline.
Time and position are linearly interpolated to the fixed lab plane inside that
substep. This is finer than the legacy field-step chord; it is not an exact
continuous-time orbit/plane intersection. Endpoint gamma and three-momentum
are interpolated **affinely as a four-vector**. The resulting point can be
off shell. Do not reconstruct gamma from the interpolated three-momentum;
the dedicated scalar gamma and kinetic fields preserve the endpoint energy
quadrature without that convexity bias.

The preparation observer only previews events in per-particle batch-local
storage. Events are committed on the main thread, after the real terminal
time is known. Conservative-source discarded field/particle trials do not
have the detector manager and cannot publish events. Frozen pre-activation
segments use their actual frozen endpoints; activation and MPI transfer
segments retain their original particle treatment. The additional records
cannot feed particle or Maxwell updates.

The new option deliberately rejects checkpoint/restart and filename collisions
with the `-substeps` stem. Existing checkpoint behavior is unchanged when the
option is false.

`substep_checks.cpp` is a software diagnostic fixture, not a physical
Maxwell/particle benchmark. It uses timelike, piecewise linear lab worldlines
and a nontrivial Lorentz boost. The known substep crossing is 2.2 ps while the
legacy whole-segment chord gives 2.0 ps. It checks the scalar affine gamma
against independent endpoint values; frozen and active events; multiple
discarded previews versus one commit; actual terminal-time truncation; MPI
root collection; and byte-identical legacy records with the new option off/on.
Compile against the candidate core library and run on fresh output directories
with 1/2/4 MPI ranks.
