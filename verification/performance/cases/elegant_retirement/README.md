# Single ELEGANT retirement run

This case uses the accepted optimized executable and the complete migrated
ELEGANT input. It enables the existing C2 particle-current retirement layer,
continues Maxwell propagation, and records a downstream fixed laboratory plane.
The user explicitly selected one production run, with no K=0 transient control.
Treating this downstream field as radiation is the user's working assumption,
not a validation result established by this run.

Start from the historical 160×160 card, not the newer 104×104 performance
template. Keep the beam, total represented charge, transverse cells and spacing,
boost, longitudinal spacing, undulator, particle substeps and diagnostic cadence.
The simulation core is unchanged.

Retirement is incompatible with the historical transverse ballistic-injection
mode. Simply disabling that mode would discard many particles already outside
the boosted initial aperture. This case instead uses the existing held frozen
input-plane startup, with activation at the HDF5 input plane; all input x/y
coordinates must fit the physical aperture. The initial CPML-source inclusion
option is disabled with the incompatible injection mode. This is an explicit
additional startup change, not a claim that only the terminal layer changed.

The field detector has no time-gating option. Its usable sampling window depends
on both reference-center stopping and the plane's passage through the boosted
longitudinal box. A small gap between detector z and stop z is not sufficient.
Run geometry_preflight.py on the complete input on a compute node and retain
its derived stop, lower guard, expected time window and mesh estimate.

Run the small held-input/retirement/MPI preflight first. Then use:

    prepare_case.py --workspace /absolute/workspace \
      --geometry /absolute/preflight/geometry.json \
      --preflight-result /absolute/preflight \
      --output-dir /absolute/new/production-directory

Submit the frozen scripts/run.sbatch with the new directory as its argument
and explicit Slurm output/error paths. The allocation is two exclusive 128-CPU
nodes, 256 MPI ranks, one thread per rank, and a 144-hour limit. The script checks
the input SHA256, frozen configuration and binary, and binding before launching.
It does not submit any control, parameter scan or postprocessing job.

completion-check.json only certifies a configured stop and complete output
transactions. Particle accounting, detector-window coverage and radiation
interpretation still require subsequent review. For this request, confirm normal
startup and leave the long job running until the user asks to inspect it.
