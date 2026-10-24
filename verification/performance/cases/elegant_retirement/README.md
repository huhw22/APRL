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

The completion check verifies output transactions, program/mesh identity,
activation and retirement accounting, and direct-path arrival coverage from the
new run's pre-retirement particle records. Radiation interpretation still needs
subsequent review. For this request, confirm normal startup and leave the long
job running until the user asks to inspect it.

For the user's subsequent 0.3 m replacement, compact-geometry.sbatch runs
compact_geometry.py and the complete-input geometry check using one CPU on an
exclusive node. The detector is retirement exit plus the original gamma times
full transverse diagonal clearance plus two lab-equivalent z cells: one required
by the original guard, one numerical margin. The old 7.5 m detector and 50 m
stop floor are not retained.

The former 100 fs quiet tail and 20 fs dynamics allowance are withdrawn. They
forced the moving laboratory detector far beyond the initial bunch's left end,
expanded the left guard from 64 to 18688 cells, and defeated the original compact
left-anchored layout. Keep dz and boost unchanged, use the original 64/128 guard
floor, and align only total Nz to 256. Do not round each guard to 256.

`observed_arrival.py` streams the complete historical 4.03 m exit cohort, uses
its evolved laboratory crossing times and momenta, propagates the carriers to
the 4.33 m retirement exit, and computes the farthest direct path to the physical
detector aperture. It also combines the latest retirement time with the
independent opposite-corner source/detector aperture bound, covering a changed
transverse distribution within the physical aperture. Keep any later projected
out-of-aperture particle in the bound. The endpoint uses the later of this
observed-cohort envelope and the input free-drift prediction, plus one solver/detector sample for discrete
bracketing. No fixed quiet tail or decimal reference-stop padding is added.
The historical startup differs from the new held startup, so this is a
data-backed planning estimate, not a rigorous bound on new self-consistent
trajectories; completion must check the actual new exit-cohort coverage.

`preflight.sbatch` now performs this compact geometry scan and the independent
512-particle held/retirement functionality test on the selected two-node layout.
The small functional test's internal post-retirement check does not set any
production window. Preserve all cancelled runs and prepare one new production
directory only after this preflight passes. There is no K=0 control or additional
old-version timing run.
