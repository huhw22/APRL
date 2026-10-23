# gnode1 fixed-16-physical-core calibration

This is a short generated-Gaussian CPU benchmark for the single unmanaged
`gnode1` host. The user explicitly authorized bounded CPU execution by SSH on
this host. The launcher rejects every other hostname and requires
`APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1`. Other cluster computation remains under
Slurm. The host has 16 physical cores and 32 SMT threads; use the 16 physical
cores as the fixed budget. CPU binding was checked separately for all four
layouts before these scientific runs.

Cards differ in independent workload size. Within each card, old and new
binaries, all layouts, and repetitions receive identical card bytes. The
Gaussian generator has a fixed seed. Self-field feedback, CPML, initial
Poisson tolerance, physical endpoint and diagnostics are enabled. Every run
must have complete HDF5 outputs, all generated particles at the detector,
zero losses, valid energy closure, and strict bitwise equality against the
same-rank baseline. Profiling is disabled for throughput measurements.

The balanced pilot passed. The subsequent layout matrices for balanced,
field-heavy and particle-heavy passed, but the user then narrowed the task.
Every run has a 180-second limit and the command holds an exclusive workspace
lock. New results are always placed in unique `results/cpu_calibration`
directories. The lock cannot reserve host resources against other users;
record load and inspect spread.

The user narrowed the next measurement to a single old/new version comparison
at identical rank and thread counts. Use `--mode pair --case particle-heavy
--repeats 5` for five paired baseline/candidate runs at 16 MPI ranks and one
thread each. Report separate total and loop wall-time medians and the complete
repeat range. The previous layout and holdout matrix is no longer in scope;
an in-progress holdout run was stopped on request and must not be used as a
completed calibration.
