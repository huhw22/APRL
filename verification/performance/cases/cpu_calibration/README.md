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

Run `tests/cpu_calibration/run.sh --case balanced --mode pilot` first. After
the pilot passes, use `--mode matrix --repeats 3` on each of balanced,
field-heavy, particle-heavy and holdout. Every run has a 180-second limit and
the entire command holds an exclusive workspace lock. New results are always
placed in unique `results/cpu_calibration` directories. The lock cannot
reserve host resources against other users; record load and inspect spread.

Fit a hardware-specific rule using only the three training cards and validate
the layout selection on the holdout. Do not transfer its coefficients to the
128-core EPYC 7742 nodes or extrapolate to more nodes. If fit or rank choice
fails the holdout, report the measured decision table instead.
