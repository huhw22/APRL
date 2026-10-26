Held ELEGANT, independent field-plane output, no retirement — 2026-10-27

User-authorized production: full original41732775-particle input, held startup,
2exclusive128-CPU nodes,256MPIx1,accepted optimized/MPI-sampling-fixed executable.
Field detector is at4.5m. All particle planes, ballistic-reference subtraction,
trajectory output and particle retirement are disabled. The physical stop is
after-last-element, which checks particle state internally and includes the
field location without requiring any paired particle-plane file.

One complete state is written atstep3466, near the prior held run's leading
particle crossing at4.027669m, before the4.5m field plane. This is an approximate
crossing time from historical trajectories, not a new event trigger. It does
not require the entire bunch to leave the undulator. SnapshotI/O timeout3600s;
existing600/60s MPI shutdown protection remains. Future restart needs the same
256-rank decomposition and a new output directory. reset_diagnostics:true permits
new downstream particle planes; crossings before the snapshot cannot be replayed.

prepare_case.py freezes a fresh run directory and never submits jobs. It copies
versioned scripts and saves the exact binary/input hashes. The preflight_result
must contain the passed small paired/field-only held regression for that binary.
run.sbatch ROOT validates hashes, performs256-rank binding, then runsgeometry.py
ROOT as a bounded1CPU Slurm step before allocating the large field arrays. Only
after geometry passes does it launch the fixed absolute executable with256ranks.
There is no node pinning or requeue.144h is the scheduler ceiling, not an estimated
runtime. The default cluster.env executable is not used or changed.

The completion checker requires normal configured dynamic stop and MPI close,
actual manifest grid matching full-input geometry, complete lab-frame fieldH5
with corrected spatial_sampling metadata and every field scalar finite. Field
times cover every completed Maxwell step. Only one detector file is permitted.
The complete finite energy ledger and committed checkpoint block extents/prefix
metadata are also checked. No particle-plane output or exact full-cohort detector
count is required. Boundary losses are reported rather than confused with HDF5
transaction failure. A new whole production-state restart is not performed by
the completion checker; accepted checkpoint regression remains its supporting test.

Logs: aprl.log for production progress; geometry.log for input/layout preflight;
binding.log,execution-environment.json andallocation.txt for launch provenance.
Completion writescompletion-check.json andEXIT_CODE. Fields:detectors/exit-field.h5.
Snapshot:checkpoints/step-000000003466.h5. The physical output remains total native
E/B, not separated radiation; quiet-tail/near-field/time-stagger questions remain.
