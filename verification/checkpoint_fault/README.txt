Checkpoint fault acceptance (explicit execution only)

This directory provides a test-only HDF5 interposer and a bounded campaign.
It does not modify production source, MPI transport, cluster settings, or
existing results. Importing run.py and preparing a campaign launch no job.

Build fault.c as a shared library with the same HDF5 headers as APRL, into a
fresh build directory. For example (replace the header prefix explicitly):

  mpicc -std=c11 -O2 -fPIC -shared -I /chosen/parallel-hdf5/include \
    fault.c -ldl -o /fresh/build/libcheckpoint_fault.so

The interposer resolves H5Fcreate/H5Dwrite from APRL's already-loaded HDF5;
do not link another HDF5 library into it. Nothing activates unless the exact
APRL_CHECKPOINT_FAULT_ENABLE=1 flag and rank/phase/marker variables are set.

Prepare a NEW result directory:

  python3 run.py prepare --directory /fresh/results/checkpoint-fault \
    --binary /absolute/build/aprl --library /absolute/libcheckpoint_fault.so \
    --ranks 4 --target-rank 1 \
    --launcher-json '["mpirun","--bind-to","core","-np","4"]'

After binding and ordinary correctness checks, explicitly execute inside a
Slurm allocation, or on the separately authorized gnode1 with <=16 ranks:

  python3 run.py run --directory /fresh/results/checkpoint-fault

Launcher and rank count must agree. The script does not submit resources.
Each case uses 8 deterministic generated particles and a small mesh, one
thread/rank, checkpoint at completed step 1, checkpoint timeout 3 seconds.
The normal shutdown/finalize timeouts remain 600/60 seconds.

Cases:
1. create: only the selected rank sleeps 20 seconds before snapshot H5Fcreate.
2. write: only the selected rank sleeps 20 seconds before the first snapshot
   H5Dwrite, after successful collective file creation.
3. bad-parent: checkpoint directory has a regular file as parent; require an
   explicit directory failure, with no published or partial snapshot.

For injected cases, all ranks must independently reach their existing
process-local watchdog and call _exit(124). The test-only _exit interposer
records its status and monotonic time using async-safe write/clock_gettime,
then immediately invokes Linux exit_group. No MPI/HDF5 calls or extra waits
are inserted on the exit path. Per-rank .enter and .exit binary markers are
decoded into RESULT.json. Every rank must have an actual 124 marker; a
launcher returning 124 by itself is insufficient.

Some launchers kill peers immediately after the first nonzero rank exit.
Such a run is conservatively FAIL unless every rank's watchdog marker exists;
do not weaken the check or describe launcher cleanup as watchdog proof. A
launcher-only test option that preserves nonzero peers may be supplied in
launcher-json if appropriate, with that exact option retained in PLAN.json.

GNU timeout independently bounds each launcher to 15 seconds (+3 seconds
kill grace), with a second Python process-group guard. A timeout cannot pass.
Published step-*.h5 must not exist after either injected failure. Partial
files are preserved as failure evidence and are never reused or removed.

The campaign refuses existing result directories, a repeated run, changed
binaries/configurations, and inherited LD_PRELOAD. PLAN.json freezes binary,
library, source and generated card hashes. Each case and the aggregate have
machine-readable RESULT.json. No acceptance is claimed before actual runs.
