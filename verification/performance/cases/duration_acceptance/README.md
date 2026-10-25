# Cowan and bounded shutdown acceptance

Use a fresh result directory and an independent Release build with the same
MPI, parallel HDF5, compiler and OpenMP settings as the accepted control.
Configure CTest's MPI launcher to this directory's `mpi-launch.sh`. Submit
`run.sbatch RESULT_DIRECTORY BUILD_DIRECTORY` only after preparing cards with
`cases.py prepare` and freezing the build, launchers and scripts in
`build-frozen.sha256`. The job reserves two exclusive 128-CPU nodes for at most
30 minutes, without a node-name constraint. Small required checks use one
node; scientific comparisons use 256 MPI ranks and one thread per rank.

The gates are sequential:

1. Validate binding and frozen inputs. Reproduce the pre-fix output-close race
   with a test-only three-second delay in the root HDF5 close. This is a
   targeted regression reproduction, not an old-version performance run.
2. Run every required check, including delayed close and missing-peer shutdown.
3. Run the missing-peer test on two nodes, one rank per node. Root deliberately
   withholds the new duration barrier; a waiting rank must exit through the
   existing local watchdog. The test uses a five-second watchdog and a longer
   independent driver deadline. A driver timeout is a failure.
4. Run accepted-control and fixed normal physical-stop cases, then two fixed
   duration-limited cases with the delayed HDF5 close. These use the same small
   held/retirement configuration, 512 particles and a 16 x 16 x 4096 grid.
5. Require readable, committed output including empty detector planes and the
   energy ledger. Compare normal control/fixed outputs bitwise, both bounded
   repeats bitwise, and bounded physics against the exact normal-run prefix.
   Bounded runs still return 1 and mark output incomplete.

The fault injection is confined to test-only shared libraries and does not
simulate a powered-off host or failed filesystem. It verifies the actual
duration barrier's local timeout path. Production default watchdog limits
remain 600 seconds for general teardown and 60 seconds for MPI finalization;
rank-local unexpected failures retain their original abort path. No per-step
calculation, precision, grid, MPI layout or production executable is changed
by the shutdown fix.

`ACCEPTANCE_PASS` requires every gate. `required/ctest.log`,
`two-node-fault/run.log`, per-case `aprl.log`, and `acceptance-result.json`
provide the evidence. This campaign checks shutdown and scientific
equivalence; it makes no new performance, radiation convergence or full-scale
output-validation claim. The previously truncated full-scale files remain
failed evidence and must not be relabeled as accepted.
