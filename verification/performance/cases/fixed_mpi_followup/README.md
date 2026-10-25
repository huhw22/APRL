# Fixed-MPI follow-up implementation checks

Each candidate is independent and compared with the already accepted Cowan
build. Run a separate campaign for row-local Cowan smoothing reuse, particle
segment intersection reuse, and quantity-specific diagnostics. Keep separate
builds, required-test results and campaign directories; never use one candidate
as another candidate's control or add their speedups together.

This experiment keeps **4 MPI processes × 4 OpenMP threads** on gnode1's 16 physical cores for every measured simulation. Each rank spans two NUMA domains. It compares the currently accepted optimized executable with a separate trial build; it never changes rank counts, precision, cards, production inputs, or the production executable.

The three cards (`grid016`, `grid128`, `particles262144`) are copied byte-for-byte from the completed memory-scaling campaign. Every run keeps that case's particle weights, mesh, feedback, boundary, initial field tolerance, substeps, physical endpoint, diagnostic cadence, and output cadence. Different cases intentionally differ in work and are not implementation speedups.

Preparation is metadata-only and can run on the login node. Supply the full revision compiled into the trial executable:

```bash
/public/home/wfollowm/micromamba/envs/physics/bin/python3 \
  /public/home/wfollowm/APRL_Workspace/source/aprl/verification/performance/cases/fixed_mpi_followup/run.py \
  --stage prepare \
  --root /public/home/wfollowm/APRL_Workspace/results/fixed_mpi_followup/NEW_DIRECTORY \
  --trial /public/home/wfollowm/APRL_Workspace/build/TRIAL/aprl \
  --trial-revision FULL_40_CHARACTER_COMPILED_REVISION \
  --required-pass /public/home/wfollowm/APRL_Workspace/results/REQUIRED_REGRESSION/PASS
```

The trial must first pass the required correctness suite. `--required-pass` names that completed suite's PASS evidence. Its directory must also contain `EXIT_CODE` equal to zero, `build.sha256` with exactly one matching absolute-path entry each for the current trial executable and its CMake cache, a matching saved `CMakeCache.txt`, and current `driver.log` plus `Testing/Temporary/LastTest.log` proving all 72 required tests passed (including duration-close and duration-barrier-watchdog). A `FAIL` marker rejects the evidence. Preparation freezes all six evidence files, and pilot entry rechecks the successful exit and exact executable identity. Keep the full required-test logs in the same directory. A stale `LastTestsFailed.log` alone is not a failure: current complete matching driver and detailed logs are authoritative. Use `node_survey/run_required.sh` with a fresh independent build/result directory and its required explicit gnode1 MPI mapping.

Control defaults to `build/cowan-accepted-r2-20261025/aprl`, compiled revision `77f1ef95267d037dcd629584e42aedc1b04668f3`, and must match SHA256 `dc537da20eacc15648cc59360a1f61bdc9ef1176950ac0cdfda21714e2ba58bb`. Source HEAD 856e172 is not the compiled revision; the production default in cluster.env is not this control. Both builds must have identical recorded compiler, Release flags, MPI, parallel HDF5, and OpenMP selections in their CMake caches. Cards, scripts, probe and executables are SHA256-frozen. Each actual manifest must confirm its own compiled revision, mesh and four ranks, and each runtime must report four threads.

Run only on the explicitly authorized gnode1, using the shared calibration lock and without fabricated Slurm variables:

```bash
APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 bash \
  /public/home/wfollowm/APRL_Workspace/results/fixed_mpi_followup/NEW_DIRECTORY/scripts/run.sh \
  --root /public/home/wfollowm/APRL_Workspace/results/fixed_mpi_followup/NEW_DIRECTORY --stage all
```

At campaign start, a two-second CPU-idle sample must show at least 95% aggregate idle and no other APRL/MPI process; a residual load average from just-completed tests alone does not block an idle host. The script does not wait for resources or change host policy.

The host/topology/binding checks precede six control/trial pilots, covering all three cases. Every trial and subsequent control is compared **bitwise** to its case's first completed 4×4 control; a mismatch stops the campaign. Checks include complete committed diagnostics, finite values, cohort and particle conservation, physical endpoint, and reported threads. No cross-rank difference is accepted or hidden by relaxed tolerances.

After pilots: six warmups, five paired throughput repetitions per case (30 runs), then six separate profiling runs. Case ordering is deterministically shuffled and the control/trial order alternates. The primary timing is `(elapsed(step 80) − elapsed(step 20))/60`; initialization, the first 20 steps and final close are excluded. The physical endpoint remains unchanged and all runs within each case must complete the same step count. Throughput has profiling off and no process-page sampling; profiling captures one NUMA page snapshot after step 20. Every launch is limited to 180 seconds and the campaign to two hours, including separately invoked stages. No task is launched by `--stage prepare` or `--self-test`.

Outputs: `driver.log` (if the launcher redirects there), `report.json`, `summary.csv`, `comparison.csv`, `profile-summary.csv`, and per-run logs, comparisons and NUMA snapshots. `PILOT_PASS` means only pilots passed; `PASS` means the entire sequence and all strict comparisons completed. Failed or partially executed directories cannot silently restart duplicate runs. Use a new result directory after reviewing a failure.

Interpretation limits: five repeats describe observed timing variability, not a production speedup guarantee; NUMA maps include runtime/temporary pages and do not measure remote accesses or memory bandwidth. MPI timings include waiting and are nested within parent stages. Do not add parent/child times or maxima from different ranks. No old historical executable, production process, GPU, system policy, or MPI network setting is changed.

The original three cards keep their 600-second local shutdown watchdog and
60-second finalize guard. They must reach their normal physical endpoints with
exit zero and complete committed outputs. The duration-exhaustion barrier is
exercised by the required suite, not by changing these throughput cards. A
180-second launch timeout, watchdog exit, incomplete output, or unexpected
nonzero process exit remains a failure; no shutdown text or status is ignored.

The source template's helper fallback is its sibling `../memory_scaling/run.py`.
Preparation freezes that helper as `scripts/memory_driver.py`; only its shared
execution, binding, inspection and comparison functions are reused. Its old
standalone preparation constants are not used by this follow-up driver.

Both variants explicitly use the two read-only FFTW libraries in
`build/cowan-accepted-r2-20261025/runtime-libs`: `libfftw3.so.3` and
`libfftw3_threads.so.3`. Preparation checks their accepted hashes and freezes
both files in the plan; run.sh adds that exact directory to LD_LIBRARY_PATH.
Missing or changed libraries fail before any simulation launch. No system
library installation or search-policy change is made.
