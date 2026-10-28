# Pending admission optimization validation

This campaign compares the frozen R3 executable with an independently built candidate. It does not submit jobs or change solver/MPI system configuration.

Prepare with the physics Python environment:

```bash
python3 prepare.py --workspace /public/home/wfollowm/APRL_Workspace --root /absolute/new/campaign --candidate /absolute/candidate/aprl
```

Each executable must retain its own `runtime-libs/` directory. Preparation records SHA256 for executables, libraries, binding probe, cards, scripts and the small input. The original41732775-record production input remains read-only; its full hash is checked on the allocated node before the full stage. Kernel source identity/build acceptance remains the task owner's responsibility.

The coordinator may submit the frozen `scripts/run.sbatch campaign science`, inspect `SCIENCE_PASS.json`, then submit `scripts/run.sbatch campaign full`. The full driver requires the same frozen campaign and candidate hash to have passed science. There is no automatic submission or resumption.

Science comprises eight short runs: dense31675-source undulator/free-drift ×1MPI×8threads and16MPI×8threads ×baseline/candidate. Both variants retain the original moving field geometry, reference solve, release state, solver, deposition, substeps and diagnostics. The three particle planes require the complete cohort and normal physical termination. Existing field2e-12 and particle128eps-per-component roundoff gates are frozen through a copied review helper; the candidate also undergoes cross-layout comparison. This is correctness coverage, not a layout performance scan.

Full comprises the original undulator4.5m production card with53504 scalar-reference cells,4099 active field cells and all41732775 input records. First baseline/candidate each run20 steps with interval10 profiling, then separate baseline/candidate runs100 steps with profiling disabled. Only output paths and the bounded duration differ from production. Duration is one ULP below the smaller of `N*dt` and N repeated double additions, so N full updates execute without a spurious N+1 update. The original `after-last-element` stop remains.

The R3 configuration rejects an early reference-center stop upstream of its4.5m detector. Therefore the full stage deliberately uses the existing bounded-duration diagnostic contract: exit1 with the specific duration-exhausted reason, exact N steps, closed readable HDF5, committed extents and `complete=false`. Its successful diagnostic status is **TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED**; it is never called a completed physical simulation. The existing MPI-free600/60s shutdown watchdog is preserved in both cards.

The duration checker additionally requires the ordered root trajectory/detector/ledger closes, final boundary/resource reports, and `duration-exhaustion-readiness-sync-complete`. The collective is placed after every rank's local output closes; its root completion proves all ranks entered it. Root-only tracing is normal, so sixteen trace lines are not required. The following expected exception calls `MPI_Abort(...,1)`. Consequently, `srun --kill-on-bad-exit` can print a generic **STEP CANCELLED** and peer **Killed/exit1** during cleanup. These exact banners are accepted only after the recorded readiness completion, expected duration exception and MPI abort with errorcode1, with launcher exit1 before its timeout. Whole-job cancellation, cancellation before that proof, explained scheduler cancellation, timeout/signal/crash/watchdog exits, unexpected rank errors and truncated/nonfinite output still fail. HDF5 transaction checks remain mandatory and never change `complete=false` to true.

Both stages reserve one exclusive128CPU node; the science1MPI case intentionally uses only8 cores for correctness coverage. All compared full runs use16MPI×8threads. Both variants set `OMPI_MCA_btl=self,vader,tcp` and the prior openib warning variable, plus close/cores thread binding. Every used layout is checked with the binding probe. No node is pinned.

Slurm has a40min hard limit; the driver has an active1800s stage budget and each full run a600s limit. A full run starts only if its entire600s allowance still fits. Budget exhaustion preserves `FAILED_OR_INCOMPLETE.json` and stops; it never extends the job or converts partial work into PASS. Small runs have240s per-case limits. Profiling and throughput are reported separately, as are reference preparation and stepping. The throughput result is one early100-step pair (including50–100-step interval), not repeated timing statistics or a full-run speedup forecast.

Scripts do not remove reference/data files and do not modify default binaries, live production cards or cluster configuration.

Energy-ledger validation is declared in `THRESHOLDS.json` before execution and applies to both science and full pairs. Integer counts and clocks are exact; `active_macroparticles` includes waiting particles, and active+removed must equal the original cohort at every sample. Source sums use `128*eps64 + 2*gamma_N(u_longdouble)` on noncancelling source scales (N fixed at31675 or41732775). Residual absolute bounds propagate the energy subtractands; ratios propagate numerator/denominator uncertainty, with an explicitly marked cross-product gate when the denominator is unresolved. Sigma-gamma is checked through central variance, mean uncertainty and the existing zero-resolution clamp. Axial chirp and uncorrelated sigma are reported without equality acceptance because pending z moments use an affine diagnostic trajectory. Their values never control admission or particle evolution; actual admission uses exact old repeated additions. E/B and per-ID x/u/t gates remain unchanged. Reports include error/source-scale and, where resolved, error/particle-energy-change metrics.

An already accepted science stage can be reused explicitly when preparing a fresh full-only campaign:

```bash
python3 prepare.py --workspace /public/home/wfollowm/APRL_Workspace \
  --root /absolute/new/full-campaign --candidate /absolute/candidate/aprl \
  --reuse-science-from /absolute/accepted/pair-campaign
```

This verifies the original eight cases and six comparisons, both exact executable hashes, small-input hash, numerical thresholds and every frozen prior artifact. Fresh full cards must match the prior cards beyond output/reference paths. The new plan records `science_reuse`; it does not fabricate a new `SCIENCE_PASS.json` or alter the original campaign. Submit only the new frozen `scripts/run.sbatch new-campaign full`. The four full runs execute sequentially on that allocation without a node pin; this keeps the new baseline/candidate timing pair on the same host. Earlier results remain independent evidence. The driver never automatically reuses a timing sample from another allocation or resumes a partial run.

The candidate's actual two-node replay is prepared separately after the one-node science gate:

```bash
python3 mpi_review.py --prepare --parent /absolute/accepted/pair-campaign \
  --root /absolute/new/mpi-campaign
```

The coordinator may submit `mpi-campaign/scripts/mpi.sbatch mpi-campaign`. Only the two candidate undulator/free-drift cards run; their references are the existing one-node16MPI×8 candidate results. Each of two exclusive128CPU nodes runs8ranks×8threads, using128 cores in total while reserving256. Binding must show two hosts,64 distinct physical cores per host, and ranks0–7/8–15 on separate hosts. Physical cards and input remain identical; only output/reference directories change. Acceptance checks normal physical stop, all31675 particle IDs at each plane and at one terminal exit, complete finite HDF5, queue activity, final moving-window geometry/counts, whole trajectories, the original E/B2e-12 and particle128eps gates, and the existing ledger contract. The limit is240s per case,720s active total and15min Slurm. `MPI_PASS.json` establishes functional cross-host replay, not a performance gain from reserving two nodes.

The original short science card initially queues only15 of31675 particles. A separate majority-pending lifecycle pair exercises repeated admission and exhaustion of the queue:

```bash
python3 long_queue.py prepare --parent /absolute/accepted/pair-campaign \
  --root /absolute/new/long-queue-campaign --stretch 8
```

The coordinator may submit `long-queue-campaign/scripts/long_queue.sbatch long-queue-campaign`. It derives a new small input and reference-domain/stop geometry for the longer bunch while preserving the fine cells, boost,64/128 guards, solver and diagnostics in both variants. At least half the cohort must initially wait, and every particle must subsequently appear at all three planes and at one terminal exit. Both runs use the same16MPI×8 single-node allocation,300s each,840s active total and15min Slurm. The original production input is untouched. `LONG_QUEUE_PASS.json` covers that controlled longer-bunch lifecycle; it is distinct from full-input throughput and radiation convergence.

Recorded use on2026-10-28 is under workspace `results/pending_admission_20261028/`:

- `pair-r1/SCIENCE_PASS.json`:18977 completed the eight-case science gate. Its immutable scripts/results remain preserved.
- `mpi-r1/MPI_PASS.json`:18980 completed both257-step candidate cases on node4/node8; the raw E/B arrays were observed bitwise equal to the one-node16rank candidate. Each case initially queued15 particles.
- `long-queue-r1/LONG_QUEUE_PASS.json`:18981 completed both716-step cases; the candidate initially queued31669/31675 particles. Relative field L2 differences were about4.53e-14, passing the2e-12 gate; these arrays were not bitwise equal.
- `contract-r2/18979_BASELINE_REVIEW.json`:18979's baseline profile completed20 updates,21 committed field frames and3 ledger records, all with `complete=0`. Its subsequent expected MPI abort produced a Slurm step-cancel banner that the old broad regex misclassified. The candidate had not started. The corrected read-only review accepts only the bounded timing/data contract, not physical completion.
- `contract-r2/test.log` (identical saved output in `CONTRACT_TESTS.log`):17 offline regression/negative-control tests passed. They cover the observed log, removed close/readiness markers, wrong steps, early cancellation, unexpected rank/abort errors, timeout/signal/crash/OOM, incomplete/truncated/nonfinite HDF5 and incorrect completion flags. Run with `python3 test_bounded_contract.py -v`; no MPI is launched.
- `full-r2/TIMING_REVIEW.json`:18983 passed all four fresh full-input bounded cases and both comparisons on node3, explicitly reusing `pair-r1` science. All outputs remain `complete=0`, with `simulation_complete=false` and status `TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED`. The separate profile-disabled100-step pair measured2.792682s/step baseline and0.310372s/step candidate, an88.886% reduction for this early sample; the50–100 interval measured2.821342 and0.312441s/step. Initialization is separate, and these results do not establish a whole-production speedup or radiation acceptance.

Independent post-run reviews (no simulation launch):

```bash
python3 review_long_queue.py --root /absolute/completed/long-queue-r1
python3 review_performance.py --root /absolute/completed/full-r2
```

The long-queue audit is deliberately specific to the recorded8x-arrival fixture and also reads sibling `pair-r1`/`mpi-r1` evidence. The performance audit keeps progress and final resource timestamps separate: the original driver numbers above include final ledger/output close time, whereas progress50→100 isolates the logged stepping interval. That interval measured2.81606758→0.3120014792s/step (88.92067% reduction,9.02581x) on18983. The100-step output field is nonzero but extremely weak before the detector signal arrives; its bitwise equality is an early replay check, not a radiation-spectrum acceptance. Complete-lifecycle particle/field evidence comes from the short and long-queue science pairs.
