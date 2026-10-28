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

The R3 configuration rejects an early reference-center stop upstream of its4.5m detector. Therefore the full stage deliberately uses the existing bounded-duration diagnostic contract: exit1 with the specific duration-exhausted reason, exact N steps, closed readable HDF5, committed extents and `complete=false`. Its successful diagnostic status is **TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED**; it is never called a completed physical simulation. Cancel/timeout/watchdog errors, unexpected rank errors and truncated output fail. The previously accepted full-scale duration-exit validator is copied and frozen, including root output-close/readiness checks. The existing MPI-free600/60s shutdown watchdog is preserved in both cards.

Both stages reserve one exclusive128CPU node; the science1MPI case intentionally uses only8 cores for correctness coverage. All compared full runs use16MPI×8threads. Both variants set `OMPI_MCA_btl=self,vader,tcp` and the prior openib warning variable, plus close/cores thread binding. Every used layout is checked with the binding probe. No node is pinned.

Slurm has a40min hard limit; the driver has an active1800s stage budget and each full run a600s limit. A full run starts only if its entire600s allowance still fits. Budget exhaustion preserves `FAILED_OR_INCOMPLETE.json` and stops; it never extends the job or converts partial work into PASS. Small runs have240s per-case limits. Profiling and throughput are reported separately, as are reference preparation and stepping. The throughput result is one early100-step pair (including50–100-step interval), not repeated timing statistics or a full-run speedup forecast.

Scripts do not remove reference/data files and do not modify default binaries, live production cards or cluster configuration.

Energy-ledger validation is declared in `THRESHOLDS.json` before execution and applies to both science and full pairs. Integer counts and clocks are exact; `active_macroparticles` includes waiting particles, and active+removed must equal the original cohort at every sample. Source sums use `128*eps64 + 2*gamma_N(u_longdouble)` on noncancelling source scales (N fixed at31675 or41732775). Residual absolute bounds propagate the energy subtractands; ratios propagate numerator/denominator uncertainty, with an explicitly marked cross-product gate when the denominator is unresolved. Sigma-gamma is checked through central variance, mean uncertainty and the existing zero-resolution clamp. Axial chirp and uncorrelated sigma are reported without equality acceptance because pending z moments use an affine diagnostic trajectory. Their values never control admission or particle evolution; actual admission uses exact old repeated additions. E/B and per-ID x/u/t gates remain unchanged. Reports include error/source-scale and, where resolved, error/particle-energy-change metrics.
