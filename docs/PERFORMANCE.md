# APRL computation-only performance checks

Use the same physical configuration to compare executables and CPU layouts. Do not change particle number, mesh spacing/extent, solver, current deposition, particle substeps, field feedback, boundaries, physical stopping condition, detector cadence, or initial field accuracy to claim an implementation speedup. The first implementation changes preserve the order of particle deposition; same-rank output comparisons therefore require bitwise equal numeric values.

## Implementation scope

Grid/CPML updates and independent particle preparation may run on multiple CPU threads. Current deposition and particle-event handling retain their original serial ordering to preserve same-rank bitwise results. The implementation also caches unchanged held-start velocities, reuses work buffers, skips work for empty/disabled features, reduces unnecessary per-step statistics collectives, and overlaps a nonblocking magnetic halo exchange with the internal electric update. These changes keep the physical input and discretization unchanged. MPI and HDF5 calls remain on the initialization thread.

Profiling is off unless requested. Stage reports distinguish `field_total`, `particle_prepare`, `particle_deposit`, `particle_migration_total`, `diagnostics`, `control`, `field_communication` and `particle_communication`. The two communication measurements are nested within the corresponding totals. `particle_deposit` also includes per-particle boundary and detector-event handling; `diagnostics` covers event collection/HDF5 and scheduled sampling. Communication time can include waiting for a slower rank. Do not sum parent and child stages, subtract maxima taken on different ranks, or interpret MPI-call time as pure network bandwidth. A synchronization after reporting prevents the report's logging cost leaking into the next measurement interval. Use profiling to locate work; use profiling-off runs for throughput conclusions.

## Open transverse reference preparation

Moving-window runs with `open_transverse_reference: true` automatically distribute
the Coulomb boundary-potential targets across the simulation's MPI communicator.
Each rank uses its configured OpenMP threads. This includes the coarse potential
table, independent interpolation checks, and full native-node fallback when a
tile fails the existing interpolation tolerance. No input-card switch is needed.
For example, 16 MPI ranks with 8 threads each can use all 128 allocated cores in
this stage. Keep the existing per-rank CPU affinity and `OMP_NUM_THREADS`; do not
give every rank the node's full thread count.

Each target retains the original source summation order. Ranks replicate the four
arrays of nonzero deposited mesh charges (32 bytes per source node) and return
potential values to rank zero in target order. There is no floating-point MPI
sum over sources. Deposition, interpolation, the blocked modal solve and file
publication remain on rank zero; only the expensive direct potential targets
are distributed. Each worker subsequently opens its own read-only reference.
Every run still prepares a fresh reference for its actual geometry and beam.

The blocked reference validates positive definiteness using the actual finite
grid: after the transverse sine transform, each longitudinal block has constant
diagonal `d` and off-diagonal `o`; its least eigenvalue is
`d - 2*abs(o)*cos(pi/(n_z+1))`. Every block must have a finite positive minimum
before files are created, and every elimination pivot is still checked. The
stronger per-direction sufficient condition used by the generic CG initializer
is not necessary for this modal solver and rejects some positive coarse-grid
operators. This changes admissibility checking, not the potential equation,
source shape, field reconstruction, or boundary data. The generic CG path keeps
its existing criterion. Coarse-grid regressions separately enumerate the full
three-dimensional sine spectrum and check native `div(E)=rho/epsilon0`.

`reference/preparation.json` reports `openTraceRanks`,
`openTraceThreadsPerRank`, `openTraceSourceNodes`, `openTraceCoarsePoints`,
`openTraceFallbackTiles`, and the minimum/maximum per-rank potential evaluations.
These counters distinguish useful potential work from MPI waiting. The existing
`openTraceSeconds` includes source extraction, communication and boundary-file
completion; `preparationSeconds` also includes the remaining preparation stages.
The blocked solve has not become a distributed solve, and whole-initialization
speedup must be measured separately.

The required `reference.parallel-*` checks compare binary source, potential,
checkpoint and boundary files plus all native E/B values against singleton
preparation. They also cover two thread counts, all source shapes used here,
signed-source fallback, empty sources, failures and subcommunicators.
`aprl_reference_trace_benchmark` is an optional-use timing executable built with
the moving-window verification targets. A shortened target interval is useful
for a bounded kernel check; its timing must not be reported as a full production
reference speedup.

## Small repeatable SASE-loading harness

`verification/performance/run.py` makes one seeded HDF5-v2 input with the existing `sase_bunch_generator`, then runs an unchanged baseline and candidate in separate directories. All input cards are byte-identical. This is a short coupled-particle/field computational check with a small undulator, **not** a realistic SASE gain or radiation-convergence experiment. It does not enable held startup, which requires HDF5-v4 plane events; the dedicated frozen-plane regressions cover that path.

Use Python with `numpy` and `h5py`. Keep the baseline binary before rebuilding. For example, from the APRL source root:

```bash
python3 verification/performance/run.py \
  --baseline /absolute/baseline/build/aprl \
  --candidate /absolute/candidate/build/aprl \
  --generator /absolute/candidate/build/sase_bunch_generator \
  --output-dir work/performance/local-plan \
  --ranks 1,2 --threads 1,2 --repeats 3 --plan
```

`--plan` creates only cards and `plan.json`, with no input generation or simulation. Use a **new** output directory and omit `--plan` to run. Existing output directories are never reused. Defaults are 1/2 MPI ranks, 1/2 candidate threads, 2 repetitions, seed 17 and 16 beamlets per wavelength. Baselines always use one thread. The generator creates 2,560 records with the default settings. `--beamlets` can create a larger compute test, but the chosen value stays fixed throughout each comparison. It must remain consistent with the fixed 100,000 physical electrons and the generator's shot-noise constraints.

The default `mpiexec` path uses Open MPI options `--bind-to core --map-by slot:PE=N`. Select `--launcher none --ranks 1` for a one-process direct launch, or `--launcher srun --srun-mpi pmix_v3` inside a Slurm allocation. For an unbound local correctness check, `--binding none` disables MPI CPU binding and defaults OpenMP binding to false, preventing multiple ranks from being pinned to the same first CPU; it is not a controlled cluster timing configuration. Extra launcher arguments are explicit repeated `--launcher-arg=VALUE`. Do not oversubscribe timing runs. No network-related environment variables are changed by the harness. BLAS helpers are restricted to one thread, and OpenMP uses cores/close binding unless already explicitly configured.

For Slurm, load the MPI environment used for both binaries, then submit from the source root:

```bash
export APRL_BASELINE=/absolute/baseline/build/aprl
export APRL_CANDIDATE=/absolute/candidate/build/aprl
export SASE_GENERATOR=/absolute/candidate/build/sase_bunch_generator
export PYTHON_BIN=/path/to/python-with-h5py
sbatch --partition=YOUR_PARTITION verification/performance/run.sbatch
```

The wrapper requests one node, two MPI slots and two CPUs per slot by default. Explicit `sbatch` options and `RANK_COUNTS`, `THREAD_COUNTS`, `REPEATS`, `BEAMLETS`, `SRUN_MPI` and `PERFORMANCE_OUTPUT_DIR` control the resource experiment. The allocation must fit the largest requested rank/thread combination. The smoke card has 128 z cells and supports no more than 64 ranks; use a separately validated production case for larger cluster scaling. Slurm execution has to be checked at the actual site.

## What the result verifies

`report.json` records executable and input SHA-256, commands, selected thread/network environment, launcher CPU affinity, actual ledger step count, whole-run wall times, repeated-run timing distributions, final loop timing and candidate `[performance]` records. Repetition order reverses to reduce systematic run-order bias. Whole-run wall time includes initialization; `loop_mean_step_s` comes from the final resource report, whose timer starts at stepping but also includes final diagnostics/output closure. Stage profiles provide more detailed timing. The small default run is useful for correctness and integration; its timings are not a production speedup estimate.

Every run must reach the configured physical stop, produce complete committed nonempty ledger/particle/field HDF5 datasets, contain finite numeric output, lose no particles and record every generated particle exactly once at the endpoint. Particle records are matched by `source_id`. All numeric output fields, integer counts and recorded time axes are compared, not only exit status.

- `implementation_equivalence_pass` requires bitwise agreement against the unchanged baseline at the **same MPI rank count**, including different candidate thread counts and repeated runs. Any failed run also fails this result.
- `cross_rank_consistency_pass` separately compares outputs to the first baseline rank count. It uses an explicit per-field relative tolerance (`--rtol`, default `1e-9`) and absolute tolerance scaled to the reference field magnitude (`--scaled-atol`, default `1e-12`). Each field reports maximum absolute/scaled error. Zero reference fields require exact zeros.
- Overall `pass` requires implementation equivalence and an unchanged shared input. Add `--require-cross-rank` to also require the cross-rank tolerance check. Existing baseline domain-decomposition differences remain visible; they must not be silently hidden by increasing tolerances or mistaken for a new optimization regression.

`APRL_PROFILE_INTERVAL=N` enables candidate stage timing every N steps; 0/unset disables it. The harness uses `--profile-interval 16` by default and passes 0 to the baseline. Run a second repeated series with `--profile-interval 0` for throughput measurements without profiling overhead. This does not change the card or output sampling. `field_communication` is included in `field_total`, and `particle_communication` in `particle_migration_total`; do not add parent and child times. MPI call time also contains waiting for slower ranks, so it is not a pure network-transfer measurement. Across-rank maxima for different stages need not occur on the same rank and should not be summed as an exact wall-time reconstruction. A result reports actual measurements, not promised acceleration.

## Production Elegant jobs

The independent `APRL_Elegant_4m_2nodes` scripts retain their physical card and original 2-node × 128-rank × 1-thread default. `THREADS_PER_RANK` sets both `--cpus-per-task` and `OMP_NUM_THREADS`; explicitly reducing `TASKS_PER_NODE` keeps the total allocated CPU count fixed for a hybrid MPI/thread comparison. Inspect the requested/effective thread count and per-rank launch affinity in `aprl.log`. The scripts retain the existing MPI network defaults and do not automatically select UCX.

## Local acceptance example, 2026-10-23

The retained baseline is `work/performance/baseline-7a8ff35/aprl` (source HEAD `7a8ff35`). Its SHA-256 is `821e427324b5d83eb26a86bf0f17d362f5217b27812b46daf111ecd44878e364`; the tested `build-performance/aprl` candidate is `79bf78a182a230eeda2c1ef397c0fd7d62c75999c843d5c2a673d38599d7a06c`. Rebuilding the candidate invalidates this exact binary association.

`work/performance/20261023-final-throughput/report.json` contains 24 profiling-off runs: 1/2 MPI ranks, baseline at one thread, candidate at 1/2/4 threads, three repetitions each. The table gives the median program-reported elapsed seconds, **including initialization**, rather than launcher time or the separately recorded loop timer.

| MPI ranks | Baseline, 1 thread | Candidate, 1 thread | Candidate, 2 threads | Candidate, 4 threads |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.636959 | 0.637580 | 0.378260 | 0.250712 |
| 2 | 0.370419 | 0.368983 | 0.227542 | 0.171826 |

`work/performance/20261023-final-profile/report.json` contains eight additional runs of the same matrix with one repetition and candidate profiling every 16 steps. Both series passed all run-completeness checks and all same-rank bitwise comparisons; every run recorded 2,560 particles and reached the physical stop after 19 steps. The stage records use the parent/child timing names described above. Cross-rank consistency fails the explicit default tolerance already for the unchanged baseline, and remains separately reported; these checks establish implementation equivalence at a fixed rank count, not cross-rank convergence.

These are initial short local measurements on `wfollowm-ACC-LAB`, using `--binding none`, an available CPU set of 0–15, and no overlapping APRL test jobs. Threads at a fixed rank count use additional CPU resources. The table is therefore neither a fixed-resource scaling result nor a production or multi-node speedup claim. A cluster comparison must retain the physical card, reserve an equal total CPU budget, check rank/thread affinity and network settings, and collect longer repeated runs. No production simulation or Slurm submission was performed for this acceptance example.
