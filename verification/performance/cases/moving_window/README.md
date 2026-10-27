# Integer moving-window state and cost checks

This optional verification target moves storage while retaining the original
Cowan, particle, deposition and CPML update implementations. No `src` file or
production input changes are required. It is not yet a production moving-window
configuration switch.

`MovingWindowRemap` keeps native contiguous arrays. For a signed origin shift
`s`, new plane `k` takes old plane `k+s`. All overlapping E/B/J values are copied
without floating arithmetic. The caller maintains an integer world-cell offset
and supplies the new origin and only newly entering planes. Complete steps end
before migration, preserving the native E/B and CPML time staggering. Array
allocation addresses remain unchanged; precompiled block plans keep temporary
allocation and indexing work out of the normal migration path.

All twelve compressed CPML histories are handled individually. The default
contract retains history at the same physical position only where both old and
new supports are active, zeros newly active histories, and discards unsupported
history. It never transfers values between derivative terms. The original CPML
object, coefficients and recurrence remain in place. This is a stated moving
boundary convention; the old fixed-grid implementation does not uniquely define
the convention when its physical boundary starts moving. State-copy equality is
not proof of identical absorption by a moving and a stationary boundary.

The user subsequently authorized a controlled alternative: `BoundaryRelativeZ`
retains just the four longitudinal derivative histories (indices 1,2,7,8) at
their boundary-relative slots. The other eight histories and all E/B/J arrays
retain the original mapping. `WorldPosition` remains the default. Both use the
original CPML recurrence and coefficients. `historyRetainedBytes` records the
storage kept in place, separately from logical copied/zeroed payload; it is not
read/write traffic or additional allocation.

Build with `APRL_BUILD_MOVING_WINDOW_PROBE=ON`; executable
`verification/aprl_moving_window_probe OUT MODE [WIDTH]`. Modes:

| Mode | Purpose |
|---|---|
| `verify` | Each policy: independent coordinate oracle, 126 labelled storage cases, four core PIC schedules, and missing-history negative control; also a direct zero-shift identity / nonzero-history distinction check |
| `bench 64/128/256` | Same-width fixed/compact/snapshot timings; three rotated repetitions, four warm-up and 64 timed steps |
| `bench-boundary 64/128/256` | Same benchmark with boundary-relative longitudinal histories |
| `wave` | Independent 512-cell large domain vs 128-cell moving domain, 8-cell longitudinal CPML |
| `wave-static` | Nearby stationary upper CPML control; 296 cells, origin −192 |
| `wave20` | Existing production longitudinal thickness 20; moving width140 retains upper CPML entrance120 |
| `wave-static20` | Fixed control with thickness20 and width308 retains upper entrance96 |
| `wave-boundary` / `wave-boundary20` | Same respective moving wave, only longitudinal-history attachment changed |
| `wave-long` / `wave-long20` | Same candidate window for 384 steps; reference1024, origin−512 |
| `wave-boundary-long` / `wave-boundary-long20` | Same extended observation with boundary-relative longitudinal histories |

The benchmark coordinate oracle uses independent full snapshots and independently
derived CPML indices. Its buffers are preallocated. Modes allocate only their
own auxiliary structures; preflight objects are released before timing. Eight
steps per width are compared bitwise before timing; the timed 64-step trials
are checked for finite states. The PIC schedules compare B half steps, complete
E steps, all histories and particle position/momentum bitwise, but do not cover
the complete `Simulation` held/MPI/checkpoint workflow.

The wave test initializes both domains from the same compact discrete vector
potential and advances them independently using original kernels. The moving
domain receives zero incoming vacuum layers and never reads later large-domain
fields. Diagnostics use the common trajectory `z/dz=100−0.8*n`, fixed in the lab
for the stated boost. Native E and B times are retained; these are not simultaneous
lab six-component fields or a direct lab Poynting measurement.

The 20-cell cases use the existing production longitudinal thickness, not a
parameter search. They retain the two-dimensional common transverse fixture
`x6/y0`, so they are not the complete production configuration. Fixed and moving
upper boundaries have different trajectories. Their finite-time late responses
are reported separately, not subtracted as isolated reflection coefficients.

`run.py prepare --root FRESH_ROOT [--binary PATH] [--case NAME ...]` freezes
source, executable and runner. Run its frozen copy on authorized gnode1 with
`APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1`. Shared-lock, idle and physical-core binding
checks precede one MPI × one thread; cases are limited to 180 seconds and the
campaign to two hours. Default cases are verify and the three benchmarks.

`review.py --cost COST_ROOT --wave WAVE_R1 --controls CONTROLS_R2 --output FILE`
independently checks records and recomputes timing and wave comparisons.
`plot.py FILE --controls CONTROLS_R2 --output DIR` exports PNG/PDF figures.

2026-10-27 results are in workspace `results/moving_window_20261027`. Migration
checks passed. Same-width extra wall time was 2.46%, 3.67%, 4.24% on one physical
core; native-kernel time varied by at most 0.17%. Helper storage was 8,808 bytes.
These are measured synthetic overheads, not a distributed production speedup.
The original 20-cell CPML moving wave matched the large-domain sampled fields for
steps0–99. Its late difference peak was 0.832% of the main pulse plane peak;
the nearby fixed control was 0.153%. The moving boundary still needs physical
assessment. The original result files and executables remain frozen. The new
two-policy comparison is separately recorded in
`results/moving_window_history_20261027`; neither policy changes production.

The extended comparison keeps `WorldPosition` as the preferred tested baseline.
For the existing 20-cell longitudinal layer, its 384-step relative L2 difference
is 2.317% and peak plane difference is 0.832% of the main pulse; the alternative
gives 13.463% and 4.131%, respectively. The alternative's peak arrives at step198,
outside the original176-step observation. Its smaller short-interval L2 did not
mean suppression. The 8-cell fixture also favors the original policy. Both
rules passed their distinct 126 state cases and four PIC equivalence schedules;
these checks do not assert physical equivalence between the policies.

The longer reference was chosen using a conservative two-z-cell complete-step
dependency bound; with source/interpolation margins its nearest possible end
return is after step403, beyond the384-step end. Every extended CSV reproduces
its original177-frame prefix byte for byte. This reference expansion does not
enlarge the moving candidate or a production window. The prototype remains
single-rank; charged-beam/reference-inflow integration is still separate work.

`review_history.py --root SHORT --previous OLD_CONTROLS --output REVIEW.json`
audits the11-case dual-policy campaign. `review_history_long.py --root LONG
--short-review REVIEW.json --output LONG_REVIEW.json` audits four extended waves.
`plot_history.py REVIEW.json --long-review LONG_REVIEW.json --output DIR` exports
the comparison figures. Use an explicit binary for frozen historical campaigns;
the runner's current default points to `build/moving-window-history-long-20261027`.
