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

All twelve compressed CPML histories are handled individually. The current
contract retains history at the same physical position only where both old and
new supports are active, zeros newly active histories, and discards unsupported
history. It never transfers values between derivative terms. The original CPML
object, coefficients and recurrence remain in place. This is a stated moving
boundary convention; the old fixed-grid implementation does not uniquely define
the convention when its physical boundary starts moving. State-copy equality is
not proof of identical absorption by a moving and a stationary boundary.

Build with `APRL_BUILD_MOVING_WINDOW_PROBE=ON`; executable
`verification/aprl_moving_window_probe OUT MODE [WIDTH]`. Modes:

| Mode | Purpose |
|---|---|
| `verify` | Independent coordinate oracle, 126 labelled storage cases, four core PIC schedules, and missing-history negative control |
| `bench 64/128/256` | Same-width fixed/compact/snapshot timings; three rotated repetitions, four warm-up and 64 timed steps |
| `wave` | Independent 512-cell large domain vs 128-cell moving domain, 8-cell longitudinal CPML |
| `wave-static` | Nearby stationary upper CPML control; 296 cells, origin −192 |
| `wave20` | Existing production longitudinal thickness 20; moving width140 retains upper CPML entrance120 |
| `wave-static20` | Fixed control with thickness20 and width308 retains upper entrance96 |

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
assessment. Further changes to history attachment/reinitialization conventions
require the user's confirmation; none are silently introduced here.
