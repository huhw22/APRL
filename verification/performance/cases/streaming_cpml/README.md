# Charged-bunch streaming window with original longitudinal CPML

This verification-only integration combines the already tested integer
`world-position` history remap with the existing fixed-ID particle-entry queue
and once-computed Cowan reference cache. It preserves the original Cowan field
advance, Boris particle advance, current deposition, precision, prescribed
release, and CPML recurrence. It does not add a production simulation option.

The original transverse domain remains 24 x 24 cells with PEC. CPML is enabled
only along z, with cells `[0,0,20]`, to isolate longitudinal-window integration.
The full control keeps the original 384-cell nonabsorbing domain `[-64,320] um`
and adds 20 absorbing cells outside each end: 424 total cells, origin `-84 um`.
The candidate has `16 + 2*guard` nonabsorbing cells plus 40 absorbing cells,
giving total widths 88/104/136 for guards 16/24/40. A guard is a nonabsorbing
distance; it is not the CPML thickness. These toy guards are not production-card
recommendations.

All cases use 1024 generated Gaussian particles, seed 710271, total charge
`-0.1 pC`, the original dx=dy=2 um, dz=1 um, and 176 steps. Particle box beta is
0 or .23. The window follows the fixed lab release/detector planes at `-0.8c`
in box coordinates. Each particle's prescribed release position and momentum
remain identical across modes. Existing post-detector ballistic motion and the
two-cell deposition collar are retained; all admitted particle charge/current
uses the original deposition rules, including passage through the absorbing
layer. Later omitted exterior return fields remain an approximation.

The boundary is total-field CPML with frozen-reference-driven outer faces and
halos. The reference is not exempted from absorption. This is not residual-field
CPML. Both control and candidate use the same CPML parameters and Cowan
initialization guards. The cache contains only the initial full-domain solve;
the independently evolving full control never supplies its later fields to the
candidate. `zero` retains the common initial slice but removes later reference
supply as a negative control.

Run preparation and execution
-----------------------------

Enable `APRL_BUILD_STREAMING_CPML_PROBE=ON`, build `aprl_streaming_cpml_probe`, and
pass its absolute executable path explicitly. The executable accepts
`OUT beta guard full|cpml|zero`, requires one MPI rank, and refuses to overwrite
scientific output files.

```
python run.py prepare --root FRESH_ROOT --binary /absolute/path/to/aprl_streaming_cpml_probe
```

The default is only the two initial preflight cases:
`b0-full-g24:0:24:full` and `b0-cpml-g24:0:24:cpml`.
After these are independently reviewed, prepare a fresh campaign with explicit
`--case name:beta:guard:mode` arguments for the remaining eight cases: dynamic
full and CPML guards16/24/40, static CPML guards16/40, and zero-supply guard24 at
both betas. The runner does not automatically advance to another campaign.

Execute the frozen `FRESH_ROOT/run.py run --root FRESH_ROOT` on the authorized
gnode1 CPU host with `APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1`. Shared-lock, idle,
binding, immutable-source/binary, process-group cleanup and timeout safeguards
match the moving-window runner. Each case has a 180-second watchdog; a campaign
has a two-hour total limit. A rejected duplicate invocation cannot overwrite
another invocation's active or completed result.

The runner gates finite output, complete cohorts, the original interior
continuity criterion (`1e-10`), and open/remapped charge budgets (`1e-10` of
total charge). Gauss residuals, including initialization, are recorded without
a numerical acceptance threshold. They are not used to justify changing the
original discretization or other numerical compromises.

Ordinary Gauss diagnostics cover only the nonabsorbing interior vertices.
The legacy CSV name `new_layer_gauss_relative` measures the first nonabsorbing
interior layer, not the newly inserted CPML plane. Reference-Gauss diagnostics
also restrict to cache-valid overlap. These observations do not impose the
ordinary vacuum divergence constraint inside the absorbing layer.

Independent review and figures
------------------------------

```
python review.py --runs PREFLIGHT_ROOT [FOLLOWUP_ROOT ...] --output REVIEW.json
python plot.py REVIEW.json --output PLOTS
```

The runner freezes these self-contained Python scripts, README, all relevant
verification C++ files, CMake configuration, production source hashes, binary,
runtime libraries and binding probe. The reviewer rechecks provenance and
process completion; sorted fixed-ID initial/release values must match in both
numeric bit patterns and canonical CSV text. It independently recomputes stable
energy changes from proper velocity, event Lorentz transforms, detector moments,
energy spread, expansion, and arrival-current bunching. It also checks queue
accounting, field cadence/pixels, and integer origins. It reports measured field
and particle differences without an invented physical acceptance threshold.

Mean particle delta energy is signed: positive means gain; mean energy loss is
its negative. Bunching uses fixed-plane arrival times at lambda=1/5/10 um, not a
simultaneous spatial snapshot. Field data retain native E(t_n), B(t_n-dt/2),
both sampled at the detector position at t_n. They are not simultaneous lab
six-component data and cannot directly establish lab Poynting flux or pure
radiation. Gauss observations remain available in the report without gating.

Timing is one run per different-sized domain, including the initial full solve,
reference cache, deposition collar, particle queue and expensive research
diagnostics. Report these as workload observations, not production implementation
speedups. This task does not establish asymptotic scaling, multi-node behavior,
weak-radiation accuracy, or full production checkpoint/diagnostic integration.
