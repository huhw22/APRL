# Charged-bunch moving-window research probe

This verification-only target tests an independently supplied, translating short
PIC window against a complete-domain simulation. It does not add a production
window option or change the accepted Cowan initializer, particle release rules,
Maxwell solver, or charge-conserving depositor.

Enable `APRL_BUILD_STREAMING_BUNCH_PROBE=ON` and build
`aprl_streaming_bunch_probe`. The executable accepts
`OUT beta guard full|clamp|mur|mur-gauss|zero`, requires one MPI rank, and refuses
to overwrite scientific CSV files. Generated input is fixed at 1024 particles,
seed 710271, charge -0.1 pC, sigma_xy=3 um and sigma_z=12 um. The 384-cell complete
domain and 48/64/96-cell windows use the same dx=dy=2 um, dz=1 um and 176 steps.
The window moves at -0.8c in simulation coordinates, following the fixed lab
release/detector planes; the boost has beta .8 and gamma 5/3. Particles have
box beta 0 or .23. Each particle's prescribed release position
and momentum are identical across cases.

The short window reads only a deep copy of the once-computed full Cowan initial
field. It does not read the complete simulation's later fields. Reference values
use native lattice and time staggering, with linear translation at beta_box.
Total fields evolve with the original Maxwell kernels and actual deposited J;
the reference supplies external halos, end conditions and new grid planes.
Particles enter through a sorted queue. A two-cell deposition collar preserves
the physical trajectory segments needed for interior current continuity. Moving
the grid has a separate geometric charge budget. Outgoing particles become
ballistic only after detection and are removed outside the collar; their later
return fields are therefore an approximation, not a claim of exact equivalence.

`clamp` supplies reference end fields. `mur` applies first-order Mur to the
tangential residual. `mur-gauss` additionally reconstructs only the new normal
electric edge to satisfy the newly interior vertex's Gauss constraint. `zero`
is a deliberately inconsistent negative control that omits subsequent reference
inflow, while retaining the common initial field slice.

Run `run.py prepare --root FRESH_ROOT [--case name:beta:guard:mode ...]` locally,
then the frozen `FRESH_ROOT/run.py run --root FRESH_ROOT` in an authorized gnode1
CPU session with `APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1`. The runner checks a shared
lock, host idleness, binding and immutable inputs; every case is limited to 180
seconds and the campaign to two hours. It stops on nonfinite metrics, incomplete
cohorts, initialization failure, continuity or charge-budget errors. It does not
convert those checks into physical acceptance.

`review.py --runs ROOT [ROOT ...] --output REVIEW.json` independently checks the
fixed-ID inputs/releases, event transforms, energy formulas, output completeness
and charge budgets. It reports detector particle moments and native field
differences. `plot.py REVIEW.json --output PLOT_DIR` exports PNG/PDF figures.

Particle momenta at plane events are linearly interpolated on their crossed
trajectory segment. Bunching describes fixed-plane arrival-current times with
lambda=c/f; it is not a simultaneous spatial bunch snapshot. Field CSVs contain
box E(t_n) and B(t_n-dt/2), both sampled at z_plane(t_n). They are not aligned lab
six-component data and must not be used directly for lab Poynting flux.

The 2026-10-27 experiment comprises 18 bounded runs, recorded in workspace
`results/streaming_bunch_20261027`. All listed cohort/charge consistency checks
passed. The moving beta=.23 Mur variants exhibited large growing field errors
and were rejected. Direct reference clamps remained bounded over this test and
gave close particle expansion and energy spread, but native E/B differences
remained about 4--5%/17--21%; average energy exchange must also be checked rather
than hidden by its cancellation. These results do not accept weak radiation or
production deployment. The same reference translation has about 0.23 normalized
Gauss mismatch with re-deposited CIC density even in the complete-domain control's
diagnostic, identifying a separate reference/source consistency issue.

Costs include the full initial solve/cache, active fields, deposition collar,
particle input and reserved capacity, queue and flags, supply queries, remaps and
diagnostics. The full probe also allocates research-only reference/deposition
arrays and runs expensive full-domain diagnostics. A total wall-time ratio is
not production speedup. Only one bunch length is tested; fixed-error window-size
scaling, weak-wave preservation, distributed boundaries and production checkpoints
remain unverified. Do not reuse these toy guard counts as production-card values.
