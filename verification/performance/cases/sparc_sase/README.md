# SPARC SASE comparison

This harness compares explicit Boris/linear, conservative-source/linear, and
conservative-source/quintic-compensated longitudinal coupling against the same
imported GENESIS-4 physical shot-noise realization. It does not change Cowan,
particle/deposition kernels, the quick reference generator, or CPML.

The reconstructed six-segment SPARC lattice uses the 2017/2019 published beam
parameters. Quadrupoles are -9 T/m, length 0.053 m, beginning 0.105 m after each
preceding undulator exit. This explicitly resolves an ambiguity in the paper;
it is not the authors' unavailable original input deck. Each segment has a
one-period internal potential ramp and an independently calculated phase.

The quoted SPARC entrance distribution is mapped ballistically upstream by
0.01 m to satisfy the unchanged common-lab-time placement constraint. Both
programs propagate this preparation drift. Arrival events, transverse momenta,
weights, and optical phases are imported from the same realization; no second
noise realization or radiation seed is added. Initial noise is checked against
the physical electron count, including transverse quadrants.

`prepare.py` and `run_preflight.py` generate and check bounded input; only closed
GENESIS initialization dumps are read. `verify_preflight.py` can independently
validate closed timing data using ceil(duration/dt). Such timing outputs must
remain marked physically incomplete; they are not physical benchmark results.
The post-close MPI readiness contract rejects crashes, timeouts, watchdog exits,
and genuine initialization errors.

`prepare_main.py PREFLIGHT OUTPUT` freezes complete six-segment cards. Submit
`main.sbatch OUTPUT` and `reference.sbatch OUTPUT` separately. The PIC array has
three tasks with at most one active; the staged node count follows benchmark.json
(16 MPI ranks x 8 threads per node), with its own 144 h allocation limit. GENESIS uses one node, 64 MPI ranks,
and a 4.25 h allocation limit. No node is pinned and no automatic retry is used.
Original collective readiness and MPI-free shutdown watchdogs are preserved.

All self-consistent E/B components are exported at each segment exit and at two
free-drift planes 0.15/0.30 m after the last segment. Particle planes include the
upstream input and actual undulator entrance. Boost energy/operator ledgers and
independent laboratory transport planes remain enabled. Raw E/B include bound
fields and cross terms: filtering near resonance is not a near-field proof.

One seed is an exploratory paired comparison. It cannot establish ensemble
SASE agreement at 5%, and a complete file is not a physical acceptance result.
Published differences, model differences (including GENESIS longitudinal space
charge versus PIC), grid convergence, aperture, and bound-field contamination
must be reported separately before a precision claim.

For the user-approved cost screen, prepare with --slippage-multiple 3,
--transverse-cell-um 40, --transverse-cells 80, --nodes 1, and --reuse-tests PREVIOUS.
This preserves peak current, emittances and energy spread; charge follows the
shorter parabolic profile. The immutable benchmark.json configures staged
scripts. It is a short-bunch derivative, not the published long-bunch experiment.

Direction convention: local GENESIS positive slippage advances toward larger
slice index, whereas fixed-laboratory-plane arrival times increase toward the
tail. Active slice groups must be reversed on import, preserving each group's
optical phases and particle data. The 2026-11-04 radiation review corrects the
initial reference mapping and retains its earlier outputs as historical evidence.

The same review also found that BeamSolver replaces ku=0 in drifts with
ks/(2 gamma_ref^2). The reference phase shifter must include this existing
carrier compensation. The corrected nominal phi is about -0.9577053 rad;
the old +1.70418 rad produced about -152.5 degrees entrance mismatch per gap.
Closed-reference verification now checks the actual nominal phase cycle, in
addition to checking that the requested parameter was applied. This phase test
does not certify the full FEL physics or the 0.67 m approximate gain prediction.
