Energy-transport diagnostic acceptance (2026-10-29)

This directory is an independent successor to moving_energy_ledger. It tests
diagnostic additions without changing Cowan/Boris/current deposition, reference
supply, release events or the original particle-plane files. No production input
or long simulation is launched by these tools.

The fixed baseline is the accepted v4 energy-ledger binary
build/moving-energy-ledger-r2-20261029/aprl (SHA490e4e6663418efcdce0bb84fb8c77248fbe0f7c38a9c8d220a5d4874a124fd1).
The candidate supplies v5 bookkeeping and independently optional fixed-lab
energy_planes. The old pending-admission v3 binary is not the baseline here.

Prepare a new reviewable plan before a candidate exists:
  prepare.py --workspace /public/home/wfollowm/APRL_Workspace --root NEW_ROOT
Once the exact candidate and frozen runtime-libs exist:
  prepare.py --root NEW_ROOT --finalize --candidate ABSOLUTE_NEW_BINARY
Alternatively specify the existing candidate during initial preparation to
freeze immediately. An unfrozen plan cannot run. A frozen plan is never edited.

Science uses one exclusive node, fixed16 MPI x8 threads, and the existing
31,675-record short input, original fine dz/boost/64-128 guards/xyCPML. One-period
undulator and free drift each run old-on/new-on/new-off/new-audit-on: eight runs,
257 physical steps apiece. The audit adds entry(-0.05m) and exit(0.06m) energy
planes while retaining all legacy E/B, three particle planes and trajectories.
All legacy numeric payloads must match exactly within a physical case/layout.

The performance stage requires the same-binary science gate. It retains the
same source/geometry; trajectory and particle-file output are disabled, raw E/B
and ledger output use a10-step cadence. Energy planes integrate every step and
initially also write at full cadence, for independent numerical integration.
Three alternating-order repeats of all four variants for both physical cases
give24 runs. Report old-on/new-on, new-on/new-off and new-audit-on/new-on costs
separately, with every paired sample, medians and extrema. Profiles stay off.
The disabled-path cost is measured against new-on; there is no old-off pair,
so do not claim a precise new-off versus old-off overhead from this campaign.

MPI requires the same science gate and uses two exclusive nodes, the same
total16x8 layout (8 ranks/node,64 actually used cores/node). Each physical case
runs new-on/new-audit-on, giving four within-layout invariance checks. Do not
hide an on/off mismatch behind cross-layout tolerances. No timing comparison
between the one-node and two-node stages is an implementation speedup claim.

Each case has a180s execution limit. Stage active budgets are1500/2700/840s;
Slurm hard limits are50/50/15min (science/performance/MPI); the unit job has20min for required tests plus explicit2/4-rank modules. Root alone submits run.sbatch ROOT science or
performance, mpi.sbatch ROOT, and unit.sbatch BUILD NEW_RESULT_DIRECTORY. First
check the current queue/resources. Keep the600/60s local watchdog and job-local
self,vader,tcp configuration. No automatic retry, requeue, node pin or submission.

The standalone module tests cover pure magnetic energy/crossing diagnostics,
kick-work decomposition and deliberate malformed accounting. Boris checks run with
the required suite. Energy-plane checks run explicitly at2/4 ranks with fresh
module directories; original2/4-rank full-volume remap checks are retained.
No new BF physical card is included in this functional campaign.

Acceptance distinguishes complete/unchanged physics, diagnostic identities,
and the actual nonzero physical imbalance. Full-time peak and terminal residuals
are kept, normalized to actual exchanges and net particle energy change, not
only initial total kinetic energy. closure_valid is never promoted by review.
The fixed-lab inlet/outlet sum is a transport comparison, not a closed budget
unless stored energy, lateral fluxes and sources for that same space-time volume
are all available. All fields are total self-consistent fields, not isolated
near fields. Changing a diagnostic is not permission to change the physics.

R2 retains the same physical cards and v4 baseline. It requires the candidate
field boris_absolute_work_J and tests zero partial_work_substeps before the
first removed particle. R1's failed net-sum identity is preserved; the revised
identity bound uses independently measured historical absolute work, not an
increased coefficient or a residual-derived correction. New preparations
must use a fresh campaign directory and the separate r2 candidate binary.

R3 uses globally unique detector names audit-entry/audit-exit; plan.json maps
logical entry/exit to the actual energy filenames. --plan-only keeps an
existing candidate unfrozen. --reuse-science OLD_ROOT can retain completed
old-on/new-on/new-off science cases after same-binary/input/physics/layout,
transaction/finite/lifecycle and exact output hashes are checked. All original
files remain read-only, including CSV. New case directories contain only the
new planned card and new review; reused_from points to original data and the
real original node remains recorded. Failed/unfinished cases are never reused.
Timing and MPI stages always launch fresh cases. Reused scientific runs are
not presented as newly measured timing on the new allocation.
