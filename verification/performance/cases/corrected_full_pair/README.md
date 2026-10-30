# Corrected full-input pair

User-authorized replay of the completed undulator 4.5 m and free-drift 0.5 m
cases. Each retains the original 41,732,775 input events, source, grid,
64/128 guards, field sampling and physical stop. Enable the accepted
reference_surface_current and inlet_charge_match. Use the accepted physics
binary directly; no kernel changes or build.

Only two raw particle planes: entry -0.05 m, exit coincident with the raw field
plane (4.5 m / 0.5 m). Independent energy-transport planes at the same locations
preserve scalar endpoint four-momentum statistics, since gamma of interpolated
raw three-momentum is a different diagnostic. Keep the boost v5 ledger and
strict discrete_operator_audit. The latter has a real runtime cost.

prepare.py freshly hashes the entire original input, checks the existing
same-binary 78-required acceptance, reuses completed full-input geometry and
compares every changed card field against an explicit whitelist. Diagnostic
worldlines must lie inside the active window to preserve pending scheduling.
It freezes copied scripts, cards, geometry, binary and acceptance identities.

run.sbatch uses one available exclusive 128-core node, 16 MPI x8, per case,
without a node list. Preserve the original 144h outer limit, MPI-free 600/60s
shutdown watchdog and self,vader,tcp. Read-only heavy probes are unset.
The full post-close checker runs on that allocation, scanning E/B and particle
outputs, unique IDs, energy-plane counts and the discrete ledger. Completion
is not physical or pure-radiation acceptance.

Campaign:
  results/corrected_full_comparison_20261030/
Do not edit copied frozen files while jobs are running. Confirm initialization,
both inlet options, exact source count, pending queue, binding and initial
progress from logs/manifests without opening running H5, then leave the long
jobs to finish until the user requests analysis.
