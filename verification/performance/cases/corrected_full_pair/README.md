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

Identity rule: particle_id is the one-based runtime input-row index; source_id
is copied from the original SDDS particleID and may be sparse after cropping.
Check each runtime row exactly once, including duplicates across chunks, then
compare its source_id with that original input row. Never require source_id
to lie in 1..N. Missing rows, duplicate runtime IDs, wrong source mappings,
incomplete transactions, nonfinite records and wrong directions still fail.

Campaign:
  results/corrected_full_comparison_20261030/
Do not edit copied frozen files while jobs are running. Confirm initialization,
both inlet options, exact source count, pending queue, binding and initial
progress from logs/manifests without opening running H5, then leave the long
jobs to finish until the user requests analysis.

Explicitly authorized post-close-checker repairs can be applied while the
solver is running without changing its binary, card, launch or diagnostics.
Archive the exact old checker and keep both original frozen manifests. The
fixed checker accepts only its own one-file replacement documented by
postcheck-amendments/source-id-rule.json: old/new SHA256, original-manifest
SHA256, archived original, case path and validator-only role must all match.
Any other frozen input change still fails. check_id_repair.py exercises both
the identity failures and this narrow guard; it is not a physics test.

For an already closed run, recheck_closed.sbatch runs the full fixed checker
on a compute node and writes completion-check.corrected.json via --output.
It does not overwrite the original failed check, batch EXIT_CODE or Slurm
history, and never launches PIC. Completion/status indexes may point to the
new result while retaining the historical post-check failure explanation.

Two-node undulator replacement: prepare_two_nodes.py freezes only the stopped
4.5 m case into a new directory. run_two_nodes.sbatch allocates 32 MPI x8,
16 ranks per node, with 256 distinct physical cores. The global 4099-cell
window and all physical/diagnostic settings stay identical; the minimum local
z slab becomes 128 cells. Completion checks use the frozen rank count, while
the existing single-node free-drift copied scripts remain untouched.
