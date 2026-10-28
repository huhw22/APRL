# Full-input moving-window pair

This campaign is explicitly authorized for the original 41,732,775-particle input.
It repeats the physical cases of jobs 18932 and 18938: the original undulator with
a raw laboratory E/B plane at 4.5 m, and no prescribed magnet with the plane at
0.5 m. Both retain the accepted fast blocked Cowan reference and moving field
window. No simulation snapshot, particle-plane diagnostic, extra absorption layer,
kernel change, or default-binary replacement is introduced.

Each case requests one exclusive 128-CPU Slurm node and uses 16 MPI ranks × 8
threads. No node is pinned. This is a feasible layout, not a measured optimum.
The actual active longitudinal grids are 4099 and 666 cells. The original
64/128-cell guards, boost, fine cell sizes, particles, release states, magnet
parameters, stepping policy and raw field cadence are retained. The held flags
still encode the source release events; they do not select old long-domain held
field pre-evolution in this enabled moving-window implementation.

`geometry.py` derives a finite disk-backed scalar reference domain that covers
every inlet query. The duration cap uses the parent's complete-input ballistic
bound and its existing factor of 1.1. The normal stop remains all valid particles
past the last detector. No extra femtosecond tail is added. The reference grid is
not a full resident E/B grid. `--verify-plan` rescans and hashes the original input
inside Slurm, including an initial-owner histogram for legacy MPI integer counts.

Preparation and execution:

1. Run `prepare.py --workspace <workspace> --root <new-campaign>` to freeze the
   cards, tested binary identity, scripts and source provenance. Existing campaign
   plans cannot be overwritten.
2. Submit `scripts/preflight.sbatch <new-campaign>` to an available node. It checks
   the full input geometry and runs both bounded short cases in the exact 16×8
   layout. Fixed cross-layout field/particle gates and complete diagnostics must
   pass before `preflight/PASS.json` is written.
3. Submit each case's `run.sbatch <case-directory>` after reviewing the preflight.
   Each job verifies its frozen inputs and binding before starting APRL. The
   unchanged local 600/60-second shutdown protection is retained. The 144-hour
   allocation limit is a hard cap, not an estimated runtime.
4. Confirm initialization and normal stepping, record job IDs and log paths,
   then stop monitoring until requested. Completion checking runs inside the
   existing job allocation; spectra are not computed automatically.

`check_completion.py` checks normal physical termination and diagnostic integrity.
Its PASS is not a claim of pure radiation, near-field separation or full-scale
physical equivalence. The first full-scale moving-window pair is intended to
provide the field data for the subsequent comparison. All pending particles are
still held on the inlet rank and scanned; particle streaming is outside this run.
