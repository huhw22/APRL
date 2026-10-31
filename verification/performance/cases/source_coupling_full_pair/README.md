# Conservative-source full-input pair

Replay the latest 19099 undulator 4.5 m and 19098 free-drift 0.5 m cards with
the accepted source-coupling/longitudinal-handoff binary 88a5a8d. The former
uses two nodes, 32 MPI x 8; the latter one node, 16 MPI x 8. No node pinning.

Keep all 41,732,775 original input events, original reference/active meshes,
64/128 guards, fine dz/dt, fast reference generation, input release events,
undulator, raw E/B sampling, physical stop, entry/exit particle and independent
energy planes, v5 energy ledger and strict discrete operator audit. No
retirement layer or checkpoint. The only physics selector changed is
mesh.particle_coupling=conservative-source. No kernel or binary rebuild.

prepare.py checks the exact binary hash, compiled-source closure, existing
78-required and multi-node acceptance, freshly hashes the original input,
whitelists every card change, and freezes scripts/cards/evidence. It reuses
existing 16-total-rank cross-node tests; the actual 32-rank production layout
is checked at startup rather than represented as a new small scientific test.

Run preparation from the workspace:

```bash
python3 source/aprl/verification/performance/cases/source_coupling_full_pair/prepare.py \
  /absolute/workspace /absolute/new/campaign
```

Submit the copied per-case run.sbatch with --nodes=2 --ntasks=32 for the
undulator and --nodes=1 --ntasks=16 for free drift. Both use exclusive nodes,
16 ranks per node, 8 threads, 144-hour existing wall limit and no requeue.
The original self/vader/tcp transport and 600/60-second shutdown watchdog
remain; all short diagnostic environment switches are explicitly cleared.

The runner distinguishes the simulation exit code from post-close checks.
The reused checker includes the repaired particle_id/input-row vs arbitrary
source_id rule. An additional completed-source-log check verifies all rank
steps and fixed-point tolerances. These are execution checks, not proof of
physical radiation accuracy. Never inspect running HDF5 outputs.
