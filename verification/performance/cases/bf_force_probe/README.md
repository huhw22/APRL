# Actual converged-interval force diagnostic

This optional observer records selected complete 16-particle phase families in
the conservative-source pusher. It reads the committed states and actual used
grid fields after convergence; it does not alter trial or replay arithmetic.

Enable with both environment variables:

```
APRL_PARTICLE_FORCE_PROBE=/absolute/existing/output/directory
APRL_PARTICLE_FORCE_PROBE_IDS=/absolute/selected-source-ids.txt
```

The ID file contains positive, unique source IDs and must include every member
of each selected family `(sourceId-1)/16`. Select families from the input before
examining outcomes. With the variables absent, no writers or ID lookup tables
are created; one cached disabled branch remains per committed interval.

Each rank/thread writes its own 288-byte interval records. At normal simulation
stop, each rank closes its files and writes a completion marker and manifest
with an FNV-1a checksum, counts, layout and selected IDs. This close adds no MPI
collective. Interrupted runs are not valid completed analyses.

For an interval with constant actual grid E/B, the independent impulse and work
budgets are

\[
\Delta\mathbf u_{E,g}=\frac{q\Delta t}{mc}\mathbf E_g,\qquad
\Delta\mathbf u_{B,g}=\frac{q}{mc}\Delta\mathbf x\times\mathbf B_g,\qquad
\Delta\gamma_g=\frac{q}{mc^2}\mathbf E_g\cdot\Delta\mathbf x.
\]

Prescribed-device contributions are independently summed at the actual saved
substep events. Their electric work is checked against the pusher's saved work.
Lab energy is reconstructed through
`delta_gamma_lab = gamma_b * (delta_gamma_box + beta_b * delta_u_z_box)`.
Lab grid work is also split directly into longitudinal and transverse electric
work. A static lab magnetic field must do approximately zero lab work.

Run only after the PIC driver's final normal-stop checks:

```
python analyze.py CASE NEW_OUTPUT --geometry CAMPAIGN/GEOMETRY.json --expected-ranks 16
python analyze.py --self-test NEW_SELFTEST_DIRECTORY
```

The reader requires `PIC_COMPLETE.json`, a closed energy ledger, every rank's
completion marker and exact byte checksums. It checks actual impulse/work
identities and continuity through rank/thread migration. It refuses to overwrite
outputs. The central cohort in `run-plan.json` excludes any front startup samples.

Spatial profiles interpolate full-interval endpoint chords. Their crossing-time
error is compared with the original actual-substep particle planes. Partial
interval device contributions use linear apportioning. These profiles do not
replace the actual-substep crossing records. Harmonic work/phase association
identifies candidate driving channels; it is not a force-removal counterfactual
or proof that a frequency band is pure radiation.

`checks.cpp` tests the mean-velocity identities, static lab magnet work, invalid
selection/work rejection and binary I/O under 1/2/4 MPI ranks. The reader's two
independent fixtures and three failure controls test numerical accounting and
closed-file handling. Full-program default and switch comparisons must also
check unchanged scientific output; those tests establish observer integrity,
not FEL physical acceptance.
