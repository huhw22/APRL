# First-segment critical-mode control

This harness compares the same candidate program with `conservative-source/linear` and `conservative-source/linear-notch`, with independent actual-substep particle-plane records enabled in both. It does not change Cowan coefficients, grid spacing or the z-matched timestep.

The SPARC input is the unchanged central 231 original slices (577:808), complete 16-phase clone families, without reweighting, velocity changes, source-ID changes or arrival-phase changes. It contains 236,544 records, about 19.982 pC, spanning three first-segment slippages. Only the original first 77-period undulator remains. The crop has hard edges; it is a controlled pair, not an identical copy of the old full bunch. Report the whole crop and central 77 slices separately. If the unfiltered crop does not reproduce fast startup, it cannot settle the original full-bunch question.

Numerical settings remain boost gamma 168.80411067346705, 80×80 transverse mesh, 40 µm transverse spacing, 8.293564128124092 µm z spacing, `c dt = dz`, original CPML, 64/128 guards, strict ledger enabled. The first-segment active window is 1742 z cells; the derived source/reference mesh is 8448 cells. Both modes explicitly use the same 8-cell particle collar. Full controls require about 3879 Maxwell steps. The crop has a different initial grid phase from the old full bunch; moving the mesh center cannot correct that because the box translation moves with it. The pair itself has identical geometry and phase.

Prepare an independent campaign:

```bash
python prepare.py WORKSPACE NEW_CAMPAIGN \
  --binary WORKSPACE/build/energy-critical-control-20261105/aprl \
  --build-record WORKSPACE/results/energy_critical_control_20261105/build-r1
```

Preparation freezes input, cards, scripts and run plans in separate `off` and `on` directories. Candidate program identity is frozen at launch from the accepted build record. Full controls have a 3300-second case bound and are gated on an accepted same-binary smoke record plus explicit `--allow-full`. No node is pinned.

The old campaign-r3 attempted a 600-step reference-center stop upstream of the last magnet, which is forbidden by the real YAML validator. Job19323 therefore failed before PIC; dependent19324 was cancelled. Valid campaign-r4 keeps the full physical geometry and uses the early progress of the full off run for cost estimation. No duration-limited incomplete run or core diagnostic stop was introduced.

Run one case in a one-node, 16-MPI×8 allocation with the prepared `scripts/run_case.sbatch`. Pass `--allow-full` only after accepted build and smoke validation. Inspect early off progress before choosing how to schedule on. Two full controls may run on different nodes with the same layout; such timings are not a strict paired speed comparison. An afterok dependency guarantees ordering, not the same node.

Three raw field planes (`entry`, `und-entry`, `p1`) retain every original Maxwell sample for mode checks. Particle planes add midpoint z=1.088 m, with separate legacy and `particle_substep_plane` files. Energy planes remain entry/und-entry/p1. `analyze.py CAMPAIGN` reads only completed full controls, emits full/core BF tables and a four-panel plot. Substep four-momenta are endpoint-affine, not an on-shell snapshot; the analyzer uses their times for BF and never reconstructs gamma from the stored interpolated momentum.

`mode_review.py CAMPAIGN` additionally requires both final `PIC_COMPLETE.json` markers before opening HDF5. It compares selected raw E/B samples with the staggered cutoff signature and calculates the native-energy-weighted old/new field overlap from every strict-ledger step. Projection compatibility is not a unique mode energy. `postprocess.py` can run both reviews from a separately frozen postprocessing directory after both Slurm jobs terminate; missing or unaccepted cases are recorded without opening their HDF5. PNGs from real controls still require manual review.

The moving-window CSV follows `energy_ledger.sample_interval_steps`, then always appends a final row. Geometry validation reconstructs the initial origin as `recorded_origin + recorded_step * dz`, checks every sampled step and origin, and accepts an identical duplicated final sample. It does not require the first row to be step one. The original 19327 review stopped on that invalid assumption; `postprocess-r2` resumes the mode-only analysis with corrected validation, keeping the original failure, completed science data and existing BF plot.

`prepare_smoke.py` builds four original 31,675-particle / 257-step free-drift cards: old program, candidate default, candidate substep output, and candidate notch with substep output. It verifies default science and the diagnostic switch bitwise, as well as notch execution and complete substep cohorts. This is software/diagnostic acceptance, not a SPARC radiation or FEL physics test.

The first smoke allocation (19322) passed all four cases and both 23-file legacy science comparisons. In its single same-node sample, steps 40–240 took 0.27765 s/step for candidate default, 0.27904 for the substep diagnostic, and 0.43406 for notch plus substeps. Thus the diagnostic added about 0.5%, while the prototype notch added about 56% relative to substeps alone. These are one small free-drift sample, not production timing estimates or an accepted speed target. The full SPARC pair measures its own cost; if the physical control is useful, [PERFORMANCE_NOTES.md](PERFORMANCE_NOTES.md) identifies mathematically equivalent kernel reductions to investigate next.

There is no automatic physics acceptance. A source-filter switch is not proof that the critical branch was removed. Check its actual E/B fingerprint and native step-flip measure, then compare BF with legacy and substep times. The decisive cases are:

- Only legacy BF grows: plane-time reconstruction is implicated.
- Both time definitions grow without suppression, but true-substep BF drops with measured mode suppression: the mode has a causal role in this controlled crop.
- The mode is suppressed but true-substep BF remains: another low-order driver must be identified.
- The unfiltered crop fails to reproduce startup: repeat with the original full bunch and first segment; do not infer the original cause from this crop.

The scripts never submit jobs themselves or replace default programs, physical cards or existing results.

`power_review.py` and its bounded batch wrapper compare the completed exit fields using an unchanged snapshot of the earlier SPARC `analyze_explicit.py::pic_plane()` reducer. They retain rectangular/energy-normalized Hann windows, the fixed ±5% resonance band, signed native E×B and transverse E² proxy, physical aperture and all stored magnetic time staggering. Only the independent comparison to the online energy plane converts the installed SciPy mu0 to the explicitly compiled PIC constant; the reported older FFT convention is preserved. `plot_power_review.py` renders the completed numerical output. These are closed-only analyses, never new PIC runs or pure-radiation acceptance.
