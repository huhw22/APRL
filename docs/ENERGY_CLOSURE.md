# Particle/field energy-closure diagnostic

`postprocess/energy_closure` compares the energy lost by macro-particles
between two fixed laboratory particle planes with the forward radiation energy
integrated by `field_plane_analysis`. It is a read-only, standalone C++ tool;
it adds no simulation communication, memory, or output when not run.

The report is explicitly a laboratory-frame observer. It does not mix the
particle planes with the boosted-frame runtime ledger. Format version 2 also
reports laboratory mean gamma and rms energy spread at both planes, the
zero-magnet control exchange, the signal-minus-control device-associated
particle change, and the part not explained by the collected forward band.
See `LAB_FRAME_ENERGY_DIAGNOSTICS.md` for the interpretation and the fast
50 A analytical scale tool.

## Recommended experiment

Use an entrance plane before the first magnetic interaction boundary and an
exit plane after its final fringe. The field detector must retain the complete
pulse. Run the physical case and a matched zero-radiation case with identical
mesh, bunch, self-field, boundary, detector timing, and random seed. For an
undulator-only test, setting `strength_parameter: 0` supplies that baseline.

Clean each field plane with its own ballistic-reference reconstruction, then
pass the cleaned signal and cleaned baseline to `field_plane_analysis` through
`zero_radiation_baseline`. Finally provide both pairs of particle planes and
that analysis file to this tool. See
`postprocess/energy_closure/example.yaml`.

For a small-particle validation, `field_analysis` may instead name a complete
`trajectory_radiation` output containing `/far_field`. The report records
whether its scalar radiation energy came from the field plane or the
independent trajectory route. The latter is a diagnostic cross-check, not a
replacement for the production field detector.

Build and run:

```bash
cmake -S postprocess/energy_closure -B build-energy-closure
cmake --build build-energy-closure -j
./build-energy-closure/energy_closure closure.yaml
```

## Stable high-gamma calculation

The tool joins downstream crossings by `particle_id`. Duplicate downstream
crossings are rejected. With `require_all_particles: true`, a missing entrance
or exit record is also a hard error. Macro-particle `mass_kg` already contains
the represented-particle scaling; the diagnostic `weight` is not applied
again.

For each particle it avoids subtracting two rounded gamma values:

```text
gamma_in - gamma_out
  = (|u_in|^2 - |u_out|^2) / (gamma_in + gamma_out),
gamma = sqrt(1 + |u|^2).
```

It then evaluates `mass_kg*c^2*delta_gamma` in extended precision and uses a
compensated sum. If a matched baseline is present, signal minus baseline is
also formed per particle before the final sum. This is materially safer than
subtracting two total beam energies near gamma 1000.

## Meaning of closure

With a static magnetic device and no injected laser, the prescribed magnetic
field does no laboratory-frame work. After matched numerical-baseline removal,
particle kinetic-energy loss should therefore be comparable to emitted
electromagnetic energy. A single field plane measures only its resolved
forward band and aperture, so it should normally be no larger than the
particle loss.

The matched `K=0` particle change is reported as a space-charge plus numerical
control. It is not renamed as stored near-field energy: an exact near/radiative
split needs a closed lab surface or equal-lab-time 3D fields. The remaining
`device_associated_particle_loss - collected_forward_radiation` is therefore
reported as unresolved rather than assigned to one mechanism.

The comparison is not a universal exact identity. The residual also contains
radiation through transverse/backward boundaries, frequencies outside the
chosen band, finite time/aperture loss, differences in initial/final bound
self-field energy, CPML absorption before the detector, and discretization
error. The ratio must be converged with particle count, mesh, aperture, time
window, detector distance, subcycling, CPML padding, and field-reconstruction
padding before it is used as a production physics claim.

For a dense bunch, matched `K=0` subtraction does not guarantee that the
signal and baseline have the same bound/self-field energy at the exit plane.
If the baseline kinetic change itself is comparable to or larger than the
putative radiation loss, report it rather than interpreting the difference as
radiation. The runtime `energy_ledger` now performs the next-level
boosted-frame audit

```text
particle kinetic change
+ domain electromagnetic-energy change
+ flux through every physical/CPML boundary
- prescribed-source work = 0
```

on equal-box-time slices. See `ENERGY_LEDGER.md`. Laboratory particle-plane
energy and laboratory detector radiation remain a separate equal-location
comparison: relativity of simultaneity prevents directly adding those lab
particle terms to boosted-frame stored-field energy. This standalone tool
therefore intentionally does not manufacture missing global terms from one
downstream plane.
