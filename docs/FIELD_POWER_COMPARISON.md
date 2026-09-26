# Particle retirement and matched power comparison

## Purpose and limits

Particle retirement is a diagnostic route for obtaining a cleaner downstream
field plane without modifying the E/B propagation kernel. A physical particle
ends at a fixed laboratory entrance after all magnetic interaction regions. A
small ballistic carrier then deposits its remaining straight-line current with
a C2 quintic taper and disappears before the detector.

This is not an exact charge-conserving removal of a net charge. Smoothness
controls the spectrum of the numerical impulse; it does not prove that the
impulse is negligible. The route is acceptable only when its zero-radiation
power is small and converges downward relative to the wanted radiation.

The first baseline is a second simulation with identical mesh, CPML, bunch,
random seed, retirement entrance/length, detector, stop, and sample times, but
with the radiating magnetic strength set to zero. The comparison tool subtracts
this baseline at the E/B-amplitude level. This measures straight-beam
near field plus retirement/initialization/boundary artifacts. It is a screening
baseline, not a mathematically exact counterfactual when the magnet changes the
particle distribution at the retirement entrance. Entrance-state replay is
not implemented in version 1; such a baseline would isolate the retirement
artifact more exactly.

## Build and run

```bash
cmake -S postprocess/field_power_compare -B build-field-power-compare
cmake --build build-field-power-compare -j
./build-field-power-compare/field_power_compare \
  postprocess/field_power_compare/example.yaml
```

The input card names the signal field plane, the matched zero-radiation field
plane, an optional cleaned result from `field_reconstruction`, an optional
trajectory-far-field file, the photon-energy band, a small transverse read
batch, and the output HDF5 file. The complete example is in
`postprocess/field_power_compare/example.yaml`.

All card mappings reject unknown keys. Existing output is refused unless
`output.overwrite: true` is explicitly selected; see
`docs/RUN_PROVENANCE.md`.

The two field files must have identical plane geometry and identical stored
laboratory sample times. Actual detector times need not be perfectly uniform:
the tool retains them for broadband power and linearly resamples each field
component to one uniform grid before its finite-record DFT.

Memory is bounded by the time series for a small transverse batch plus four
filtered components for one point. The complete detector plane is never copied
into memory. This is an analysis-server tool and adds no simulation-time MPI or
I/O path.

When `analytic_electron_reconstruction` is present, it must have the same
plane geometry, sample count, and laboratory-time extent. Its uniform time axis
may differ from the raw detector's nonuniform interior samples. The tool
resamples both routes to the same uniform band-analysis axis and computes the
field-level residual between the retirement estimate and the old analytic
straight-line-electron subtraction estimate.

## Definitions

For signal `s` and zero-radiation baseline `0`, the radiation estimate is
formed before any quadratic quantity:

```text
E_delta = E_s - E_0
B_delta = B_s - B_0
P_delta(t) = integral max(((E_delta cross B_delta)_z)/mu0, 0) dx dy.
```

The tool also writes signed Poynting power. It never uses
`P_signal - P_baseline`, because that would discard interference terms and is
not the power of the difference field.

The requested photon-energy band is applied to every transverse point by
zeroing out-of-band DFT bins and inverse transforming E and B. The output
contains broadband and band-limited signed/forward power curves, trapezoidal
time-integrated energies, and peak forward powers for signal, baseline, and
difference fields.

There are two band-power definitions. The legacy instantaneous curve uses the
real filtered fields directly. The strict trajectory comparison instead forms
the Hilbert quadratures after band filtering and evaluates

```text
P_cycle(t) = integral 0.5 Re((E + i H[E]) cross
                             conjugate(B + i H[B]))_z / mu0 dx dy.
```

This cycle-averaged analytic-signal power has the same convention as
`trajectory_band_power_W`. In particular, an instantaneous real-field peak can
be almost twice its cycle average even when both calculations are correct. The
legacy `field_to_trajectory_*` attributes remain for compatibility; validation
must use the `cycle_averaged_*_to_trajectory_*` attributes.

If `trajectory_far_field` is supplied, it must have exactly one shot, linear
frequency spacing, the same photon-energy limits, and at least two points on
both angular axes. The tool reconstructs the cycle-averaged analytic-signal
power from `electric_field_spectral`, integrates it over the configured solid
angle, centres the periodic time origin on the peak, and normalizes its time
integral to `band_energy_J`. It then reports field/trajectory ratios for both
peak power and band energy.

The optional analytic reconstruction is kept as an independent route. The
tool reports its broadband and band-limited power, its ratios to the retirement
difference and trajectory result, and a `method_residual` formed as

```text
E_residual = E_retirement_difference - E_analytic_reconstruction
B_residual = B_retirement_difference - B_analytic_reconstruction.
```

This distinguishes disagreement of the cleaned fields from a misleading
difference between two separately squared powers.

The trajectory angular grid must represent the same radiation acceptance as
the finite field plane. That relation is not guessed by the program because a
near-field plane does not in general map to one source point and one angular
rectangle.

## HDF5 output

`/power_comparison` contains:

- `time_s` and `band_time_s`;
- `{signal,baseline,difference}_{signed,forward}_power_W`;
- `{signal,baseline,difference}_band_{signed,forward}_power_W`;
- `{signal,baseline,difference}_band_cycle_averaged_{signed,forward}_power_W`;
- corresponding forward-energy and peak-power attributes;
- `baseline_to_difference_band_energy_ratio` and
  `baseline_to_difference_band_peak_power_ratio`;
- optional `trajectory_relative_time_s`, `trajectory_band_power_W`, trajectory
  energy/peak attributes, legacy instantaneous field/trajectory ratios, and
  strict `cycle_averaged_field_to_trajectory_*` ratios;
- with an analytic reconstruction, its raw, band, and cycle-averaged power;
  band-limited `method_residual` curves and metrics; retirement/analytic ratios;
  and both instantaneous and cycle-averaged analytic/trajectory ratios;
- the detector area, plane z, photon-energy limits, Nyquist energy, definitions,
  and a completion marker.

## Acceptance sequence

No single universal percentage is imposed. For each physical problem:

1. Check that the detector Nyquist energy exceeds the requested maximum and
   that the useful pulse is fully inside the saved time window.
2. Require both zero-radiation baseline peak power and energy to be small
   relative to the difference field.
3. On a small-particle run, compare difference-field peak power and energy with
   the trajectory result over a genuinely matched angular acceptance. Use the
   cycle-averaged field quantities and repeat at increasing macro-particle count
   with fixed total charge; a sparse coherent bunch can move the integrated
   energy by order unity while giving an apparently plausible peak.
4. Run the same radiating case without retirement through
   `field_reconstruction`, then inspect the retirement/analytic ratios and the
   field-level method residual. Vary analytic transverse smoothing separately;
   this exposes sensitivity to the assumed electron core rather than radiation.
5. Repeat at two or more retirement lengths and detector gaps. The wanted
   power/energy should stabilize while baseline artifacts decrease or remain
   safely subdominant.
6. Repeat longitudinal resolution/cadence, transverse resolution, aperture,
   CPML thickness, and frequency/angular-grid convergence.

For an observation angle `theta`, transverse sampling must resolve the
transverse wavelength approximately `lambda/sin(theta)`. A z-resolved grid can
therefore agree near axis while missing wide-angle power if dx/dy are too
large. This failure must not be attributed to retirement.
