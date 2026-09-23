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
with the radiating magnetic strength set to zero. This measures straight-beam
near field plus retirement/initialization/boundary artifacts. It is a screening
baseline, not a mathematically exact counterfactual when the magnet changes the
particle distribution at the retirement entrance. A future state-replay
baseline initialized from that entrance would isolate the retirement artifact
more exactly.

## Build and run

```bash
cmake -S postprocess/field_power_compare -B build-field-power-compare
cmake --build build-field-power-compare -j
./build-field-power-compare/field_power_compare \
  postprocess/field_power_compare/example.yaml
```

The input card names the signal field plane, the matched zero-radiation field
plane, an optional trajectory-far-field file, the photon-energy band, a small
transverse read batch, and the output HDF5 file. The complete example is in
`postprocess/field_power_compare/example.yaml`.

The two field files must have identical plane geometry and identical stored
laboratory sample times. Actual detector times need not be perfectly uniform:
the tool retains them for broadband power and linearly resamples each field
component to one uniform grid before its finite-record DFT.

Memory is bounded by the time series for a small transverse batch plus four
filtered components for one point. The complete detector plane is never copied
into memory. This is an analysis-server tool and adds no simulation-time MPI or
I/O path.

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

If `trajectory_far_field` is supplied, it must have exactly one shot, linear
frequency spacing, the same photon-energy limits, and at least two points on
both angular axes. The tool reconstructs the cycle-averaged analytic-signal
power from `electric_field_spectral`, integrates it over the configured solid
angle, centres the periodic time origin on the peak, and normalizes its time
integral to `band_energy_J`. It then reports field/trajectory ratios for both
peak power and band energy.

The trajectory angular grid must represent the same radiation acceptance as
the finite field plane. That relation is not guessed by the program because a
near-field plane does not in general map to one source point and one angular
rectangle.

## HDF5 output

`/power_comparison` contains:

- `time_s` and `band_time_s`;
- `{signal,baseline,difference}_{signed,forward}_power_W`;
- `{signal,baseline,difference}_band_{signed,forward}_power_W`;
- corresponding forward-energy and peak-power attributes;
- `baseline_to_difference_band_energy_ratio` and
  `baseline_to_difference_band_peak_power_ratio`;
- optional `trajectory_relative_time_s`, `trajectory_band_power_W`, trajectory
  energy/peak attributes, and field/trajectory ratios;
- the detector area, plane z, photon-energy limits, Nyquist energy, definitions,
  and a completion marker.

## Acceptance sequence

No single universal percentage is imposed. For each physical problem:

1. Check that the detector Nyquist energy exceeds the requested maximum and
   that the useful pulse is fully inside the saved time window.
2. Require both zero-radiation baseline peak power and energy to be small
   relative to the difference field.
3. On a small-particle run, compare difference-field peak power and energy with
   the trajectory result over a genuinely matched angular acceptance.
4. Repeat at two or more retirement lengths and detector gaps. The wanted
   power/energy should stabilize while baseline artifacts decrease or remain
   safely subdominant.
5. Repeat longitudinal resolution/cadence, transverse resolution, aperture,
   CPML thickness, and frequency/angular-grid convergence.

For an observation angle `theta`, transverse sampling must resolve the
transverse wavelength approximately `lambda/sin(theta)`. A z-resolved grid can
therefore agree near axis while missing wide-angle power if dx/dy are too
large. This failure must not be attributed to retirement.
