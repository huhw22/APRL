# Laboratory-frame energy diagnostics

This is the physical interpretation layer for an undulator-only run. It is
deliberately separate from the boosted-frame runtime ledger. A prescribed,
static laboratory magnetic device has

```text
W_magnet,lab = 0.
```

Its transformed electric field can do nonzero work in the computational
frame because energy and longitudinal momentum mix under a Lorentz boost.
That boosted work is necessary for the numerical ledger, but it is not an
extra laboratory energy source.

For a longitudinal boost,

```text
Delta E' = gamma_b (Delta E - beta_b c Delta p_z).
```

A static laboratory magnet has `Delta E_magnet=0` but can exchange momentum,
so its boosted-frame energy component need not vanish. The tools below keep
the two observers separate instead of trying to cancel that mixing by
subtracting unrelated scalar energies.

## What can be separated at low cost

Two fixed laboratory particle planes give a numerically stable particle
kinetic-energy change. A laboratory field plane gives the radiation contained
by its forward aperture, time window, and frequency band. A matched `K=0`
run gives a useful control for collective self-field evolution and numerical
drift:

```text
L_signal  = particle kinetic loss in the undulator run
L_control = particle kinetic loss in the matched K=0 run
L_device  = L_signal - L_control
E_forward = collected forward radiation in the selected band
R_lab     = L_device - E_forward
```

`L_control` is a **space-charge plus numerical control**, not a direct
measurement of stored near-field energy. `R_lab` can contain differential
bound-field energy, side/backward radiation, missed aperture or frequency,
and numerical error. The tool never silently labels all of it as near field.

An exact split would additionally need either a closed laboratory radiation
surface or equal-laboratory-time three-dimensional E/B data. The latter is
intentionally not a default output because it is expensive in a boosted-frame
simulation.

The `postprocess/energy_closure` report is therefore a laboratory observer:
it reads only laboratory particle planes and laboratory radiation products.
It does not combine them with equal-boost-time stored-field energy. Report
format version 2 adds:

- entry and exit mean gamma and rms energy spread, evaluated with represented
  macro-particle mass weights;
- signal, zero-magnet control, and their per-particle matched difference;
- collected forward radiation and the unresolved laboratory remainder;
- explicit flags that near-field change was not directly measured and that
  energy spread is not a separate energy reservoir.

Energy spread measures redistribution within particle kinetic energy. The
conservation equation uses the sum of particle energies; `sigma_gamma` is a
beam-quality diagnostic and is not added as another term.

## Fast analytical estimate

`lab_frame_energy_estimate` supplies a pre-run scale estimate without loading
particles or fields:

```bash
cmake -S . -B build
cmake --build build -j --target lab_frame_energy_estimate
./build/lab_frame_energy_estimate config/lab_frame_energy_estimate.yaml
```

For a long round Gaussian beam it uses

```text
lambda_q = I / (beta c)

U_bound / L_bunch ~= (1 + beta^2) lambda_q^2 / (4 pi epsilon0)
                      * [ln(b / (2 sigma)) + EulerGamma / 2].
```

The factor `1 + beta^2` includes both laboratory electric and magnetic
self-field energy. The absolute value needs the user-supplied outer scale
`b`; this may represent a pipe radius or a finite-bunch cutoff. For a uniform
fractional transverse-size change, the long-beam approximation gives

```text
Delta U_bound / L_bunch
  ~= -(1 + beta^2) lambda_q^2 / (4 pi epsilon0)
      * ln(sigma_final / sigma_initial),
```

so the leading change is independent of the arbitrary outer cutoff. Positive
configured size change means expansion and therefore negative stored-field
change. The bound-field scale is proportional to `I^2`, whereas generalized
perveance is proportional to `I/gamma^3`. The configured radius sweep holds
current fixed; longitudinal compression or current-profile evolution must be
evaluated separately.

The undulator reference uses the angle-integrated single-electron result

```text
dE / dz = r_e m_e c^2 gamma^2 K^2 k_u^2 / 3
```

and the on-axis fundamental

```text
lambda_1 = lambda_u (1 + K^2 / 2) / (2 gamma^2).
```

The default coherent enhancement is one. This is an incoherent reference,
not a prediction for a pre-bunched beam. Current, gamma, and rms size do not
determine the longitudinal form factor; a coherent enhancement must be
provided manually from a separately justified bunching calculation.

## Requested 50 A scale test

The committed example uses `I=50 A`, `sigma_x=sigma_y=60 um`, `gamma=1174`,
`lambda_u=2.34 cm`, `K=1.093`, and `b=10 sigma`. Because no bunch duration or
device length was specified, the fundamental results are normalized per metre
of flat-top bunch and per undulator period. A 100-period device is included as
an explicit illustration.

| quantity | estimate |
|---|---:|
| line charge | `1.66782e-7 C/m` |
| electrons per metre of bunch | `1.04097e12 1/m` |
| particle kinetic energy per metre of bunch | `99.9695 J/m` |
| generalized perveance | `3.62574e-12` |
| fundamental wavelength / photon energy | `13.5595 nm / 91.4373 eV` |
| incoherent loss per electron per period | `1.33337 eV` |
| incoherent radiation per metre of bunch per period | `2.22382e-7 J/m` |
| approximate stored bound E+B energy per metre of bunch | `9.49023e-4 J/m` |
| bound-field store / particle kinetic energy | `9.49313e-6` |
| 100-period incoherent radiation per metre of bunch | `2.22382e-5 J/m` |

The absolute bound-field store is about 42.7 times the 100-period incoherent
radiation reference, but that ratio is **not** the energy exchange. Only its
change matters. In the same model:

| rms-size change | bound-field change per bunch metre | magnitude / 100-period radiation |
|---|---:|---:|
| `+1%` expansion | `-4.97517e-6 J/m` | `0.2237` |
| `+10%` expansion | `-4.76551e-5 J/m` | `2.1429` |

Thus the average envelope space-charge force is strongly relativistically
suppressed, but a percent-level change of the large bound-field store can
still be relevant to an energy-loss diagnostic. For a coherently pre-bunched
case the radiation may be much larger than the incoherent reference, and the
simulation-derived lab budget is the quantity to trust.

## Existing-output regression

The extended laboratory report was also run on the previously retained
gamma-1000 validation outputs. In the weak-collective one-electron-equivalent
case, the zero-magnet control was zero at reported precision, particle loss
was `5.93986e-21 J`, and the 50--100 eV, +/-2 mrad trajectory radiation was
`3.25979e-21 J`, or 54.88%. The unresolved 45.12% is expected to include the
finite band/angle acceptance and cannot be assigned entirely to near field.

In the `10^6` represented-electron case, the signal particles gained
`1.82555 nJ` while the zero-magnet control gained `65.2575 nJ`. Their matched
difference was a nominal `63.4320 nJ` device-associated loss, but the
0--300 eV forward plane contained only `3.63257 nJ` (5.73%). The absolute
control exchange was 1.029 times the matched difference and 94.27% remained
unresolved. The report therefore correctly identifies this old compact run as
collective-field/numerical-control dominated rather than claiming that its
matched particle difference is radiation.

## Recommended production interpretation

1. Place identical laboratory particle planes before and after the magnetic
   interaction region in the signal and `K=0` control runs.
2. Use identical charge distribution, initial self-field, mesh, time window,
   detector aperture, and random seed.
3. Subtract the field baseline at E/B amplitude level before forming the
   forward radiation energy.
4. Run `energy_closure` and inspect `L_control` before interpreting
   `L_device`. If the control is comparable to or larger than the device term,
   the result is collective-field dominated and requires convergence.
5. Treat `R_lab` as unresolved until aperture, band, time, particle count,
   grid, initial-field padding, and CPML have converged.
