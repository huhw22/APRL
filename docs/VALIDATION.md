# Numerical validation status

This page records what has actually been tested. It distinguishes numerical
convergence from comparison with a simplified analytical model.

## Two-plane orbit validation

An eight-macroparticle, ten-period, `gamma=2`, `K=0.1` regression run was
executed with two MPI ranks. All eight reference/field-plane pairs matched:

```text
maximum transverse straight-line error       4.26354e-14 m
maximum arrival-time error                   4.24429e-22 s
maximum relative proper-velocity change      1.53732e-7
```

The built-in HDF5 records and summary attributes agreed with an independent
C++ two-plane calculation.

## Trajectory far-field peak intensity

A separate central one-electron test used the same ten-period analytical
fringe-field undulator and one trajectory sample per completed field step.
The central-axis peak of `d^2W/(domega dOmega)` converged as follows:

| particle/field steps per period | peak energy | peak intensity (J s/sr) | relative to 160 |
|---:|---:|---:|---:|
| 40 | 7.96 eV | 1.798003131e-36 | -2.3144e-5 |
| 80 | 7.96 eV | 1.798037796e-36 | -3.8649e-6 |
| 160 | 7.96 eV | 1.798044745e-36 | reference |

Thus the earlier test established more than the peak position: for this
fundamental, low-K regression, the 40-step result is numerically converged in
peak intensity to about `2.4e-5` relative to the 160-step run.

A constant-longitudinal-velocity, small-K hand formula gives
`2.65099e-36 J s/sr`, 32% above the converged result. It is not an acceptable
absolute reference for this low-gamma orbit with the program's smooth fringe
field and exact longitudinal dynamics. The table therefore supports numerical
convergence of the implemented orbit/radiation calculation, not a 32% claim
against a physically matched analytical solution. A matched high-gamma
benchmark is reported below. It diagnoses the power-comparison path but does
not replace a production high-harmonic convergence study.

## Prescribed-device particle subcycling

An eight-particle, `gamma=4` planar-undulator test deliberately used only
`1.667701552` Maxwell steps per one-micrometre period. Requesting eight
particle samples selected five Boris substeps and realized `8.33850776`
analytical-device samples per period. A reference run reduced Cowan-z `dz`
and its time step by five, retained the same physical box and detector, and
therefore needed one particle step per Maxwell step.

At the fixed laboratory particle plane, all eight particle IDs matched. The
coarse-grid five-substep result differed from the five-times-refined reference
by at most:

```text
position norm                         5.90055e-13 m
arrival time                          2.79510e-23 s
relative proper velocity              6.33332e-5
relative gamma                         2.01306e-7
```

With the same coarse Maxwell step but no particle subdivision, those maxima
were `4.72221e-8 m`, `6.44754e-18 s`, `2.73793e-3`, and `2.39028e-5`,
respectively. This validates the intended analytical-device orbit use, not a
relaxation of Maxwell sampling: grid E/B, current deposition and detector
cadence were still evaluated at the field step.

## Field reconstruction

Two checks have been completed:

- For one analytical uniformly moving charge, the reconstructed subtraction
  leaves a 3.35% RMS electric-field residual away from the central
  one-micrometre macro-particle core on the finite test grid. The core is
  intentionally CIC-regularized and is not expected to match a singular point
  charge.
- On the eight-particle regression, using the left ballistic reference or the
  actual colocated particle detector changes the reconstructed forward energy
  by less than `2e-9` relative. Both give a raw-to-cleaned forward-energy change
  of about `-20.63%` for that deliberately small numerical aperture and time
  window.

The `-20.63%` number is a diagnostic result for that test window, not a
universal correction. Production cards must repeat padding, smoothing,
aperture, time-step, and transverse-resolution convergence studies.

## Field-plane spectrum and coherence normalization

A uniform forward plane wave was generated on an `8x8` aperture with 256
samples, exactly 16 samples per optical cycle. The analytical accepted energy
is

```text
epsilon0 c E0^2 A T / 2 = 2.174499822e-14 J.
```

The full-record angular-spectrum tool returned `2.174499822e-14 J`, a relative
difference of `5.8e-16`. The angular coherent fraction, spatial reference CSD
and one-frequency temporal CSD were all one; the fluctuation spectrum was
zero.

The same signal was divided into four nonoverlapping 64-sample Hann windows.
After the documented `sqrt(N/sum(w^2))` normalization and integration over all
positive-frequency side-lobes, the mean window energy was
`5.436153658e-15 J`, versus `5.436249555e-15 J` analytically (relative
`1.76e-5`). The residual is the finite Hann/end-bin discretization. Coherent
fraction and both normalized CSD checks remained one, while maximum numerical
fluctuation density was `2.1e-25 J/eV/sr`.

A second four-window test alternated two equal-energy, exactly orthogonal
transverse FFT modes. Both the Gram-matrix global degree of transverse
coherence and the integrated ensemble-mean-field fraction were `0.5`, as
expected. This separately checks that the global metric is not hard-wired to
the rank-one reference tests.

Two I/O regressions also passed: subtracting an identical field/baseline pair
gave exactly zero band energy, and an existing `/reconstructed_field` file was
read and analyzed without conversion. These checks validate FFT
normalization, discrete solid-angle Jacobian, time-window normalization,
amplitude subtraction and both accepted HDF5 layouts. They do not validate a
finite-aperture FEL result; detector distance, aperture and grid convergence
remain problem-specific.

## Experimental retirement power closure

A two-macroparticle CPML run used `gamma=4`, boost gamma 2, a three-period
`K=0.1` undulator, a 2 micrometre C2 retirement layer, and a 32x32 laboratory
field plane. A matched `K=0` run used the same bunch seed and all numerical
settings. The field cadence gave a 49.59 eV Nyquist limit; comparison used
25--48 eV.

| quantity | baseline | field difference | trajectory near-axis | baseline/difference | field/trajectory |
| --- | ---: | ---: | ---: | ---: | ---: |
| peak power | 3.52197e4 W | 4.12602e5 W | 3.85292e5 W | 0.0854 | 1.071 |
| band energy | 5.96082e-12 J | 7.50761e-11 J | 8.17674e-11 J | 0.0794 | 0.918 |

The trajectory reference above used `theta_x,theta_y` in +/-0.05 rad, which the
coarse transverse grid can resolve. Expanding the trajectory integral to
+/-0.35 rad changed the field/trajectory ratios to 0.359 in peak power and
0.222 in energy. This is consistent with unresolved wide-angle transverse
wave number in the test mesh and is not evidence that retirement loses 78% of
the near-axis radiation.

A separate straight-beam decomposition check compared one and two MPI ranks.
The two field files differed by `5.28e-5` in 1--20 eV forward energy and
`7.67e-4` in peak power relative to the one-rank result; both runs reported two
retirement entries and completions with the same residual exit charge.

These tests establish the power-comparison plumbing, event handling, and MPI
carrier migration. They do not yet approve retirement for production: taper
length/gap convergence, a physically representative high-gamma grid, and an
entrance-state replay baseline remain required.

## High-gamma paraxial power closure

A dedicated diagnostic used laboratory and boost `gamma=1000`, a three-period
`K=0.5`, 30 mm undulator, a coherent transverse Gaussian sheet, a 200 micrometre
scale field plane, and the 50--100 eV fundamental band. The trajectory integral
used the exact direction-dependent formula over a paraxial +/-0.26 mrad angular
rectangle. Total bunch charge was fixed while the number of macro-particles was
increased. The ratios below use the field-side cycle-averaged analytic-signal
power, so their peak definition matches `trajectory_band_power_W`.

| macro-particles | field/trajectory peak | field/trajectory band energy |
| ---: | ---: | ---: |
| 2 | 1.0581 | 1.6368 |
| 8 | 0.9973 | 1.5436 |
| 32 | 1.0085 | 1.0898 |
| 128 | 1.0096 | 1.0067 |

The matched `K=0` retirement baseline stayed at roughly `1e-17` or less of the
radiating peak and at the `1e-18` scale in band energy. The earlier
apparent peak excess was primarily a definition mismatch: for eight particles,
the instantaneous real-field peak was 22.83 MW, while its cycle average was
16.84 MW and the trajectory peak was 16.89 MW. Once identical power definitions
and sufficient macro-particle sampling were used, both peak and energy closed
to about one percent.

Changing the paraxial trajectory rectangle from +/-0.18 to +/-0.35 mrad moved
the accepted trajectory energy, as expected, but no single cutoff repaired the
sparse eight-particle energy discrepancy. At these angles, exact nonparaxial
geometric corrections are far too small to explain a 50% error. The observed
convergence with macro-particle count therefore identifies coherent sampling
statistics as the dominant energy-error source in this diagnostic; transverse
FDTD dispersion remains a separate percent-level convergence item.

## Retirement versus analytic electron subtraction

The old non-retired route was retained as an independent check: a second plane
fits each particle's downstream straight line, and `field_reconstruction`
subtracts the resulting analytic uniformly moving electron field. In the
32-particle high-gamma case all particles matched. The maximum transverse
straight-line error was 0.401 micrometres, the maximum arrival-time error was
`1.42e-18 s`, and the maximum direction error was `2.13e-6 rad`.

With no transverse core smoothing, the two cleaned field routes agreed in
cycle-averaged peak to `2.1e-5` relative. Retirement energy was 3.422 nJ and the
analytic-subtraction result was 3.320 nJ, a 3.1% difference; their field-level
residual contained 5.0% of the analytic result's instantaneous band energy.
Increasing the analytic core smoothing to 5 micrometres reduced that residual
energy ratio to 0.86%, while changing the analytic cycle-averaged energy to
3.328 nJ. The peak is robust. Across the three sampled smoothing values, the
analytic energy itself spans about 1.1%, while the two cleaning routes differ
by 2.8--3.9%; the old route therefore retains a model-dependent electron-core
parameter without attributing the full route difference to smoothing alone.

The same analytic route was then repeated with 128 macro-particles and no
additional transverse smoothing. All 128 two-plane records matched. The
maximum straight-line transverse error was 0.279 micrometres, the maximum
arrival-time error was `1.12e-18 s`, and the maximum direction error was
`1.48e-6 rad`.

| 128-particle route | cycle-averaged peak | 50--100 eV energy | peak/trajectory | energy/trajectory |
| --- | ---: | ---: | ---: | ---: |
| trajectory radiation | 20.7420 MW | 2.90963 nJ | 1.0000 | 1.0000 |
| retirement minus `K=0` | 20.9416 MW | 2.92926 nJ | 1.00962 | 1.00675 |
| analytic electron subtraction, 0 micrometre | 20.9399 MW | 2.89734 nJ | 1.00954 | 0.99578 |

At this particle count, retirement and analytic subtraction differ by only
`8.3e-5` relative in peak and 1.10% in energy. They bracket the trajectory
energy: retirement is 0.675% high and analytic subtraction is 0.422% low.
The field-level instantaneous band residual contains 1.54% of the analytic
field energy; as a quadratic residual-field metric, this is not the scalar
energy difference between the two routes.

This result establishes that the straight-line analytic route remains useful
after macro-particle convergence and does not require a modified propagation
kernel. No 128-particle smoothing scan has yet been performed, so its core
parameter sensitivity at this count remains an open convergence item.

The resulting first-version policy is therefore:

- use retirement minus the matched `K=0` field at the E/B-amplitude level as
  the primary, model-light estimate after macro-particle convergence;
- retain analytic electron subtraction as an independent diagnostic and report
  its smoothing sensitivity rather than tuning it silently;
- compare either route with trajectory radiation only through the
  cycle-averaged power definition and a documented angular acceptance.

This benchmark is deliberately a three-period coherent fundamental test, not a
production approval for high harmonics. Time step, longitudinal and transverse
resolution, retirement length/gap, aperture, CPML, and long-undulator
convergence are still required for a physical run. The current simulator also
requires box-padding convergence for its Gauss-consistent electrostatic
initial self-field.

## Particle/field energy closure

A dedicated closure test used the same `gamma=1000`, three-period `K=0.5`,
30 mm undulator and 128 transverse Gaussian macro-particles. Laboratory
particle planes bracketed the complete compact fringe support. Energy loss was
formed per particle from

```text
(|u_in|^2-|u_out|^2)/(gamma_in+gamma_out)
```

in extended precision and then accumulated with a compensated sum. This avoids
subtracting two rounded total beam energies.

The full-charge (`10^6` represented electrons) test is not a two-term energy
balance. With the Gauss-consistent initial field enabled, the signal particles
gained `1.82555 nJ`, while the matched `K=0` bunch gained `65.2575 nJ` from
its evolving collective/near field. Their per-particle difference is therefore
a nominal `63.4320 nJ` loss, but the matched, background-cleaned forward field
contains only `3.63257 nJ` over 0--300 eV. The difference includes a large
change in bound/self-field energy and cannot be labelled radiation loss. It
shows why a dense charged bunch requires a global ledger containing kinetic
energy, stored E/B energy and flux through every boundary; a single downstream
plane is insufficient.

That retained compact case did not include a converged entrance/exit field
window, so it diagnosed a failure mode but could not localize the missing
energy. A new test kept the charge and 128 macro-particles but used colocated
laboratory particle/field planes at 0 and 159 mm and integrated the **signed**
raw Poynting flux. Its window sweep was:

| bunch start / lab stop | longitudinal cells | matched raw field transport | fraction of `5.919605 nJ` particle loss |
|---|---:|---:|---:|
| `-1 mm / 175 mm` | 400 | `3.230129 nJ` | `54.57%` |
| `-50 mm / 175 mm` | 400 | `3.696551 nJ` | `62.45%` |
| `-50 mm / 175 mm` | 800 | `3.696547 nJ` | `62.45%` |
| `-50 mm / 350 mm` | 400 | `4.935941 nJ` | `83.38%` |
| `-50 mm / 540 mm` | 400 | `5.874187 nJ` | `99.23%` |

Doubling `Nz` at fixed start/stop did not add detector samples or recover the
pulse; extending the simulated lab interval did. In the longest run, the
signal particle loss was `5.902545 nJ`, the `K=0` particle loss was
`-0.017060 nJ`, and their matched difference was `5.919605 nJ`. The matched
raw entrance/exit field transport left only `0.045418 nJ` (0.77%) unresolved.
This remainder still permits transverse flux, stored-field change between the
planes, residual tail truncation and numerical error because two z planes are
not a fully closed surface.

The independent amplitude-level signal-minus-`K=0` field analysis over the
full transverse plane and 0--300 eV gave `5.430281 nJ`, or 91.73% of the
particle loss. The `0.443906 nJ` difference between raw signed transport and
that radiation product contains the bound/non-propagating component,
frequencies or angles rejected by the analysis, and field-decomposition cross
terms. It should not be called missing energy. This test identifies detector
time-window truncation as the dominant cause of the earlier large discrepancy.

A weak-collective control reduced the total represented charge to one electron
while retaining 128 fractional macro-particles and the Gauss-consistent initial
field. Its matched `K=0` kinetic change was zero at the reported precision. The
undulator case lost

```text
5.939857680e-21 J.
```

The independent trajectory far field from the same run gave
`2.858701852e-21 J` in 50--100 eV and the paraxial `+/-0.26 mrad` rectangle,
or 48.13% of the particle loss. Expanding the rectangle to `+/-2 mrad` gave
`3.259794048e-21 J`, or 54.88%. The field-plane/trajectory benchmark above is
already closed to about one percent for the `+/-0.26 mrad` acceptance. Thus
the controlled particle loss has the correct sign and order of magnitude, and
the deliberately finite forward band is smaller as expected. The residual may
contain frequencies outside 50--100 eV, non-forward flux, bound-field change
and discretization error; this test does not claim exact global conservation.

Turning off the initial field was useful only as a diagnostic and is not a
physical closure test. In that path, the undulator run gained `30.4231 nJ`
while `K=0` stayed unchanged. Reducing total charge from `10^6` electrons to
one changed the gain to approximately `3.044e-20 J`, demonstrating the
expected quadratic charge scaling of the grid-mediated term. Increasing
prescribed-device sampling from about 75 to 270 samples per period did not
remove it, because particle subcycling does not refine the Maxwell/current
cadence. Raising macro-particles from 128 to 512 reduced the gain by only about
12%, not as `1/N`. These runs confirm that disabling the Gauss-consistent field
creates a large startup/self-field artifact and must not be used to approve
energy conservation.

The reusable `postprocess/energy_closure` tool generated these reports. It can
read either a complete field-plane analysis or a complete trajectory far-field
file, and can optionally combine four reconstructed field planes into the
longitudinal laboratory control volume above. Production approval still needs
aperture, mesh, CPML, particle-count and longer-device convergence; the 99.23%
result is a small three-period conservation test, not a production FEL claim.

## Runtime full-energy ledger and energy-spread decomposition

The optional runtime ledger was exercised with the same gamma-1000,
52x52x400-cell, 128-macro, three-period configuration at total represented
charges of one electron and `10^6` electrons. Every run used the
Gauss-consistent initial field, no particle crossed the CPML entrance, and
`K=0` was paired with `K=0.5`.

For the one-electron driven control, the boosted-frame terms were
`delta K=5.5087e-24 J`, `delta U=4.8912e-24 J`, cumulative outward flux
`-7.3148e-26 J`, and prescribed-device work `9.9570e-24 J`. The remaining
`3.6970e-25 J` was 1.81% of the summed exchange scale. The near-zero `K=0`
case had a smaller absolute residual, `1.9658e-25 J`, but a 32.1%
exchange-normalized residual because its denominator is nearly zero.

For `10^6` electrons the driven terms were `delta K=1.1357e-14 J`,
`delta U=4.8669e-12 J`, outward flux `-7.3231e-14 J`, and prescribed work
`4.4430e-12 J`. The residual was `3.6204e-13 J`, or 3.85% of the exchange
scale and 7.996% of the initial boosted kinetic-plus-field energy. This is a
successful accounting/sign test, not yet a production conservation result.
The code now warns at shutdown when the latter relative residual exceeds the
configured tolerance.

The large-charge driven group ended at mean laboratory gamma `999.9196671`,
projected `sigma_gamma=0.0437625`, and linearly detrended
`sigma_gamma=0.0109717`; 93.7% of projected variance was associated with the
linear gamma-z correlation. Its matched `K=0` group ended at mean gamma
`1000.0001609` and `sigma_gamma=0.00014297`. These quantities are accelerator
phase-space diagnostics, not extra terms in the energy equation. They are
evaluated on an equal-box-time slice, so fixed laboratory particle planes are
still required for rigorous entrance/exit slice analysis.

The same smoke case completed with one and two MPI ranks while only rank zero
opened the ledger HDF5. The two-rank output contained seven committed records
and `complete=1`, and the standalone C++ reader recovered its initial and
final records. Disabled cards execute none of the new boundary-power,
statistics, reduction, buffering, or file paths.

## Input-weight, numerical-preflight, and resource smoke tests

A two-rank HDF5 v2 input with two records of relative weights 1:3 and a YAML
total of ten electrons produced represented macro-particle counts 2.5 and 7.5.
Startup reported a global total of ten and preserved the electron
charge-to-mass ratio. The existing HDF5 v1 example remains readable and uses
equal weights.

The undulator temporal-resolution guard was exercised with a deliberately
under-resolved card. It measured 1.6677 particle steps per period against a
requested 32 and rejected the run before field allocation. For that card it
reported `dt <= 3.47678e-17 s` and Cowan-z `dz <= 1.04231e-8 m` as corrective
limits. This is an input safety check, not a convergence result for a selected
harmonic.

On the two-rank 40x40x60-cell generated-Gaussian smoke test, the one-step
startup calibration predicted approximately 0.0304 seconds per step after the
1.25 time safety factor and a 0.456-second duration upper bound. The physical
loop completed 13 steps in approximately 0.298 seconds; complete wall time was
approximately 0.327 seconds. Modeled peak memory was 51.6 MiB per rank versus
roughly 41--43 MiB measured resident memory. These numbers only establish that
the report is conservative for this local smoke test. They are not transferable
performance promises; a representative target-machine calibration is still
required.

The relativistic transform was also isolated with both beam and boost
`gamma=1000`. The direct product-difference form produced approximately
`1.0004e-10` relative momentum/gamma round-trip error. After changing the
longitudinal four-velocity transform to light-front components, the same
round trip measured `2.27e-16` in both quantities. A zero-field two-rank run
then completed 22 steps and its trajectory output used the same inverse
transform. The reported individual-electron lab-energy roundoff scale was
`1.13e-7 eV`. This validates the transform/output precision path; it does not
by itself validate a small radiative energy loss, which still requires
per-particle gamma differences and compensated or extended-precision
accumulation.

## Integer mesh and unequal-slab regression

The committed mesh cards now use `cells: [40, 40, 60]` and
`cell_size: [1.0, 1.0, 0.2]`. The global physical extent is produced by
multiplication; no length/spacing quotient is rounded into a cell count.

A seven-rank generated-Gaussian run deliberately divided 60 z cells into
unequal slabs. Startup reported integer slab counts in `[8, 9]` with remainder
4, and the run reached its configured stop after 13 steps. The same migrated
card also completed the two-rank HDF5 field/particle-detector smoke test.

Additional negative tests rejected a fractional cell count, rejected the old
`lengths`/`resolution` keys with a direct migration message, and rejected a
Cowan card with `dx < dz` while recommending the required x `cell_size`.
These checks validate deterministic grid construction and input diagnostics;
they are not a claim that every MPI decomposition has identical performance.

## Gauss-consistent initial particle field

The generated eight-particle card was initialized on the 40x40x60 grid with
`relative_tolerance: 1e-10`. One MPI rank and seven unequal z slabs both
required 311 CG iterations. Their post-deposition discrete Gauss residuals
were respectively `9.497829439e-11` and `9.497828774e-11` in relative L2
norm; the maximum absolute residuals were `1.039594744` and `1.039594375`
`V/m^2`. The represented charge agreed at `-1.281741307e-18 C`.

The seven-rank decomposition retained 8--9 z cells per rank with remainder 4
and completed the same 13-step physical run. Its maximum reported temporary
Poisson memory was 0.5643 MiB per rank, compared with 3.1806 MiB for the
single-rank solve. A separate two-rank HDF5-input run converged in 276
iterations to `8.896098713e-11`, then completed both field- and
particle-detector output.

These tests validate CIC interface summation, unique ownership of shared z
vertices, matrix-free scalar halos, and the true cross-rank `Ez` divergence
post-check. They validate the discrete Gauss constraint, not the infinite-
space accuracy of the zero-potential outer boundary or an exact moving-bunch
Lienard-Wiechert initialization; box padding and boost-frame choice remain
physical convergence studies.

A two-rank negative test capped the solver at one iteration and aborted before
field advance with explicit tolerance, iteration, and aspect-ratio remedies.
