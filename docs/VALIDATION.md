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
warns that a Gauss-consistent initial self-field has not yet been implemented.

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
