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
benchmark remains required before production use.

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
