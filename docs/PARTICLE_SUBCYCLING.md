# Particle subcycling

Particle subcycling separates analytical-device orbit accuracy from the
Maxwell time step without changing the field solver or its storage. It is
intended for long or high-harmonic undulator studies where the analytical
laboratory magnetic field should be sampled more finely than the E/B grid is
advanced.

## Automatic count

For the fastest loaded particle and shortest planar-undulator period, startup
computes the number of Maxwell steps per laboratory period. It then selects
the smallest integer `N_sub` satisfying

```text
N_sub * Maxwell_steps_per_period >= particle_steps_per_undulator_period.
```

The chosen count and realized samples per period are written to the ordinary
rank-zero log. `mesh.maximum_particle_substeps` limits unexpected cost. If the
required count exceeds that limit, initialization stops before field
allocation and reports two alternatives: raise the explicit limit or refine
the Maxwell step (including a maximum Cowan-z `cell_size` z suggestion).

## What one substep does

At the start of a Maxwell step, the code samples staggered grid E/B once at
the particle position. That grid sample is held fixed while a sequence of
relativistic Boris pushes repeatedly evaluates the analytical laboratory
devices at the current substep position and time. This choice keeps particle
ownership and MPI exchange at one communication round per Maxwell step.

After the last substep, the existing charge-conserving depositor deposits the
single chord from the field-step start to the final position. Detector and
trajectory cadence also remain tied to the Maxwell step. No substep arrays,
files, or extra MPI exchanges are created.

## Numerical limits

Subcycling improves integration of prescribed undulator and external-device
forces. It deliberately does **not**:

- increase Maxwell temporal bandwidth or a field detector's Nyquist energy;
- resample the self-consistent grid E/B within the field step;
- expose intra-step current curvature to the Maxwell solver;
- replace field-grid, macro-particle, detector-aperture or CPML convergence.

Consequently a high value can converge an analytical-device orbit while the
radiation field remains under-resolved. Production cards must refine the
Maxwell grid/time step against the highest radiation band independently. A
substep convergence scan should hold the Maxwell grid fixed first; the field
grid must then receive its own convergence scan.
