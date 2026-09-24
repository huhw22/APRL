# Initial particle self-field

## Purpose

The Maxwell state must satisfy the discrete Gauss constraint for the loaded
particle charge before the first physical step. Starting from zero E with a
nonzero bunch creates a numerical startup pulse that is not particle
radiation. The default initialization therefore constructs the longitudinal
electric field once and then discards every potential/solver work array. The
evolution state remains direct SI E/B; A/phi never enters the time advance.

## Discrete construction

1. Each boosted-frame particle deposits charge on Yee vertices with the same
   first-order CIC shape used by the charge-conserving trajectory depositor.
2. Contributions on duplicated MPI z-interface vertices are summed. Global z
   planes have one owner, so dot products never count a shared plane twice.
3. A matrix-free distributed conjugate-gradient solve computes the scalar
   potential for

   ```text
   - discrete_laplacian(phi) = rho / epsilon0
   ```

   on interior vertices. The six outer potential surfaces are zero, matching
   the finite-domain field boundary.
4. E is formed by the matching Yee edge gradient. A second CIC deposition and
   a true cross-rank divergence calculation verify
   `div(E) = rho/epsilon0` before the run can continue.
   The solver also verifies that the complete CIC charge remains on interior
   vertices; a cloud touching a zero-potential surface is rejected with an
   instruction to add at least one cell of mesh/CPML padding.

The scalar potential is only a temporary numerical device for satisfying the
constraint. Four vertex slabs are live during CG (`phi`, residual, direction,
and operator image); the charge lattice is allocated only during deposition
and post-check phases. Startup reports iterations, relative L2 residual,
maximum absolute residual, total charge, and temporary memory per rank.

## Input

```yaml
initial_self_field:
  enabled: true
  relative_tolerance: 1.0e-10
  maximum_iterations: 10000
```

The block is optional and defaults to these values. An enabled solve is a hard
preflight: non-convergence or a failed post-check aborts with a correction
suggestion. Disabling it is allowed only for controlled zero-field/regression
comparisons and prints a prominent warning.

## Physical scope

This step supplies the Gauss-constrained, curl-free electric component. It is
not an exact Lienard-Wiechert initialization for a bunch with a large residual
mean velocity or velocity spread, and it does not synthesize an initial
current-dependent magnetic field. The boosted frame should therefore remain
close to the bunch mean rest frame.

Zero potential at the outer box is numerically consistent with the finite
domain, but it is not an infinite-space boundary condition. Quantitative runs
must converge transverse/longitudinal padding and CPML separation. This is
especially important for space-charge-sensitive tracking: a small Gauss
residual proves the discrete constraint, not that the chosen box is physically
large enough.
