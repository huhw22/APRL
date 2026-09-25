# Initial particle self-field

## Purpose

The Maxwell state must satisfy the discrete Gauss constraint for the loaded
particle charge before the first physical step. Starting from zero E with a
nonzero bunch creates a numerical startup pulse that is not particle
radiation. The default initialization therefore constructs a rigid-beam E/B
field once and then discards every potential/solver work array. The evolution
state remains direct SI E/B; A/phi never enters the time advance.

## Discrete construction

1. Each boosted-frame particle deposits charge on Yee vertices with the same
   first-order CIC shape used by the charge-conserving trajectory depositor.
2. Contributions on duplicated MPI z-interface vertices are summed. Global z
   planes have one owner, so dot products never count a shared plane twice.
3. The represented-mass-weighted mean axial velocity `beta_bar` is evaluated
   in the simulation frame. A matrix-free distributed conjugate-gradient solve
   computes the scalar potential for

   ```text
   -[Dx^2 + Dy^2 + (1-beta_bar^2) Dz^2] phi = rho / epsilon0
   ```

   on interior vertices. This is the axial rigid-beam relativistic-Poisson
   equation used in the Vay formulation. It reduces continuously to ordinary
   Poisson when the boost follows the beam mean velocity.
4. Fields are formed on their native Yee locations as

   ```text
   E = -grad(phi) + beta_bar (beta_bar dot grad(phi))
   B = beta_bar cross E / c
   ```

   E is stored at time zero. B is evaluated at the Yee leapfrog time `-dt/2`
   by the rigid-translation shift `z -> z+beta_bar*c*dt/2`. A second CIC
   deposition and a true cross-rank divergence calculation verify
   `div(E) = rho/epsilon0` before the run can continue.
5. For CPML runs, the inner CPML surface is the static zero-potential surface;
   CPML field cells and auxiliary variables start at zero. This avoids the
   startup flow produced when a static Coulomb tail is loaded into CPML while
   its memory variables remain zero. For PEC regression, the mesh outer surface
   remains the zero-potential surface.
6. The solver verifies that the complete CIC charge remains strictly inside
   the static surface. It also reports bunch RMS resolution and centroid
   padding in RMS units; under-resolved axes and padding below four RMS emit
   explicit warnings.

The scalar potential is only a temporary numerical device for satisfying the
constraint. Four vertex slabs are live during CG (`phi`, residual, direction,
and operator image); the charge lattice is allocated only during deposition
and post-check phases. Startup reports iterations, Gauss residuals, charge,
mean velocity, electric/magnetic energy, RMS resolution, padding and temporary
memory per rank.

## Input

```yaml
initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-10
  maximum_iterations: 10000
```

The block is optional and defaults to these values. The legacy comparison is
selected with `model: electrostatic-poisson`. An enabled solve is a hard
preflight: non-convergence or a failed post-check aborts with a correction
suggestion. Disabling it is allowed only for controlled zero-field/regression
comparisons and prints a prominent warning.

## Physical scope

This step supplies the Gauss-constrained field of a rigid distribution with a
common axial velocity. It is not an exact retarded Lienard-Wiechert solution
for an arbitrary velocity spread, and it does not make an arbitrary sampled
phase-space distribution a Vlasov equilibrium. A finite unconfined electron
bunch has no static free-space equilibrium: physical space-charge expansion
must be separated from initialization error by a matched K=0 drift and by
macro-particle/grid/padding convergence.

Zero potential at the static surface is a finite-domain approximation, not an
infinite-space boundary condition. Quantitative runs must converge the
transverse and longitudinal distance from the bunch to the CPML entrance. A
small Gauss residual proves the discrete constraint; it does not prove that
the box or the sampled distribution is physically adequate.

## Entrance-drift validation

A finite-length K=0 test used `gamma=boost_gamma=1174`, Gaussian peak current
`50 A`, `sigma_x=sigma_y=60 um`, and `4096` macro-particles representing
`10^6` electrons. The reconstructed laboratory bunch reference remained at
`-50 mm`; the first magnet's physical entrance was `34 mm` and its compact
fringe interaction began at `10.0238 mm`. Head-anchored synchronization placed
the common-time laboratory head at `-49.9987 mm`, leaving `60.0225 mm` to the
interaction boundary versus an `11.74 mm` one-cell lab margin. The transformed
bunch occupied only `3.2666 mm` in boosted z and fitted the configured box.

The corresponding laboratory events span `3.8349` light-metres because of
relativity of simultaneity. This is logged as a diagnostic only: it is not a
physical entrance drift and no longer forces the bunch reference to `-2 m`.
For HDF5-v3/v4 Elegant input, the separate physical requirement is that every
record advances forward from `input_plane_z` to the reconstructed snapshot.
When magnets exist, that snapshot must also fit before the first magnetic
interaction region; detector-only and element-free drifts have no artificial
magnetic-entrance requirement.

The updated `-50 mm` CPML run initialized `0.202850 pJ` of field energy and
reached the configured `200 mm` reference stop without a placement failure.
The final field energy was `0.203168 pJ`; the integrated CPML-interface flux
was `-3.92e-6 pJ`. These values retain the earlier conclusion that the residual
sub-percent variation is a convergence target rather than a physical entrance
drift requirement.

Before moving the static boundary to the CPML entrance, the same CPML test
changed field energy by about `17.5%` and showed `-0.035 pJ` flowing back from
the absorbing layer. The large startup effect was therefore a CPML/static-tail
initialization mismatch, not ordinary free propagation. The remaining
`~0.1%` variation is comparable to the coarse-grid ledger residual and is a
convergence target, not a claimed physical energy exchange.

## References

- J.-L. Vay, "Simulation of beams or plasmas crossing at relativistic
  velocity", *Physics of Plasmas* 15, 056701 (2008),
  <https://doi.org/10.1063/1.2837054>.
- WarpX relativistic electrostatic solver documentation,
  <https://warpx.readthedocs.io/en/26.06/theory/models_algorithms/electrostatic_pic.html>.
