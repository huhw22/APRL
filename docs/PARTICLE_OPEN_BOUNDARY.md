# CPML-aware particle boundary

## Experimental transverse boundary injection

The default policy below remains unchanged. For a Lorentz-synchronized bunch
whose transverse extent is much larger at boosted time zero than at the lab
input plane, an explicit opt-in mode is available:

```yaml
mesh:
  ballistic_outside_transverse: true
```

This requires at least two CPML cells in both x and y and is incompatible with
particle retirement. The reconstructed common-time lab snapshot is checked
against the transverse physical aperture with a warning only. The boosted
snapshot may extend arbitrarily far beyond x/y CPML and the outer mesh; z must
still initially fit the physical region. Every particle remains in the MPI
particle list. While outside the physical aperture it advances at its frozen
velocity without grid or prescribed device fields. On entry, only the segment
inside the physical aperture is deposited as Maxwell current, and subsequent
in-aperture steps use the regular Boris pusher. Crossing an MPI z-slab preserves
the same orbit and clipped current path. The final report gives entry crossings
and the number still outside, not unique-particle entry counts.

By default, the initial Poisson solve includes only particles at least one full
cell inside all physical faces. With
`initial_self_field.include_transverse_cpml_particles: true`, particles inside
the transverse CPML strips also contribute to the one-time physical-region
initial E/B; particles beyond the outer mesh remain excluded. The temporary
Poisson grid is padded in x/y, while the evolving CPML fields and memories
still start at zero. If no particles qualify, E/B starts from zero. Omitted
external charge is not represented by the initial field, and boundary-entry
current does not provide an exact
charge-continuity/Gauss-law closure for the whole original bunch. This is a
finite-window approximation, **not** a validated open-boundary PIC injection
condition: perform transverse-aperture, CPML-distance and initial-field scans,
and treat energy-ledger closure and absolute radiation amplitudes cautiously.
When the energy ledger is enabled it reports the raw energy residual, an
independent approximate work source from the in-aperture part of incoming
particle current, and the remaining corrected residual. All particles stay in
the kinetic-energy total. `closure_valid` remains false in this mode; see
[runtime energy accounting](ENERGY_LEDGER.md).
The mode is useful only when omitted external particles are intended to become
physical particles upon re-entry; it cannot recover their prior collective
fields or electromagnetic history.

Energy-ledger format version 3 additionally reports the cumulative absolute
work of the modeled entry current and its ratio to stored physical-region
field energy. This is a source-activity diagnostic, **not** the fraction of
omitted external self-field. A small ratio alone cannot establish that the
omitted field is negligible; aperture and initial-field convergence remain
necessary.

## Purpose

The inner CPML surface is the end of the physical particle domain. A particle
that reaches this surface must no longer contribute trajectory or detector
data, because the layer beyond it is a numerical wave absorber rather than an
observation region. Deleting its charge immediately, however, would create a
static field and an outgoing numerical pulse through a discontinuity in the
discrete continuity equation.

The program therefore separates a physical particle from a numerical current
carrier. The former ends exactly at the first CPML entrance. The latter is a
   small, output-free record which transports and damps the remaining current
inside CPML. CPML-disabled faces retain the charge-conserving open boundary at
the outer computational surface.

## Physical-to-numerical transition

For every pushed trajectory from `x_old` to `x_new`:

1. Find its first intersection with an enabled inner CPML surface or an outer
   computational face.
2. Deposit the physical segment from `x_old` to the intersection with the
   existing first-order charge-conserving trajectory depositor.
3. Store an exact terminal trajectory event at a CPML entrance when trajectory
   output is enabled, and stop particle-plane participation there.
4. Replace the remainder of that same time step with a ballistic numerical
   carrier. It keeps only position, frozen coordinate velocity, charge, and current
   weight; it performs no field gather, Boris push, trajectory output, or
   detector work.
5. Deposit the carrier segment with its damped current weight. Remove the
   residual carrier at the outer face using the same virtual outward-current
   closure as the open-boundary fallback.

There is no one-step pause at the transition. At an exact edge or corner the
first face is selected deterministically, so an event and its charge are
counted only once.

## Matched current damping

Carrier velocity is frozen at the CPML entrance. Along its actual boosted-frame
path the current weight is updated as

```text
w(t + dt) = w(t) exp[-integral sigma(x(t)) / epsilon_0 dt].
```

Here `sigma/epsilon_0` is the same graded electric-conductivity rate used by
the field CPML. Rates from simultaneously active x, y, and z layers are added.
The path integral is evaluated with three-point Gauss--Legendre quadrature;
this is exact for the configured cubic conductivity grading while a segment
remains in one polynomial piece. A transferred segment retains its global time
fraction, so splitting it between z ranks does not restart the damping.

The current deposited over a segment uses the weight evaluated at the segment
midpoint. This is a numerical absorption device, not a physical model of an
electron moving through material. The CPML `kappa` and complex-frequency
shift `alpha` remain field-update parameters and are not folded into the
particle-current matching law.

## Outer cleanup and non-CPML fallback

When the residual carrier reaches the outer box, its final in-domain segment
is deposited and its terminal CIC charge is exported through a virtual normal
current. If `rho_face` is that terminal density, the virtual link supplies

```text
(rho_after - rho_before) / dt + div(J_inside) + div(J_out) = 0,
div(J_out) = rho_face / dt.
```

The virtual link lies outside the E/B lattice and is not an additional Ampere
source. On a face without CPML, a physical particle follows this outer-face
path directly. Extremely small residual carrier weights are discarded at a
roundoff-level threshold.

## MPI, memory, and output behaviour

Internal z-rank interfaces are migration surfaces, not loss surfaces. Both
physical particles and numerical carriers can continue on the neighbouring
rank. Carrier transfer includes the start of the unfinished segment and its
global time fraction so deposition and damping remain continuous.

Each carrier is eight floating-point values: 64 bytes with the current scalar
type. The simulation records only the local carrier vector, its peak size, and
six face totals for each boundary category. It allocates no dense face arrays
and performs no particle-boundary I/O. Face reductions occur once at shutdown.
Trajectory files contain the exact CPML-entry event but never contain carrier
samples.

## Stopping and interpretation

- CPML-entry particles immediately leave the valid-particle set. Numerical
  carriers therefore do not delay `after-last-element`; this mode is intended
  to stop when the physical particle calculation is complete.
- To retain electromagnetic ring-down after all physical particles have left,
  use `reference-center-z` with a target sufficiently downstream. Carriers then
  continue to drive their damped current while the Maxwell solver advances.
- Neither boundary treatment separates radiation from near field. Radiation
  reconstruction should use the retained physical trajectories or laboratory
  detector planes.
- Entry into CPML is a domain-loss diagnostic. Frequent or high-weight entries
  indicate that the physical aperture may be too small.
- The optional relativistic-Poisson initializer supplies a Gauss-consistent
  particle self-field. Boundary continuity cannot repair a disabled,
  unconverged, or otherwise inconsistent initial field.

At normal shutdown the log separately reports physical CPML entries, direct
outer-face exits, residual carrier cleanup, signed charge, and peak carrier
memory for `x-`, `x+`, `y-`, `y+`, `z-`, and `z+`.
