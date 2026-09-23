# Charge-conserving open particle boundary

## Purpose

Particles can leave the finite transverse mesh even when the beam core is
well confined. Deleting such a macroparticle at its last in-domain position
would remove charge without the matching current and would introduce an
artificial discrete Gauss-law error. Keeping the particle at the wall would
instead leave a static near-field contribution that is not part of the wanted
outgoing radiation.

The program therefore uses an open, absorbing particle boundary on all six
outer computational faces. This is separate from the electromagnetic
boundary: CPML absorbs outgoing E/B waves, whereas the particle boundary
closes the charge-continuity equation when a particle leaves the mesh.

## Discrete operation

For every pushed trajectory from `x_old` to `x_new`:

1. Find the first intersection with the global computational box.
2. Deposit the segment from `x_old` to that intersection with the existing
   first-order charge-conserving trajectory depositor.
3. Represent the terminal CIC charge by a virtual outward normal current on
   the selected face.
4. Remove the particle from the valid-particle set and record its face, charge,
   and loss count.

If `rho_face` is the terminal CIC charge density, the virtual link supplies
the missing term

```text
(rho_after - rho_before) / dt + div(J_inside) + div(J_out) = 0,
div(J_out) = rho_face / dt.
```

The sign of the particle charge is retained. At an exact edge or corner one
face is selected deterministically, so the charge is exported once rather
than duplicated across two or three faces. The local continuity diagnostic
covers all faces and exact transverse-corner exits.

The virtual normal link lies outside the E/B lattice and is not an additional
Ampere source. The resolved in-domain current up to the boundary remains in
the Maxwell update; only the charge that has actually left the open domain is
exported.

## MPI and memory behaviour

The z decomposition treats internal rank interfaces only as migration
surfaces. A trajectory is split at the interface, its tangential current is
summed by the existing halo exchange, and the particle continues on the
neighbour rank. Only an intersection with a global outer face is counted as a
loss.

Normal production runs retain six local counters and six signed charge sums.
They allocate no face-sized particle-boundary arrays and perform no particle
boundary I/O. Face-resolved MPI reduction occurs only once at shutdown for
the log. Optional dense face-current storage exists solely for explicit local
continuity tests and is disabled by the simulation driver.

## Interpretation and limits

- Absorption occurs at the outer computational face, not at the CPML entrance.
  CPML and the particle boundary therefore do not silently redefine the
  physical aperture.
- This mechanism prevents an artificial charge discontinuity; it does not
  delete radiation or near fields that were already present on the E/B grid.
  Those fields must propagate through CPML in the ordinary Maxwell update.
- It is not a radiation/near-field separator. Radiation reconstruction should
  still use the retained trajectories or suitably placed laboratory detector
  planes.
- The current development solver still lacks a Gauss-consistent initial
  particle self-field. Boundary continuity cannot repair an inconsistent
  initial field, so that item remains a separate prerequisite for final
  self-consistent production results.
- Frequent particle motion inside CPML is a warning that the transverse box
  or physical aperture is too small for the intended observation region. A
  later dedicated particle-buffer policy may absorb particles before the
  field PML, but that would require a separately defined physical surface and
  is not assumed here.

At normal shutdown the log reports loss count and net escaped charge for
`x-`, `x+`, `y-`, `y+`, `z-`, and `z+`.
