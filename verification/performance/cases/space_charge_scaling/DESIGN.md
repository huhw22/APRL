# Charge scaling before extending the shape model

Question: as the physical space-charge signal becomes weaker, does the PIC
error keep a constant fractional size, keep a fixed eV size, or follow a
different dependence? A single dense sphere result cannot establish this.

This first comparison fixes the cold simultaneous rest-frame sphere (R=100 um),
lab gamma, 32768 quiet-start samples and weights, six physical lab planes
(last 1.3261993809708743 m), 160x160x1332 active mesh, 3328-cell reference mesh,
dx/dy/dz/dt, 64/128 guards, CPML, fast reference and both accepted inlet options.
Only represented charge changes: N=1e8, 1e7, 1e6, 1e5. N=1e8 reuses the closed
same-binary uniform_full case, after bitwise checking its generated input.
The new inputs come from each charge-specific reference worldline. The common
cold initial state is fixed; the entry-plane momenta cannot also be identical
after a different amount of physical pre-entry space-charge evolution.

All new cards retain the existing physical kernel and strict relative ledger.
No C++ changes or new build. Each uses one available exclusive node, 16 MPI x8,
28-minute hard allocation and 1400-second PIC limit, with original shutdown
watchdog and self,vader,tcp. No node pinning. H5 is read only after close.

Primary physical observable: the same IDs' lab energy change between the same
fixed lab planes, compared with the theory's own plane-crossing events.
Secondary diagnostic: theory at each actual PIC event's rest time, matching
the earlier reports. Report both; do not compare different physical end times.
Record signed eV bias, its absolute value, bias/theory, rms size, and energy
spread, plus the strict energy ledger. A large ratio near a vanishing reference
signal must not be called a large absolute error.

At fixed propagation distance the cold weak-charge gain is approximately Q^2,
whereas an imperfect reference-field force can have a leading Q term. Thus
error may shrink in eV while its percentage grows; neither a constant percent
nor a constant eV floor can be assumed. Fit local exponents only over measured
points, retaining signs and all source values. No production extrapolation.

Changing length/aspect ratio is a separate question. A nonrelativistic uniform
spheroid semi-analytic Coulomb-explosion model exists (Grech et al., PRE84,
056404, https://arxiv.org/abs/1105.0409). It would need explicit internal-velocity
and retardation scope checks before serving as a nonspherical reference.
Passing the charge scan does not validate an actual long beam or radiation.
