# Matched source notch for the Cowan critical-mode control

The recommended first control keeps the existing Cowan curls and `c dt = dz`
exactly, and replaces only the longitudinal CIC source/adjoint shape by a
symmetric seven-point convolution:

\[
F(\theta)=1-\sin^6(\theta/2),\qquad
f_{-3\ldots3}=(1,-6,15,44,15,-6,1)/64.
\]

The convolution is real, even, has unit DC gain, and has an exact zero at
`theta = pi`. Its low-frequency error starts at sixth order. With `T` the
longitudinal grid shift, `F = I - [(2I-T-T^-1)/4]^3`.

For the SPARC target photon branch, `theta ~= pi/10`, its amplitude response is
0.9999853447. The product of source and adjoint gather responses differs from
one by -0.002931%. The bunch-phase spatial scale `theta ~= 0.415` has response
0.9999235491, and the response product differs by -0.015290%. These are transfer
functions of the numerical coupling, not predicted final FEL power errors.

For comparison, the five-, seven-, and nine-point families
`1 - sin^(2p)(theta/2)` have source/gather product changes at `pi/10` of
-0.119737%, -0.002931%, and -0.0000717%. The seven-point source uses eight
longitudinal charge/tangential-current nodes and seven normal-current edges;
the nine-point choice uses ten/nine. Two-point Gauss quadrature on each crossed
tensor cell is exact for convolved CIC weights (polynomial degree at most three).

## Required matched operations

Use the same `F` on charge nodes, transverse-current nodes, normal-current
edges, the initial/reference charge distribution, and inlet charge matching.
Use its transpose on the electric gather. Because `F` is symmetric this is the
same stencil. Retain the repaired full source support at admission and the
pending support lookahead; the particle center still follows the original
geometric boundary events.

In a translation-invariant region, `D_z F = F D_z`, and hence

\[
\Delta\rho_F+dt\,\nabla_h\cdot J_F
=F(\Delta\rho+dt\,\nabla_h\cdot J)=0.
\]

If the current is `J_F = F S v`, the electric gather is `S^T F Ebar`, giving
the same mesh work as the particle work. Filtering current alone would break
this pairing; filtering only initial charge or inlet density would break the
source matching. The prepared reference metadata must identify this source
kernel so it cannot be reopened as linear or quintic.

The magnetic gather does not enter the electric work identity. Nonetheless,
the suspected BF mechanism acts through the unaveraged magnetic force, so a
causal control should also filter the grid magnetic field used by the gather,
or independently demonstrate that the raw critical field has been strongly
suppressed. Applying `F` to the native B values before ordinary interpolation
preserves the original time layer and has unit low-frequency response to the
same order. MPI magnetic halos must include the extra longitudinal radius.

## Limits that the analysis must preserve

Changing sources does not change the source-free Cowan dispersion relation.
It does change the finite-grid charge/current/force response, slightly even at
the target scale; the control must report that response and measured cost.

An exact Nyquist zero in the sources does not guarantee that inlet completion,
CPML/reference mismatch, or floating-point errors cannot excite a critical
field independently. The physical control must first check the recorded
Nyquist fingerprint, the same-grid old/new field overlap, and the actual
particle gather. If the critical content was not reduced, unchanged BF cannot
be used to exclude its contribution.

No field damping or modified Maxwell kernel is implied by this proposal.
Any additional damping used later needs its removed energy reported as a
separate numerical channel, not as radiation or particle work.
