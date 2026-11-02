# Cowan lab budget against the verified Yee pairing

Read-only continuation of `yee_lab_sbp_audit`. The actual Maxwell/CPML kernel
is linked from the immutable production archive. Particle/deposition/reference
methods and production defaults are not changed. Only the previous oblique
Cowan vacuum pulse is replayed, at its original matched step c*dt=dz.

Let D(E) = deltaB_Cowan(E) - deltaB_Yee(E). D is evaluated independently from
the existing Cowan smoothing coefficients and stencil. Its prediction is
checked against the actual original-kernel magnetic increment at every stored
time. Diagnostic caches are never passed to the physical field update.

For the same paired energy H and integer-time paired cP used by the Yee audit:

```
delta H + Q_H = C_H
C_H = <B_plus, D(E_old)+D(E_new)>_B / (2 mu)

delta cP + Q_P = C_P - F_Gauss
C_P = [P(E_old,D(E_old)) + P(E_new,D(E_new))]/2
F_Gauss = [F_E(E_old)+F_E(E_new)]/2 + F_B(B_plus)
```

Thus the reference lab residual is independently predicted by
Gamma * (C_H + beta*C_P - beta*F_Gauss). Energy and momentum contributions
must remain separate. The four time/space definition differences from the Yee
audit are measured again to reconstruct the old observer residual.

These Cowan-reference operator terms are numerical channels, not newly added
physical external power. Explaining them does not prove continuum accuracy,
establish a Cowan conserved norm, or explain the full particle gain in 19174.
No matched residual is required to vanish before these terms are measured.
No term is inferred from the old residual.

The diagnostic adds a complete one-plane collar for node-centred components:
an upper node beyond a shared interface is the neighbour's node 1, whereas an
upper half-node is its node 0. Native physics MPI exchange remains unchanged.

Before the wave replay, 64 arbitrary nonuniform 3D fields check finite-step
energy/momentum identities with the intended Cowan terms, and four arbitrary
field grids check the independent stencil against the original kernel. These
are program/math tests, distinct from the actual vacuum propagation result.

Prepare and submit a fresh result directory using `prepare.py`. Allocation is
one exclusive 128-CPU node, 16 MPI x 8 threads, 20-minute hard limit, with binding
checks and the retained original 600/60-second MPI-free teardown watchdog.
All frozen sources, executables, kernel archive and input parameters are hashed.
