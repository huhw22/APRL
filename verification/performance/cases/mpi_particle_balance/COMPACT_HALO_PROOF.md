# Three communicated layers: nominal support plus an IEEE access guard

Production `linear-binomial3` balancing communicates **three layers**, compared
with the original eight. In real-number arithmetic the nonzero source/gather
support needs only two. That result is **not a guarantee about the existing
floating-point accesses**. The third layer preserves the original arithmetic,
including tiny nonzero endpoint weights and the ordinary magnetic sampler's
zero-weight corner reads. Production does not skip corners or threshold small
weights. Other shapes retain their original eight-layer path.

This changes communicated storage only: convolution coefficients, native Yee
locations, orbit quadrature, conservative velocity, and Cowan curl stay the
same. The compact path requires the original `c dt = dz` and one-cell
longitudinal displacement contract.

## Real-number two-layer support bound

On a conservative substep,

\[
\Delta\mathbf x_s=c\,h_s\,
 \frac{\mathbf u_s+\mathbf u_{s+1}}{\gamma_s+\gamma_{s+1}},
\qquad \gamma_s=\sqrt{1+|\mathbf u_s|^2}.
\]

The triangle inequality gives

\[
|\Delta z|\le\sum_s|\Delta z_s|
\le c\sum_s h_s\le c\,dt=dz.
\]

The initial fixed-point guess also obeys this bound. The grid does not shift
inside the gather: the moving-window shift follows the completed step, and
ownership is assigned on the shifted grid before the next step.

In dimensionless local coordinates, a slab has `n` cells. For a source starting
in `[0,n]` and moving at most one cell, its orbit and midpoint satisfy

\[
q\in[-1,n+1],\qquad m\in[-1/2,n+1/2].
\]

Internal ownership is normally half-open `[0,n)`. This bound does not authorize
starting a source arbitrarily far outside its owner. Pending upstream sources
and real outer faces retain their existing truncation/continuation rules.

For `j=floor(q)`, the three-point convolution of nodal CIC occupies `j-1..j+2`;
the normal-current edge basis occupies `j-1..j+1`. Nonzero support is therefore

\[
\begin{aligned}
(\rho,J_x,J_y,E_x,E_y)\text{ nodes}&:\quad [-2,n+2],\\
(J_z,E_z)\text{ edges}&:\quad [-2,n+1].
\end{aligned}
\]

Node arrays contain `0..n` and edge arrays `0..n-1`, hence nominally two ghost
layers at each end. Current and its electric adjoint use these same weights.

The magnetic gather evaluates `m+s`, `s=-1,0,1`. Its z interpolation coordinate
is `m+s` for `Bz`, and `m+s-1/2` for half-grid `Bx,By`. Their nonzero supports
are respectively `[-2,n+2]` and `[-2,n+1]`, also nominally two layers.

At exact `m=n+1/2,s=1`, however, ordinary half-grid interpolation reads the
upper corner `n+2` even though its weight is exactly zero. That is already a
third edge ghost access. A nominal two-layer special sampler can skip exact
zero corners, but this alone does not solve the floating-point case below.

## Observed IEEE endpoint: a real nonzero third layer

The frozen r4 four-rank test fails on rank 2's path 48, using production
`dz=1.2543262415359034e-6 m`, local `n=2`, and the intended endpoints `n` to
`n+1`. Its actual doubles are

| Quantity | Value |
|---|---:|
| Local world origin | `5.017304966143613e-6 m` |
| Begin z | `7.5259574492154205e-6 m` |
| End z | `8.780283690751325e-6 m` |
| Normalized begin | `2.0000000000000004` |
| Normalized displacement | `1.0000000000000004` |
| Normalized end | `3.000000000000001` |

The existing cell-splitting quadrature consequently contains a final interval
of width `8.881784197001252e-16`. Its second Gauss point is
`q=3.000000000000001`, with `floor(q)=3`. Node 5 has basis weight
`2.220446049250313e-16`; including the quadrature factor its weight is
`9.860761315262648e-32`. It is tiny but **not zero**. Normal-current edge 4 also
gets a nonzero third-layer contribution. Removing either would change the
existing deposition and its electric adjoint.

In the original four-rank/two-cell layout, some outlying contributions would
later be discarded at the true global face, but the two-layer storage check
fails before that routing step. The independent regression preserves exactly
these world coordinates while widening downstream slabs, so both third-layer
currents remain inside the global physical domain. It requires nominal depth
two to reject, depth three and eight to agree bitwise, and both retained
third-layer currents to be nonzero. With two MPI ranks the same coordinates
are reproduced using a common world origin of `4 dz`; with four/eight ranks
the original rank-2 coordinates are used directly.

## Production rule

Three layers retain the two-layer mathematical support plus one layer for
these ordinary IEEE endpoint effects and original magnetic corner accesses.
The magnetic sampler keeps its original corner loop, order, and coefficients;
production does not use the nominal-two exact-zero skip option. The source and
electric gather keep all existing nonzero weights and arithmetic.

This is a checked access contract, not a claim that one guard layer covers an
arbitrary loss of coordinate precision. A genuinely required fourth-layer
access throws. A multi-cell orbit also throws rather than silently altering
its support. Coordinate-scale tolerance only accommodates ordinary roundoff
in the one-cell displacement check; actual support remains independently
bounds checked.

At true global faces, the original electric/source truncation and magnetic
constant endpoint continuation remain. MPI interfaces use donor values.

## Independent acceptance in `compact_checks.cpp`

The production comparison is depth **three versus eight** on mixed 2/3/8-cell
slabs, with the production spacing and all original endpoint cases retained:

- 49 full one-cell source chords per rank, both directions, exact interface
  endpoints and ULP-near endpoints, including true global ends.
- The explicit IEEE regression above: depth two rejects, three/eight preserve
  the same retained nonzero third-layer currents bitwise.
- Production-grid magnetic midpoint/three-shift samples preserve the ordinary
  sampler, including zero corners. A genuinely nonzero fourth layer rejects.
- Constraint electric work versus actual mesh current, with distinct old/new
  electric fields and shared native planes.
- Electric-only, magnetic-only, and mixed initial/residual particle forces,
  particle states, and replay work; pure magnetic energy is checked separately.
  An intended endpoint can normalize just outside a real global sacrificial
  face. The test retains that seed as a negative control: when an independent
  interval check shows its entire initial orbit misses the domain, both three
  and eight layers must reject with the same original error. A rejection by
  only one mode always fails. Valid in-domain endpoints still require the
  original bitwise state, force, and replay comparisons.
- Multi-cell current and gather negative controls stay within a large slab,
  so they test the displacement contract rather than only an outer face.

A separate **nominal-two** magnetic test uses binary-exact `dz=2^-20`,
`dx=dy=4dz`, integer-cell origins, and an exact-zero skip sampler. It must
exercise the zero third-corner case. This isolates the real-number limit and
is explicitly not the production communication method. It still rejects
actual nonzero third-layer support; no tolerance treats small weights as zero.

Original planner/mixed-eight-layer and old-mode modules remain intact.
Successful checks establish a communication substitution and its access
bounds, not long-FEL physical accuracy or a measured speed improvement.
