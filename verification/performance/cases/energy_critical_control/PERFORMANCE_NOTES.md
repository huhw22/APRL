# Cost of the matched seven-point control

This is a code-level cost review, **not a profile or a speedup measurement**.
The closed 257-step smoke sample measured 0.27904 s/step with substep records
and 0.43406 s/step with the notch shape, using steps 40–240. That is about
55.6% more time in this one short paired sample. It is not an estimate of the
full first-segment or production cost.

No implementation is changed by this note. The current control is deliberately
kept frozen while its physical effect is measured.

## Work visible in the actual loops

`longitudinal_shape::weights` divides the straight orbit into tensor-cell
segments. For the linear-notch shape it evaluates two Gauss points per segment.
At each point it visits 16 Jx/Ex entries, 16 Jy/Ey entries, and 28 Jz/Ez entries:
60 callbacks, or 120 per segment. The ordinary first-order analytical orbit
gather/depositor has 12 entries per segment. The electric callback reads both
old and candidate E. These are source-level counts for generic interior
segments; zero displacement, repeated/zero weights, caches, compiler behavior
and boundary clipping change actual executed work.

`ConservativeSourceCoupling::Impl::magnetic` calls the complete trilinear
magnetic sampler seven times. Each sampler performs 24 native reads and
recomputes its cell coordinates and fractions. It therefore has 168 reads per
call, compared with 24 for ordinary CIC. It is called while establishing a
particle force guess and again in each local residual evaluation. The field
does not change during those particle residual evaluations, but their sample
positions can change.

The wide shape also changes the electric and magnetic MPI collars from one to
eight layers and uses an eight-layer current reduction. This is additional
communication and memory traffic, but its share of the 55.6% measurement has
not been profiled. The prescribed-field substep pusher and Cowan loops are
unchanged.

## First option: analytical integration and merged magnetic weights

The filter coefficients are

\[
 c_{-3:3}=\frac1{64}(1,-6,15,44,15,-6,1),\qquad
 F(\theta)=1-\sin^6(\theta/2).
\]

On one orbit cell segment, all CIC coordinates are affine in its local
parameter. Every convolved node weight is also affine:

\[
 w_k(f_z)=c_k(1-f_z)+c_{k-1}f_z.
\]

For any affine factors `a` and `b`, the exact segment average is

\[
 \langle ab\rangle=
 \frac{a_0+a_1}{2}\frac{b_0+b_1}{2}
 +\frac{(a_1-a_0)(b_1-b_0)}{12}.
\]

Thus Jx/Ex use the exact product of a transverse CIC weight and `w_k`; Jy/Ey
use the analogous other transverse coordinate; Jz/Ez use `c_k` times the exact
two-transverse-coordinate product. No Gauss loop is required. The result has
60 visits per segment rather than 120, while retaining the exact same
discrete source shape, orbit integral, and work identity in real arithmetic.
This can reuse the original first-order `integratedProduct` construction.
Floating-point results need a roundoff tolerance rather than a bitwise claim.

For magnetic sampling, compute the original midpoint's transverse fractions
once and convolve its two longitudinal weights into eight. The same result
then needs at most `3*4*8 = 96` native reads, rather than `7*24 = 168`, and avoids
seven repetitions of coordinate normalization and floor operations. Close to
a real open face, merge weights by the **actual clamped donor index** so that
the original constant end continuation is preserved. MPI faces must continue
to use donor data, not clamping.

These are the lowest-risk ways to reduce the present arithmetic. They remove
work inside the existing particle path without adding full-grid sweeps or
changing the mathematical method. The callback-count reductions alone do not
predict a whole-step speedup.

## Second option: move the convolution to the grid

Let `D` deposit ordinary first-order orbit current and `F` denote the seven-point
longitudinal convolution on each native Yee component. In the translation
invariant interior, the current control is

\[
 J_F=FDj,\qquad
 \mathcal G_F(E)=D^T F^T E=D^TFE.
\]

It can therefore be evaluated by ordinary first-order deposition followed by
a grid convolution, and by convolving the candidate time-mean electric field
before an ordinary first-order analytical orbit gather. For consistent native
inner products,

\[
 \langle FDj,E\rangle=\langle Dj,F^TE\rangle.
\]

This is an algebraic rearrangement of the current control, not a change to the
Cowan curl or to `c dt = dz`. A similarly filtered B grid can replace the seven
separate spatial magnetic samples; B remains at its original time layer.
Only the gather cache is filtered: filtered E/B must **not** overwrite the
Maxwell state, detector data, CPML memory, or strict-ledger fields.

The costs then change schematically from

\[
 T_\mathrm{direct}=T_\mathrm{kernel}
  +I\,T_\mathrm{wide\ particle}(P)
\]

to

\[
 T_\mathrm{grid}=T_\mathrm{kernel}
  +I\,T_\mathrm{CIC\ particle}(P)
  +I\,T_F(G_E+G_J)+T_F(G_B),
\]

where `I` is the actually observed fixed-point iteration count, `P` counts
participating orbit segments, and `G` counts component entries filtered. Old E
and B are fixed through a Maxwell-step coupling solve and can be cached once;
candidate E and deposited J change between iterations. Constraints evaluated
after the final solve must use the final filtered candidate.

The small SPARC active grid is `80*80*1742 = 11,148,800` cells while its cropped
input has 236,544 particles, many initially pending. An extra convolution of
every full component can cost more memory traffic than direct wide sampling
for a sparse beam. Therefore full-grid filtering is not automatically a win.
Useful implementation choices are occupied transverse columns or particle
support tiles with a three-cell longitudinal halo, or a measured crossover
between direct and grid paths. Tile validity must cover every candidate orbit
in the current trial, not only the previous iteration's particle positions.
All paths must produce the same source and adjoint gather.

## Conditions that a fast implementation must retain

1. Apply the same `F` to source charge, all current components and the electric
   adjoint; apply the same spatial B gather filter as the frozen control.
   In the interior, `F` commutes with the Yee difference, so
   `delta_t(F rho) + div(F J) = F(delta_t rho + div J)`.
2. A finite physical-domain truncation does not generally commute with `F`.
   Deposit enough **extended** ordinary current before filtering and only
   then truncate. Discarding out-of-domain CIC contributions before filtering
   would lose source tails that should enter the first physical layers.
   The electric gather must be the transpose of exactly that finite mapping.
3. Keep the repaired frozen-source admission prefix, eight-layer pending
   lookahead and common eight-cell safe-exit collar. Internally using CIC does
   not make the physical filtered source narrow.
4. Use the actual staggered node/edge ownership at MPI interfaces; sum source
   contributions once and filter with correct donor collars. Do not double
   count shared Ex/Ey node planes. Global electric and magnetic continuation
   are different and cannot share an arbitrary padding rule.
5. Reference density can be ordinary deposited density followed by the same
   convolution before the existing solve, with correct block halos and
   Dirichlet support checks. It must keep its distinct preparation metadata.
   Inlet charge matching needs the convolution of neighboring charge layers;
   filtering an already completed Ez layer would be a different method.

Before replacing this frozen control, compare direct versus fast charge,
current, E work, B force and boundary/prefix behavior on arbitrary fields and
orbits at 1/2/4 MPI ranks, followed by a same-card coupled run. Check the true
fixed-point iteration count as well as kernel time. The observed BF and actual
critical-mode suppression must agree within roundoff/convergence limits.
No runtime or physical speedup is claimed until those measurements exist.
