# Three-point source coupling acceptance

This is a separate `linear-binomial3` mode. The historical `linear`,
`linear-notch` (seven point), and `quintic-compensated` modes remain available.
These files prepare program acceptance checks; no result is claimed by this
README before the bounded test jobs complete.

The new longitudinal convolution is

\[
(c_{-1},c_0,c_1)=(1,2,1)/4,\qquad
F_3(\theta)=\cos^2(\theta/2).
\]

It gives four nodal CIC weights, or three normal-current edge weights. Charge,
current, orbit-averaged electric feedback, magnetic feedback, reference charge,
and the inlet use the same shape. The electric feedback remains the adjoint of
the actual mesh current, using the existing mean electric time layer. The
Cowan curl, timestep, particle release events, transverse CPML, and pending
admission schedule are unchanged. Existing eight-layer communication and
source collars are retained in this first implementation.

The exact Nyquist source response is zero. The original homogeneous Cowan
cutoff matrix is unchanged: this acceptance must not be described as a cure
for arbitrary pre-existing field or boundary seeds. Also, an unchanged vacuum
curl dispersion relation does not demonstrate an unchanged coupled FEL
response. For paired source and feedback the extra response is `F_3^2`,
which differs from unity by about -4.8345% at `theta = pi/10` and -8.3082%
at `theta = 0.415`. These are transfer-response changes, not FEL power errors.

## Checks

`checks.cpp` uses an independent convolution of the unchanged linear
charge/current depositor as its numerical oracle. It checks positive compact
basis weights, partition of unity, the charge centroid, the Fourier symbol,
whole-domain deposited charge, multi-cell trajectories, rank interfaces,
retained charge continuity, electric mesh/work adjointness, magnetic Nyquist
rejection, full frozen-source admission support, and reference/inlet charge.
A prepared reference must reject both ordinary and seven-point cache readers.
The unchanged axial vacuum dispersion and the driven cutoff recurrence are
reported separately from particle coupling acceptance.

`admission_checks.cpp` additionally checks six entry phases, the eleven-step
pending-to-active handoff, repeated trial/replay accounting, a deliberately
clipped-source negative control, and an unchanged ordinary linear deposition.
The independent charge oracle here is also a direct three-shift CIC sum.

`build.py` builds the frozen source snapshot, runs both new modules at one,
two, and four MPI ranks, then repeats the historical notch, substep, quintic
admission, and quintic coupling modules at those ranks. It finishes with the
repository's required regression suite. It needs a fresh build directory and
a `plan.json`/`FROZEN.json` prepared by the coordinating task. `build.sbatch`
provides a bounded single-node allocation and the existing Open MPI transport
selection. Cross-node smoke/production startup checks are separate.

A subsequent full-program drift smoke test is needed to exercise YAML,
reference preparation, pending admission, detectors, and the strict ledger
together. Successful module or ledger accounting checks are program acceptance;
they do not by themselves establish FEL physical accuracy.
