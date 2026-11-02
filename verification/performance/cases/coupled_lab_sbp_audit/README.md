# Coupled conservative-source Cowan lab diagnostic

This case retains the original short r5 production-flow cards and original
physics library. An observer-only overlay is materialized on the frozen r5
source. Only simulation.cpp orchestration with read-only hooks is rebuilt;
field, particle, deposition, coupling and boundary kernels stay in the existing
archive. Original raw field and three particle planes must remain bitwise equal.

The reference field functional is H=eps/2<E,E>+1/(2mu)<B,B>+1/(2mu)<B,K_C(E)>
and P_ref=P(E,B)+P(E,K_C(E))/2. K_C is the original diagnostic Cowan stencil
with true MPI collars and its default exterior convention for this diagnostic
functional. Real physical magnetic reference halos are unchanged; their exact
increment difference is recorded as an explicit magnetic supply source.

Sources are measured independently from actual committed operations: magnetic
reference supply, CPML, synchronization, global tangential end-node holding,
electric CPML/inlet/synchronization, and transverse reference Jz. Reference Jz
is captured before/after the committed outer source call, separate from all
trial iterations and actual particle current. Remap is an independently
evaluated jump of the two reference functionals.

Energy source pairing uses the full-step mean E or B; momentum source pairing
uses P(deltaE,Bplus) and P(Eold,deltaB). Source values are not derived from the
final balance residual. A 64-case arbitrary 3D finite-step test checks these
source identities before the two actual PIC runs.

The clean all-active interval excludes frozen/carrier/handoff/removal and is
used for the first particle-field attribution. Original diagnostics remain,
and the immutable momentum replay independently measures committed particle
work/impulses. Earlier phases are not silently accepted as particle-conserving.
Matching a discrete budget does not establish continuum accuracy or resolve
the full 19174 particle gain.

Workflow: prepare.py creates fresh output/build directories from the existing
r5 private diagnostic source and archive; run.sbatch/run.py runs only in the
bounded exclusive allocation. post.py accepts completed scalar/rank records
and closed HDF5, checks original science identity and reports both operator and
coupling differences. review.py plots all clean-interval curves with one fixed
full-record particle-change denominator per case, with the same denominator
in the numerical-detail panels. Absolute J records remain available.

The 2026-11-02 campaign is in results/lab_energy_balance_20261102/coupled_sbp_r2.
Job 19227 completed both scientific segments but its frozen reader applied
Gamma twice to an already lab-projected identity check. The corrected post.py
was run read-only by 19228, without changing the frozen scientific inputs or
loosening the threshold. The old batch failure/reader remain in that result
directory; analysis_r1 holds the corrected reader provenance.

These are dedicated diagnostic builds; no production default is replaced and
no production performance improvement is claimed. C++ modifications are
observer hooks in simulation orchestration and the private observer header.
All numerical archive members except simulation.cpp.o remain byte-identical
to the r5 diagnostic baseline. Full-record particle conservation, continuum
lab accuracy, and the 19174 root cause remain unaccepted.
