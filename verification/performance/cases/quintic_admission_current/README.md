# High-order upstream source admission checks

`checks.cpp` directly calls `ConservativeSourceCoupling::depositFrozenSourceSupport`, the same API used by pending and newly admitted high-order sources. It does not change particle trajectories or run a PIC simulation.

The module checks:

- Five outside-center phases plus the reconstructed real source-4097 admission phase. Support prefix plus ordinary inside suffix is compared against an independent full-trajectory current array.
- Charge continuity on retained old nodes k >= 1. A clipped-prefix negative control must reproduce a substantial defect.
- Constraint work against the actual deposited mesh J dot E, with a nonuniform electric field; the support prefix is accounted without center clipping.
- Three discarded current trials and one committed replay. Discarded trials record no work; committed prefix and suffix each record work once.
- A real pending queue with exact repeated additions, including eight pre-admission steps, the admission step, and following inside steps. It must admit once and reproduce the full-source current every step, with no pending/admitted overlap.
- A linear instance rejects the high-order helper without changing current. The original linear inside-source current stays bitwise identical.
- MPI module execution with 1/2/4 ranks, including collective current completion and work accounting; this complements the existing wide-current MPI module.

Compile against the **candidate** `libdirect_eb_core.a`, not the previous binary's library, using the same include/runtime options as `longitudinal_main_flow/build.py`. Run with a fresh output directory for every MPI size. Root will run candidate/old full-program cards separately; these module checks do not establish physical acceptance or prove the production millijoule anomaly has vanished.

Expected JSON: `passed: true`, `physical_acceptance: false`; full-current, retained-continuity and work discrepancies below 2e-10. The negative control must exceed 0.1 in normalized retained-node L1.
