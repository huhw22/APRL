Moving-window inlet plateau diagnosis (2026-10-30)

Physical report and successful closed comparisons:
  results/energy_ledger_test_record_20261030/physical_benchmarks/
  charge_inlet_plateau_origin_20261030_r3/README.md

These drivers freeze independent directories before Slurm submission. They
preserve the original Cowan/Boris/deposition/CPML, the input particle events,
and strict energy accounting. Guard, transverse resolution, and transverse
aperture variations are controlled geometry experiments, not default changes.
See the result PROGRAM_VALIDATION.md for all job identities and failed attempts.

Read-only force sampling:
  APRL_PARTICLE_MOMENTUM_PROBE=<new empty directory>
  APRL_PARTICLE_MOMENTUM_PROBE_DT=<actual nominal field dt>
  APRL_PARTICLE_MOMENTUM_PROBE_STRIDE=1
The stride defaults to 20. Changing the nominal dt to fake a denser stride
does not correctly cover particle release substeps and must not be used.
FIRST=1 only selects the first record per source ID per thread; the reader
deduplicates globally by earliest actual event because particles can migrate.

All review*.py scripts read completed outputs. review_dense.py compares the
actual Boris substep with the independent smooth-sphere field applied at the
same actual PIC state. The complementary remainder includes changed paths,
distribution and continuum/discrete differences; it is not one boundary term.
review_position.py gives only a sparse, approximate position-feedback check.
review_transient.py aligns complete-cohort curves at step 123, exposing the
common transient that aperture/mesh changes did not eliminate.

The 1000-step stride=1 output has bitwise-equal science ledger and all 17
particle/trajectory record datasets against the accepted sparse comparison.
Raw dense traces are large and this diagnostic is not a performance setting.
