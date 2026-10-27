# Moving-window field-interface diagnosis

The frozen reference is permitted as an injection approximation, but persistent interface-driven field errors must be diagnosed and reduced before any Simulation/MPI integration. Original Cowan, Boris, deposition, precision, physical card and src kernels are unchanged.

`run.py` freezes a bounded single-rank CPU campaign, checks resource/binding/input provenance, and requires old full/cpml scientific CSV replay. Case specifications are `name:beta:guard:mode[:section_steps[:reference_directory]]`. Each case is capped at 180 seconds on authorized gnode1; campaign cap two hours.

Modes: `full` and `residual` reproduce prior probes. `residual-refresh` restores the already prescribed electric faces after remap. `residual-native-b` reconstructs the reference magnetic sample directly from initial E, avoiding repeated interpolation; it is a candidate, not a proven complete fix. `full-record` plus the three `residual-oracle-*` modes isolate exterior/new-layer feed from CPML subtraction using actual full-domain future states: diagnostic only. `held-record` instead evolves an independent prescribed uniform source with original Maxwell and deposition. `residual-held-reference` uses this independent reference while its actual particles still undergo PIC. `residual-held-zero` also prescribes the window particles and is only a zero-response control.

Native record v2 explicitly distinguishes physical-PIC and prescribed-source donors and rejects mismatches. Records retain native E(n), B(n-1/2) and staggered locations. Readers keep only two frames and never extrapolate missing data. Generation still evolves a full domain and has the original length-dependent cost: no asymptotic or production speedup is claimed.

`analyze.py` independently reviews r1 stage sections and legacy replay. `first_step.py` proves the first CPML discrepancy algebraically from recorded fields. `oracle_review.py` reviews r2 channel isolation, field/particle metrics and the local interpolation attempt. `held_review.py` reviews independent-source experiments and clearly separates true PIC from zero-source controls. Scientific errors remain measured quantities, not an automatic acceptance threshold.

All builds and result rounds are independent and frozen. No production switch, MPI integration or remote push is part of this work.
