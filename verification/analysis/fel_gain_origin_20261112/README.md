# Same-input SPARC gain-origin investigation

The current 11-segment APRL reference uses the unchanged native GENESIS Fawley
loader: 1024 particles/cell, 16 phase clones of each of 64 transverse/energy
states. Only first-harmonic through third-harmonic and quadrant shot-noise tests
were previously accepted. They do not certify the 3D two-particle distribution.

prepare.py freezes a separate campaign using the exact current imported beam,
the same phase-matched lattice and only its first six segments. Four 32-MPI
GENESIS controls share one exclusive 128-core node: SC-off with the previously
validated ADI fix, long-range SC with that fix, long-range SC using the installed
binary, and SC-off with half the longitudinal step. Zero-length diagnostic
markers dump fields at exits 1 and 6 without changing any physical element.
No new APRL simulation or default/method change is made.

run.py checks input and program hashes, each 32x1 binding, literal 6D/current
roundtrip of every active optical cell, finite diagnostics and the tracked
endpoint. Full input mapping and conditional-family correlations are audited.
compare.py uses fixed source-ID groups for PIC BF versus the corresponding
original GENESIS cells. Full-pulse, unwindowed ±5% band energy at the same
±560 um aperture is reported separately for exits 1 and 6. Continuous native
GENESIS pulse energy is clearly labelled as full bandwidth/full aperture.
No independent amplitude matching or power rescaling is used.

The installed GENESIS BeamSolver zeros esc[] and has the shortRange() call
commented out. An input nz=1 does not activate that function. The native
importbeam path reads literal x/y/px/py/gamma/theta/current and does not invoke
the quiet loader; this differs from importdistribution phase reconstruction.
Both facts must be checked against the actual source and imported data, not
inferred from the one4one keyword alone.

sample_input.py performs only a small static read of one central cell and
plots the conditional family harmonics. This sample cannot establish a
whole-beam fraction or the causal share of output power. Historical transverse
reflection controls already demonstrate strong sensitivity of early PIC BF to
this loading family, but they used an earlier method/input and cannot be
claimed as a quantitative intervention on the current F3 campaign.

New physics controls remain evidence pending until their COMPLETE files and
postprocessing are reviewed. A successful import, bitwise diagnostic test or
energy ledger does not establish physical FEL accuracy.
