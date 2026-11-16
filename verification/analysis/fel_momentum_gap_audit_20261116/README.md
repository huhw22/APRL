# Existing-data momentum-gap audit

This workflow does not run particles or Maxwell evolution. It checks the
completed independent-w355 case 19796_1 against its frozen source, separates
energy work from longitudinal impulse pairing after the static external-magnet
cancellation, and independently integrates the recorded worldtube currents.

`audit.py` uses inverse Lorentz rays and exact Simpson quadrature for ports and
time caps, then reuses the frozen observer integrator for particle-free gates.
`side_check.cpp` independently integrates the four sides using affine bounds
and two-point Gaussian integration without polygon moments. `summarize.py`
checks the frozen source and writes the report and standalone bilingual plots.

The actual external impulse is inferred from the static-magnet identity,
not measured from an unavailable substep trace. Local electric/magnetic force
attribution and a same-domain native/continuum worldtube bridge remain unknown.
The two different laboratory residuals are never subtracted to invent such a
decomposition. Previous reports, raw data and frozen runtime code stay intact.
