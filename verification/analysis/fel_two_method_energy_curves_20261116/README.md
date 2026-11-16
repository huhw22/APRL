# Two-method six-segment energy curves

Read only the completed independent-w355 cases 19796_0 and 19796_1.
`three_point.sbatch` runs the generalized validated review reader for the first
case. `render.py` reuses the earlier momentum-case caches and creates separate
standalone bilingual figures, with both methods on each comparison chart.
`finalize.py` checks signs, decomposition, unchanged sources and output identities,
preserves the earlier momentum values, and packages only figures and records.

The time curves retain all 28200 records. The fixed laboratory curves connect
the seven measured planes only. The static external impulse cancellation is
algebraic; neither a measured substep force trace nor a same-domain decomposition
of the two laboratory residuals is invented. The quintic method's historical
momentum-matching label is preserved with its earlier qualification.
