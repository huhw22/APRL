# Independent w355 six-segment energy review

`review.py` reads only the closed `energy-momentum` run (19796_1). It checks the
boost discrete energy identity, all-rank whole-Cowan energy/momentum identities
and their laboratory Lorentz combination, and an independent fixed laboratory
worldtube using both ports, four sides and two null-time caps. Missing particle
cap energies remain unknown; no residual supplies a missing term.

The field/particle readers and worldtube quadrature are copied from the reviewed
2026-11-13 templates and frozen helper scripts. Changes are limited to reading
the actual completed 28199-step run, handling empty side-event files and recording
cap coverage from this run rather than assuming the earlier p5/p6 outcome.
Only the new independent-w355 data are analyzed. The old phase-clone numerical
results are not reused as physical evidence.

Submit a frozen `scripts` copy through `review.sbatch`, with a separate result
directory and `PROTECTED_OUTPUTS.json` inode/size/mtime guard. `present.py` makes
only single-axis standalone figures; the physical and numerical ledgers retain
different geometries and normalization definitions.

The two-method follow-up can call `review.py PARENT OUTPUT energy-binomial3`
to analyze the closed three-point case. Omitting the final case argument keeps
the original momentum-case default. Each case writes to its own result folder;
the earlier frozen scripts and numerical results remain unchanged.
