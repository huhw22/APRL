Full-size Cowan interior-path comparison

Run one exclusive two-node Slurm allocation with 256 MPI ranks, one thread
per rank. Compare the accepted current production executable with the
independently validated Cowan interior-path executable on the same nodes.
This is not a replay of the historical pre-optimization version.

Preparation freezes the compact production card, original input identity,
both binaries and scripts. The full-size cards change only mesh.duration
and output paths: each advances 600 complete steps, retaining the original
physical stop and observation-bounded mesh, held injection, full 41,732,775
particles, resolution, 0.3 m retirement layer, and diagnostic cadence.

Run binding checks and a small paired held/retirement correctness check
before either full-size run. No profiling is enabled. One allocation has a
three-hour hard wall limit; each full-size run additionally has an 80-minute
limit including initialization and shutdown. Limits are maximum durations,
not runtime predictions. Any failing gate stops the remaining sequence.

The short duration deliberately ends before the physical production stop.
Its nonzero simulator exit is accepted only when the exact exhaustion
message, complete 600-step records, preserved geometry, finite committed
readable output and incomplete flags are all verified. A wall timeout is
not a pass. No production simulation is launched after this test.

Compare elapsed-time differences over matching steps 50-100, 100-400 and
400-600, and compare same-rank scientific output exactly. This is one pair,
not a repeated estimate, and covers early injection/undulator evolution,
not the full retirement or far-field observation stage. Existing production
job 18893 and default program paths are not modified.

Read plan.json for the frozen inputs and limits, slurm-JOBID.out for phase
progress, control/aprl.log and trial/aprl.log for per-50-step elapsed times,
and comparison.json for the automatic paired result after completion.

Fixed-plane timestamp checks use the stored timestamps' binary64 rounding
resolution and the accumulated box-time arithmetic error. They check the
entire sequence against the first time plus n*dt/gamma. A relative tolerance
on adjacent differences alone was invalid for 19 ns absolute timestamps
with 4.5 as increments. This changes only the test checker; scientific
control/trial arrays and attributes are still compared bit for bit.

If a completed control was rejected solely by that checker, resume_trial.py
can prepare a new directory containing an independent copy of its closed
outputs. resume-trial.sbatch runs only the missing full-size trial, with
fresh small preflights and binding on any two available 128-CPU nodes,
as requested by the user; neither submission nor runtime pins node names.
This continuation has a 90-minute total cap and an 80-minute full-run cap.
The resulting comparison records both node sets and separate allocations.
Timing ratios may include differences between nodes and cannot be presented
as an isolated same-node implementation-speedup estimate. The original
failed job and checker logs remain unchanged.
