# Historical ELEGANT 4 m diagnostic analysis

This workflow analyzes existing output. It does not launch APRL, alter the accepted
compute implementation, change an old input, or compare implementation performance.

## Controlling source

The historical folder is `/public/home/wfollowm/APRL_Project/APRL_Elegant_4m_2nodes`.
The complete run is `runs/20261021T070618Z-861230`, recorded simulator revision
`b816b547d3239acb0fd19e893a251f149f60516d` and 256 MPI ranks. It must not be relabeled
as the currently accepted simulator. Two earlier attempts are retained for provenance
and excluded from scientific comparisons: one has no diagnostics, and another contains
truncated/corrupt HDF5 products.

`copy_snapshot.py` copies every file into a new destination, hashes while reading,
re-reads the destination to verify SHA256, checks source sizes/mtimes before and after,
and produces an explicit completion manifest. It never modifies the old directory.
`inspect_history.py` records HDF5 schemas, completion markers and small ledger endpoints;
it does not scan particle or field arrays.

## Questions and definitions

- Five fixed lab crossing surfaces: z=-0.03, 1, 2, 3, 4.03 m. The final surface is
  `/ballistic_reference`, not `/particle_plane`. Each claims 41,732,775 committed
  records. Read committed prefixes only, validate IDs and weights, and compare the
  same population. The beam is a 40 micrometre crop of the original input.
- Momentum is `u=gamma*v/c`; gamma is `sqrt(1+u.u)`. Weights come from recorded
  represented charge/mass and must not be multiplied a second time by macro weight.
  Loss is entry energy minus downstream energy, positive for energy lost.
- Projected RMS energy spread mixes longitudinal correlation and slice variation.
  A linear detrend is only a removal of linear correlation, not proof of irreversible
  heating. Slice statistics use fixed arrival-time bins with reported occupancy.
- Bunching at a fixed plane uses arrival-time phase,
  `b(lambda)=abs(sum(w*exp(-2*pi*i*c*(t-tref)/lambda)))/sum(w)`.
  All particles have the same plane z; a transform of recorded z is meaningless.
  Use a common reference wavelength for cross-plane evolution. The planar-undulator
  reference is `lambda_u*(1+K^2/2)/(2*gamma^2)`, with input K=1.092714657005,
  lambda_u=0.0234 m. The nominal input gamma gives 13.5529046 nm.
- The field detector is at lab z=5.3707685460967935 m and records raw self-consistent
  E/B, including particle self-fields, without the prescribed magnetic-element field.
  Report field amplitudes, recorded aperture, native Poynting flux and finite-window
  integrals. A forward electric-field intensity is an explicitly labeled proxy.
  Lab E/B mix boost-frame leapfrog samples, so shifting lab B by half a time sample
  is not an exact alignment. No near-field-subtracted radiation claim is made.
- The boost ledger uses common boost-time slices, distinct from lab crossing surfaces.
  Check `R=Delta(K_active+K_removed)+Delta U_inside+E_out-W_prescribed-W_boundary`.
  The six outward face integrals are signed net flux, not positive-only radiation.
  Keep changes in stored energy distinct from absolute stored energy.

## Known interpretive limits

This run has `ballistic_outside_transverse=true` and `closure_valid=0` throughout.
The final log reports 26,764,758 particles outside the transverse domain, retained
without grid or prescribed fields. The ledger particle energy still includes these
particles. A small boundary-source-work term does not bound the omitted exterior
self-field. Neither the residual nor the downstream field proves physical convergence.

The field is still substantial at the last stored sample. `complete=1` means the
configured simulation stop was reached, not that a complete radiation pulse was
recorded. Integrals cover only the stored window and aperture. Lab and boost energies
must not be added to form a single conservation equation.

## Reproducible execution

Create a new result directory, make a verified `original/` snapshot, and put frozen
copies of these scripts in `scripts/` with `SHA256SUMS`. Submit `run.sbatch` with the
absolute result directory as its sole argument. Use explicit `--output`/`--error`
locations in that directory. The scripts must not change while a job is queued/running.

The job reserves an exclusive 128-CPU Slurm node, as requested by the workspace
policy, but uses exactly one CPU and one library thread per sequential analysis step.
It validates binding and small synthetic cases before scanning production datasets.
It writes a new `derived/job-*` directory and a completion marker only after all stages
succeed. A completion marker denotes successful execution; scientific interpretation
still requires review of the generated quality checks, summaries and figures.

When no exclusive node is available, preserve the queued job and report its ID and
reason instead of waiting or moving the workload to an unauthorized host.

Physical definitions are grounded in the historical HDF5 schema, recorded input card,
and APRL detector/ledger source. The planar resonance convention is also described in
[DESY's first-operation paper](https://flash.desy.de/sites2009/site_vuvfel/content/e403/e1646/e949/infoboxContent969/first_op_32nm.pdf).
