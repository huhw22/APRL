# Runtime particle/field energy ledger

The optional runtime ledger is the global complement to the two-particle-plane
diagnostic in `ENERGY_CLOSURE.md`. It audits the self-consistent simulation on
equal **boosted-frame time** slices instead of trying to infer stored near-field
energy from one downstream laboratory plane.

## Conservation equation

For a run without an incident-wave boundary source or particle retirement, the
ledger evaluates

```text
K_accounted = K_active + K_removed

R = (K_accounted - K_initial)
  + (U_field - U_field_initial)
  + E_field_out
  - W_prescribed.
```

`K` is boosted-frame macro-particle kinetic energy, accumulated with the stable
identity

```text
gamma - 1 = |u|^2 / (gamma + 1).
```

`U_field` is cell-centred E/B energy inside the CPML entrance surfaces; CPML
cells themselves are excluded. `E_field_out` is the time integral of the
outward Poynting flux through those six surfaces. It is accumulated every
Maxwell step with trapezoidal time integration even when complete ledger
records are sampled less often. On a PEC axis, the corresponding flux is set
to zero.

When a physical particle enters CPML or otherwise leaves the physical region,
its interpolated crossing kinetic energy moves from `K_active` to cumulative
`K_removed`. This is a material-energy flux, distinct from electromagnetic
Poynting flux. It prevents particle escape from appearing as unexplained
energy loss.

Laboratory magnetic elements are not stored on the Maxwell grid. Although a
static ideal magnet does no work in the laboratory, its Lorentz transform has
an electric component in the boosted frame. The pusher therefore accumulates
`W_prescribed` along every Boris substep using a trapezoidal line integral of
the prescribed electric field. Omitting this term would make an undulator run
look non-conservative in the computational frame.

`R` is stored as `balance_residual_J`. A zero residual is not built into the
algorithm: charge deposition/current gathering, the Boris push, collocated
energy sampling and the half-step Yee B field do not form an exact algebraic
energy-conserving discretization. The residual must converge when the mesh,
macro-particle count, initial-field padding and device sampling are refined.

## Configuration and cost

```yaml
energy_ledger:
  enabled: true
  directory: output/run-name
  filename: energy-ledger.h5
  sample_interval_steps: 10
  buffer_records: 64
  compression: 0
  warning_relative_tolerance: 0.01
```

When the block is absent or `enabled: false`, no boundary-power loop, particle
work diagnostic, MPI reduction, buffer or output file exists. When enabled:

- the six surface powers and prescribed-device work are accumulated locally
  on every field step;
- the full interior-volume energy, particle moments and MPI reductions occur
  only every `sample_interval_steps`;
- only rank zero opens `filename`;
- `buffer_records` controls its small append buffer;
- interactive shutdown commits and closes a readable file with `complete=0`;
- a normal configured stop closes it with `complete=1`.

The default `warning_relative_tolerance` is 0.01. Shutdown emits a warning when
`abs(relative_balance_initial)` exceeds it, but it does not discard an
expensive run. The warning recommends grid/time-step, macro-particle and
initial-field-padding convergence. `closure_valid=1` means that all source
categories needed by this equation are represented; it does **not** mean that
the numerical residual passed the warning threshold. Incident TF/SF waves
make it zero because their boundary work is not yet in the ledger. It also
becomes zero after experimental particle retirement begins, since that taper
is not charge-continuous.

## HDF5 contract and report tool

The root group is `/energy_ledger`. Readers use only
`records[0:committed_records]`; the scalar `complete` marks normal physical
completion. Each record contains:

- step and boosted time;
- active/removed macro-particle counts and represented active electrons;
- active and cumulative-removed particle kinetic and total energies;
- interior field energy;
- cumulative outward field energy, total and in face order
  `x-,x+,y-,y+,z-,z+`;
- cumulative prescribed-source work and both residual normalizations;
- field fractions relative to kinetic-plus-field energy and to particle total
  energy including rest energy;
- weighted laboratory and boosted gamma moments;
- projected laboratory relative rms energy spread, linear longitudinal chirp,
  and rms spread after subtracting the best weighted linear gamma-z trend.

The small C++ reader prints the initial record, final record and changes:

```bash
./build/energy_ledger_report output/run-name/energy-ledger.h5
```

No Python package is required.

## What energy-spread growth means

Energy spread is a variance, not a separate reservoir in the conservation
equation. The total particle energy is controlled by the weighted mean gamma.
Particles can exchange energy through the self-consistent field so that rms
spread grows while the mean and total remain almost unchanged. Conversely,
coherent radiation can lower the mean while a correlated energy modulation
raises the spread.

The stored projected `sigma_gamma_lab` mixes at least two effects:

1. a coherent longitudinal correlation or chirp, quantified by the weighted
   least-squares slope `linear_chirp_gamma_per_m`;
2. the remaining rms after that linear trend is removed,
   `uncorrelated_sigma_gamma_lab`.

This decomposition is useful for debugging but is evaluated on one equal-box-
time hypersurface. Relativity of simultaneity means it is not a complete
equal-laboratory-time beam diagnostic. A fixed laboratory particle plane is
still required for accelerator-quality entrance/exit phase-space and slice
energy-spread analysis. In particular, do not add the reported laboratory
particle-energy change to boosted-frame field energy: those quantities live
on different spacetime hypersurfaces.

## Gamma-1000 control result

The current 52x52x400, 128-macro, three-period control used a total charge of
one electron and then `10^6` electrons, each with a Gauss-consistent initial
self-field. No particle left the physical region. Values below are boosted-
frame ledger terms at the final stop.

| charge / K | delta K (J) | delta U field (J) | field out (J) | prescribed work (J) | residual (J) | residual / exchange |
|---|---:|---:|---:|---:|---:|---:|
| 1e / 0 | 1.32e-32 | 4.05e-25 | -2.08e-25 | 0 | 1.97e-25 | 0.321 |
| 1e / 0.5 | 5.51e-24 | 4.89e-24 | -7.31e-26 | 9.96e-24 | 3.70e-25 | 0.0181 |
| 1e6 e / 0 | 1.32e-14 | 3.92e-13 | -2.08e-13 | 0 | 1.97e-13 | 0.321 |
| 1e6 e / 0.5 | 1.14e-14 | 4.87e-12 | -7.32e-14 | 4.44e-12 | 3.62e-13 | 0.0385 |

The near-zero `K=0` exchange makes its exchange-normalized percentage
ill-conditioned even though its absolute residual is smaller. The driven
large-charge case closes to about 3.9% of the summed exchange terms, not yet to
a production tolerance. It therefore supports the bookkeeping and sign
convention but also confirms that grid/macro/padding convergence is still
required.

For `10^6` electrons, the initial field fraction including particle rest
energy was `5.53e-5`; final values were `6.01e-5` for `K=0` and `1.15e-4` for
`K=0.5`. In contrast, field divided by **boosted-frame kinetic plus field** was
near one because the chosen boost nearly removes the mean beam motion. This is
not a contradiction: the two denominators answer different questions.

The driven large-charge group ended with projected
`sigma_gamma=4.376e-2`, while the linearly detrended value was
`1.097e-2`; about 93.7% of the projected variance was explained by the linear
longitudinal correlation. The matched `K=0` group produced only
`sigma_gamma=1.43e-4`. Thus most of the observed spread in this control is
coherent longitudinal phase-space structure, not a new store of energy and
not automatically an irreversible slice-energy-spread increase.
