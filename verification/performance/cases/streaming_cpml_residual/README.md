# Residual-input CPML comparison

This verification-only experiment implements the user-approved residual-CPML
comparison. It preserves the existing Cowan field advance, Boris particle
advance, current deposition, precision, world-position CPML-history remapping,
particle queue, and once-computed initial reference. The original full and
cpml modes must reproduce their earlier scientific CSVs byte for byte.

In the new residual mode, only the derivative source values read by CPML have
the frozen reference subtracted at their native field time. Sparse boundary
buffers save those original total-field values and restore them exactly after
the existing CPML calculation. The evolving fields remain total E/B; particles
see those total fields, and currents use the original depositor. There is no
reference-current subtraction, Gauss correction, replacement Maxwell update,
or evolving full-domain field replay.

Geometry and physics are unchanged from the preceding CPML integration:

- 24 x 24 transverse cells, original PEC; longitudinal CPML [0,0,20].
- Full control: 384 nonabsorbing cells in [-64,320] um, plus 20 absorbing
  cells at each end, giving 424 total cells and origin -84 um.
- Candidate: 16 + 2*guard nonabsorbing cells plus 40 CPML cells. Guards
  16/24/40 give total widths 88/104/136.
- 1024 generated particles, seed 710271, charge -0.1 pC, particle box beta
  0 or .23, dx=dy=2 um, dz=1 um, and 176 steps.
- Identical fixed-lab release and detector planes, with window velocity -0.8c
  in box coordinates. Particle beta=0 does not mean a stationary window.

Original post-detector ballistic paths and the deposition collar remain. The
initial reference cache is still a full-domain Cowan solve, translated using
the existing approximation. Changed distributions and omitted exterior return
fields remain limitations to measure against the independent full mode.

## Bounded execution and legacy replay

Enable APRL_BUILD_STREAMING_CPML_PROBE=ON, build aprl_streaming_cpml_probe and
aprl_streaming_cpml_residual_checks, and supply their absolute paths:

```
python run.py prepare --root FRESH_ROOT --binary /absolute/path/to/aprl_streaming_cpml_probe --checks-binary /absolute/path/to/aprl_streaming_cpml_residual_checks --previous-review /absolute/path/to/streaming_cpml_20261027/REVIEW.json
```

The default is three beta0 guard24 preflight cases: full, cpml, residual.
After review, a fresh campaign runs seven explicit --case name:beta:guard:mode
specifications: beta.23 full/cpml guard24, beta.23 residual guards16/24/40, and
beta0 residual guards16/40. This produces ten cases without rerunning zero-supply
negative controls. Only residual mode accepts guards other than24 in this runner.

On authorized gnode1, use APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1 and execute the
frozen FRESH_ROOT/run.py run --root FRESH_ROOT. Shared lock, host-idle gate,
physical-core binding, frozen inputs, process cleanup and 180-second
per-executable watchdog are retained. The campaign limit is two hours. The
optional checks executable is supplied in this experiment and runs after
binding; its successful output must contain Residual CPML checks PASS.

Preparation freezes the old review and the original five scientific CSVs for
every replayed full/cpml case: initial particles, release particles, detector
events, step history, and field-plane data. Any byte difference in a legacy mode
stops the campaign before subsequent cases. Residual outputs also receive
completion-time hashes. All probe/helper/check sources, CMake, core source
hashes, binaries, runtime libraries and local analysis scripts are frozen.

The ordinary cohort, continuity and charge-budget gates remain. All Gauss
residuals are recorded without quantitative acceptance thresholds. Their ordinary
vacuum form is diagnosed only in nonabsorbing interior vertices; the legacy
new_layer_gauss_relative column means the first nonabsorbing interior layer.

## Independent analysis

Use the physics Python interpreter from cluster.env:

```
python review.py --runs PREFLIGHT_ROOT MAIN_ROOT --previous-review /absolute/path/to/streaming_cpml_20261027/REVIEW.json --clamp-review /absolute/path/to/streaming_bunch_20261027/REVIEW.json --output REVIEW.json
python plot.py REVIEW.json --output PLOTS
```

The current independent full-domain run is the primary control. Inputs and
prescribed releases are checked both as numeric bit patterns and canonical CSV
text. Energy differences, Lorentz event transforms, expansion, energy spread,
arrival-current bunching, field differences, queue accounting and CPML clearances
are independently checked or recomputed. The report records absolute B errors
and cB differences normalized to the full-control E-plane peak, so a nearly zero
beta0 reference B cannot alone dominate interpretation.

Earlier same-guard total-field CPML and reference-clamp results are secondary
historical comparisons. New guard24 full/cpml outputs are exact replays of the
preceding implementation. Historical clamp results use fewer total cells because
they omit the CPML layers; their timings are not a matched-cost control.
No approximation threshold is invented. Small particle errors alone do not
accept weak radiation, fields, or net energy exchange.

residual_reference_input_wall_s is already inside field_advance_wall_s and
therefore inside evolution/overall time. residual_reference_work_bytes is
already inside the boundary/grid memory footprint. These are breakdowns, not
extra amounts to add again. Original modes must report zero buffer bytes,
input values and input time. One run per case is a workload observation, not
repeated throughput evidence or production speedup.

Positive particle delta energy means gain; energy loss is its negative. Bunching
uses fixed-plane arrival times, not a simultaneous bunch. Native E(t_n) and
B(t_n-dt/2) remain staggered, both sampled at the detector location at t_n. Plots
do not establish simultaneous lab E/B, Poynting flux, pure radiation, multi-node
behavior, asymptotic complexity, or production readiness.
