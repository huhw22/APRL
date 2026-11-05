# Independent physical-electron source

This tool generates the fixed-lab-plane crossing events used by PIC and exports
the same microscopic sample to Genesis. It does not run a particle pusher or
change field/deposition methods. The first supported analytic distribution is
the declared SPARC parabolic current profile with Gaussian transverse phase
space and Gaussian energy spread.

Each row represents **one electron**, with `macro_weight=1`. A fixed integer N
is sampled from the complete joint distribution. There are no cloned beamlets,
post-generation moment matching, prescribed bunching factors or rescaled shot
noise. For independent weighted samples, noise variance scales as
`sum(w*w)/sum(w)**2`; simply scattering fewer, heavier macroparticles would not
preserve the physical noise level.

## Generation and conversion

```
python generate.py example_37periods.json NEW_MASTER.h5
python genesis.py NEW_MASTER.h5 NEW_GENESIS.par.h5 --wavelength-m 4.913172351123507e-7
python verify.py NEW_MASTER.h5 --genesis NEW_GENESIS.par.h5 --output NEW_VERIFY.json
```

All output names must be new. Complete files are published atomically without
overwriting existing final or partial files. Consumers must check the saved
completion marker; a partially written file is not an accepted source.

`example_37periods.json` keeps the original complete parabolic pulse and selects
its central111 optical slices. It does not compress a fresh parabola into that
window. Physical fixed-N mode requires the integer count to equal the rounded
integrated profile charge divided by e. Diagnostic low-N tests have a separately
named mode; their actual charge is N*e and they are not macroparticle versions
of the high-charge beam.

Generation-time bin populations follow one multinomial draw conditioned on N.
The time and joint transverse/energy coordinates are then drawn independently
inside each bin. Counts are not forced equal per optical period. Separate random
keys are used for counts and particle coordinates. Each source ID consumes eight
Philox uint64 values; three Box–Muller pairs, one time uniform and one reserved
value make the sample independent of chunk boundaries in the frozen runtime.
Cross-platform transcendental-library bitwise identity is not promised.

The stated Twiss parameters and energy spread belong to the parameter plane.
Events are mapped to the output plane by a defined straight, magnet-free drift:

```
x_in = x_parameter - d*ux/uz
y_in = y_parameter - d*uy/uz
t_in = t_parameter - d*gamma/(c*uz)
```

Gamma for the event map is reconstructed from the actual stored on-shell proper
velocity. No space-charge evolution is implied by this preparation map. The
master metadata records both planes, the full config, seed, algorithm/version,
actual N/Q, sampled-bin counts and source-payload checksum.

## Common-source contract

The master is the existing APRL v4 64-byte event record. PIC can read this file
directly; its `beam.input.electrons` must be the actual integer N, so the existing
reader gives each record charge -e and mass m_e exactly.

Genesis output uses native HDF5 slices and `&importbeam`, not
`&importdistribution`, which resamples and randomizes phases in the local
Genesis version. Use `one4one=true`, `nbins=1`, `shotnoise=false`. Varying slice
populations and empty padding are preserved. Slice currents follow the actual
unit-electron population. The converter uses stored arrival times, not the
generation-bin labels, because the preparation drift can move events across
those bins.

The native phase is defined from `s=-c*t` with a saved integer slice origin and
common time anchor; `px,py` are proper transverse momenta p/(mc). A per-slice
`source_id` array is included only for initial auditing. Genesis itself does
not propagate that ID through tracking. The converter records its historical
current-constant difference rather than changing the source to hide it.

The audit ID adds eight bytes per electron to the usual six Genesis doubles.
For60360290 electrons the master plus native conversion therefore use about
6.75GiB payload, before container overhead. Generation is chunk bounded, and
verification additionally uses a one-byte-per-electron ID bitmap. These tools
save repeated generation; they do not reduce the subsequent cost of tracking N
real electrons.

## Checks before large input generation

```
python checks.py NEW_CHECK_DIRECTORY
python genesis.py --write-test-deck NEW_GENESIS.par.h5 NEW_IMPORT_TEST_DIRECTORY
```

`checks.py` verifies independent analytic moments/count/noise statistics over
multiple seeds, chunk-size bitwise equality, absence of copied transverse rods,
deterministic two-format event conversion, variable/empty/boundary slices and
rejection of incomplete/corrupt/weighted inputs. `verify.py` independently
reconstructs events and validates every audit ID. A real native PIC reader and
Genesis import/write test under one and multiple MPI ranks are separate required
integration checks before large generation.

Source generation/import acceptance is not FEL radiation or space-charge
dynamics acceptance. Newly generated input still needs the intended short
physics comparisons before production use.
