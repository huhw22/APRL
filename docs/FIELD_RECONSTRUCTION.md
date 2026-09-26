# Particle-background field reconstruction

## Purpose

`field_reconstruction` is an offline C++/FFTW tool. It reads:

1. one raw laboratory field-plane HDF5 file;
2. either the field plane's `*-ballistic-reference.h5` file or a laboratory
   particle detector at exactly the same z;
3. one YAML post-processing card.

It writes a new HDF5 file containing a uniformly sampled, particle-background-
subtracted E/B plane and raw/cleaned Poynting diagnostics. Neither input file
is changed. The executable is independent of the simulation and is intended
to run on the lower-cost analysis server.

## Model and scaling

The tool does not sum every particle at every pixel and sample. Each downstream
crossing is deposited once with cloud-in-cell weights onto an `(t,y,x)` grid.
Three real source grids are transformed with FFTW, after which five background
field components are obtained by multiplication with analytical Maxwell
transfer functions. The leading work is therefore

```text
O(number of crossing records) + O(N log N),  N = Nt * Ny * Nx,
```

instead of `O(particles * samples * pixels)`. FFTW threads are configured in
the YAML card. Memory is proportional to the padded field cube, not to the
number of particles; the startup and output metadata expose every padding and
model parameter.

For a particle crossing the observation plane with longitudinal velocity
`v_z`, the constant-velocity Fourier-space scalar denominator is

```text
D = kx^2 + ky^2 + omega^2/(gamma^2 v_z^2).
```

For a narrow-energy, mainly longitudinal bunch, the implementation uses the
charge-weighted mean of `1/(gamma^2 v_z^2)` in `D`. It retains separate
deposited sources proportional to `q`, `q/v_z`, and
`q/(gamma^2 v_z^2)` for the component numerators. The reconstructed field is
the uniform-velocity solution of Maxwell's equations, not an electrostatic
field pasted onto the plane.

This is still a model, not an exact replay of the PIC history. A single
crossing cannot recover earlier acceleration fields, grid dispersion, CPML
history, or the solver's initial-field error. Those are deliberately retained
in the cleaned field as radiation/numerical content. The following guards make
the approximation explicit:

- maximum relative gamma spread;
- maximum RMS transverse beta;
- maximum fraction of charge outside the padded FFT cube;
- user-controlled transverse smoothing and time/x/y zero padding.

The run exits when a guard is exceeded. Padding and smoothing must be varied
until the radiation observables converge.

## Input card

```yaml
input:
  field_file: ../../output/run/detectors/radiation.h5
  particle_file: ../../output/run/detectors/radiation-ballistic-reference.h5
  require_complete: true
  particle_read_chunk: 65536

model:
  fft_threads: 8
  padding_factor: [2, 2, 2]   # [laboratory time, y, x]
  transverse_smoothing_m: 0.0
  maximum_relative_gamma_spread: 0.05
  maximum_rms_transverse_beta: 0.05
  maximum_outside_charge_fraction: 0.001

output:
  file: reconstructed-fields.h5
  compression: 0
  overwrite: false
```

Paths are resolved relative to the card. A ballistic-reference record is
propagated in a straight line from its left entrance to the field plane before
deposition. A `/particle_plane` input must be colocated with the field plane.
For either input, only the first valid downstream crossing of each particle ID
is used. Duplicate detection uses APRL's dense 1-based particle IDs
as an exact bitmap: it costs about 1.25 MB for `10^7` particles rather than the
hundreds of MB typical of a hash table. Duplicates and invalid records are
counted.

Field-plane sample times are generally quantized to Maxwell steps and need not
be exactly uniform. The tool linearly resamples raw E and B onto a uniform grid
spanning the same committed time interval before applying the FFT. The new
time axis and the resampling contract are stored in the output.

## Output

The output group `/reconstructed_field` contains:

- `time_s`, shape `[samples]`, the uniform laboratory time grid;
- `electric_V_per_m`, shape `[samples,ny,nx,3]`;
- `magnetic_T`, shape `[samples,ny,nx,3]`;
- `raw_signed_power_W` and `cleaned_signed_power_W`;
- `raw_forward_power_W` and `cleaned_forward_power_W`;
- `complete`, set only after every field component and diagnostic is written.

Signed power is the transverse integral of `(E cross B)_z/mu0`; forward power
uses `max((E cross B)_z/mu0, 0)`. Both are integrated over stored cell
centres. Attributes `raw_signed_energy_J` and `cleaned_signed_energy_J` give
the trapezoidal time integral of the signed plane flux and are the quantities
used for an entrance/exit laboratory control-volume balance.

Attributes `raw_forward_energy_J`, `cleaned_forward_energy_J`, and
`relative_forward_energy_change` give the trapezoidal time integral and its
relative change. This is the directly evaluated intensity contamination for
the configured aperture and time window; it is not, by itself, an error bar
against an analytical undulator spectrum.

The original detector stores leapfrog B at its half-step time. The current
tool preserves that raw stagger while evaluating the analytical background at
the labelled sample time. This residual is expected to decrease with the main
time step and must be included in convergence tests.

## Application procedure

Use the ballistic-reference file for normal production because it needs only
one compact record per particle. Enable two-plane validation only on a small
bunch to verify that the straight continuation is accurate. A colocated
particle detector can then be substituted as the reconstruction input; the two
cleaned outputs should agree within the orbit-validation error.

For large files, first run with padding `[1,1,1]` to estimate memory and wall
time, then increase the factors. Zero padding reduces circular wraparound but
multiplies FFT memory by the product of the three factors. Compression reduces
disk use but increases analysis CPU time.
