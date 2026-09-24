# Field-plane spectrum and coherence analysis

## Purpose and scope

`field_plane_analysis` is the production-oriented radiation post-processor for
a fixed laboratory field detector. It runs independently of the simulation,
uses threaded FFTW and serial HDF5 on an analysis node, and does not require
particle trajectories. It accepts either:

- the simulator's raw `/field_plane` file;
- a particle-background-subtracted `/reconstructed_field` file; or
- a raw/reconstructed signal plus a geometrically and temporally identical
  zero-radiation baseline, subtracted at electric-field amplitude level.

It writes the positive-frequency forward angular spectrum, integrated energy
spectrum, polarization/Stokes data, a coherent/fluctuation split, selected
angular cross-spectral densities (CSD), and selected two-frequency CSD. This is
the main path for large macroparticle runs; trajectory radiation remains a
small-particle validation path.

## Radiation model

For every retained temporal FFT bin, the tool takes a two-dimensional Fourier
transform of the laboratory electric field over the detector aperture. With

```text
kx = 2 pi qx/(Nx dx),  ky = 2 pi qy/(Ny dy),  k = 2 pi f/c,
nx = kx/k, ny = ky/k, nz = sqrt(1-nx^2-ny^2),
```

only bins with `nx^2+ny^2<1` are forward propagating. The complex electric
field is projected onto a local horizontal/vertical basis transverse to
`n=(nx,ny,nz)`. For the continuous transform

```text
A(kx,ky,f) = integral E(x,y,t) exp[-i(kx x + ky y + 2 pi f t)] dx dy dt,
```

the stored angular spectral energy density is

```text
d^2 W/(dE_gamma dOmega)
  = 2 epsilon0 c k^2 nz^2 |A_transverse|^2 / ((2 pi)^2 h_eV)
```

for ordinary positive-frequency bins. The Nyquist self-conjugate bin uses the
one-sided factor one instead of two. One `nz` is the normal Poynting-flux
projection and the second comes from `dkx dky = k^2 nz dOmega`. The matching
discrete solid-angle weight is `dkx dky/(k^2 nz)`, so summing density times
weight recovers the accepted forward-mode energy crossing the recorded plane.

This is a downstream-vacuum, forward-radiation angular-spectrum model. It is
not a general near-field decomposition and intentionally does not use the Yee
B field, whose detector value is half-step staggered and whose analytic
particle-background part in a reconstructed file is evaluated at the labelled
time. Longitudinal/evanescent content and backward waves are excluded. Place
the detector in vacuum and repeat aperture, transverse resolution, detector
distance and time-step convergence studies.

## Input card

```yaml
input:
  field_file: ../../output/run/detectors/radiation.h5
  # Optional, for the retirement route:
  zero_radiation_baseline: ../../output/baseline/detectors/radiation.h5
  require_complete: true

analysis:
  photon_energy_band_eV: [50.0, 100.0]
  transverse_window: none       # none or hann
  angular_zero_padding: [2, 2]  # [y, x]
  time_windows:
    enabled: true
    interval_s: [5.0e-15, 25.0e-15]
    duration_s: 4.0e-15
    step_s: 2.0e-15

calculation:
  frequency_block: 8
  spatial_batch_points: 32
  fft_threads: 8
  maximum_working_mib: 4096
  maximum_output_gib: 64.0

coherence:
  spatial_photon_energy_eV: [75.0]
  spatial_reference_angles_rad:
    - [0.0, 0.0]
  temporal_photon_energy_eV: [65.0, 70.0, 75.0, 80.0, 85.0]
  temporal_reference_angles_rad:
    - [0.0, 0.0]

output:
  file: field-plane-analysis.h5
  compression: 0
```

Paths are resolved relative to the analysis card. Raw and baseline files must
have identical committed sample times and geometry. The tool linearly
resamples quantized detector times to one uniform grid and reports the largest
correction in units of a time step. A band above the resulting Nyquist photon
energy is a hard error.

`angular_zero_padding` interpolates the FFT k-grid; it does not add aperture or
physical angular resolution. `transverse_window: hann` suppresses a hard
aperture edge but changes the accepted field and is not energy-corrected.
Compare it with `none` rather than silently choosing whichever looks smoother.

`frequency_block` trades memory for repeated field reads and temporal FFTs.
The tool estimates peak working memory and HDF5 size before opening the output;
the two configured limits are hard guards. A larger block is faster when the
analysis node has enough memory. `spatial_batch_points` bounds each HDF5 read.

## Time-window ensemble

With `time_windows.enabled: false`, the complete record is one rectangular
sample. Its CSD is necessarily rank one and cannot demonstrate partial
coherence.

When enabled, every complete Hann window inside `interval_s` is one ensemble
sample. The requested duration and step are quantized to detector samples.
The window amplitude is multiplied by `sqrt(N/sum(w^2))`, so integrated
spectral energy is unbiased for a stationary signal; a narrow band can still
omit Hann side-lobes. Overlapping windows are correlated and do not represent
the same number of independent shots.

Choose the interval only after examining `sample_energy_spectrum_J_per_eV`
and `sample_band_energy_J`. Startup, saturation or pulse decay inside the
interval appears as reduced coherence; that may be a physical question, but
it is not stationary time averaging.

## Output

The root group is `/field_plane_analysis`; `complete` changes to one only after
all products are flushed. Principal axes are `photon_energy_eV`,
`frequency_Hz`, `omega_rad_per_s`, `wavelength_m`, `kx_rad_per_m` and
`ky_rad_per_m`. The angle at each frequency follows the formula stored in the
`angular_coordinates` attribute. A zero `solid_angle_weight_sr` marks an
evanescent/non-forward bin.

Main datasets are:

| dataset | shape | meaning |
| --- | --- | --- |
| `mean_angular_electric_field_spectral` | `[f,ky,kx,2]` complex | ensemble-mean horizontal/vertical aperture spectrum, `V s m` |
| `mean_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | mean angular intensity |
| `coherent_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | intensity of the ensemble-mean complex field |
| `fluctuation_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | mean minus coherent intensity |
| `mean_stokes_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx,4]` | I,Q,U,V with the stored convention |
| `solid_angle_weight_sr` | `[f,ky,kx]` | discrete `dOmega` quadrature |
| `mean_energy_spectrum_J_per_eV` | `[f]` | angularly integrated spectrum |
| `sample_energy_spectrum_J_per_eV` | `[sample,f]` | stationarity/shot diagnostic |
| `coherent_fraction_spectrum` | `[f]` | angularly integrated ensemble-mean-field/mean-intensity ratio |
| `global_degree_of_transverse_coherence` | `[f]` | basis-independent `Tr(W^2)/Tr(W)^2` over angle and polarization |
| `ensemble_gram_matrix_J_per_eV` | `[f,sample,sample]` complex | quadrature-weighted sample Gram matrix with the same nonzero eigenvalues as the full angular CSD |
| `sample_band_energy_J` | `[sample]` | rectangular-bin sum over the retained band |

`/spatial_coherence` stores selected-energy
`angular_cross_spectral_density_J_per_eV_sr[energy,reference,ky,kx,2,2]`
and the basis-invariant `spectral_degree_of_coherence_squared`. Requested
references are mapped to the nearest propagating FFT direction separately at
each frequency; actual angles and indices are stored.

The global degree and Gram matrix avoid materializing an all-to-all angular
CSD. If `a_s` is the solid-angle and polarization weighted angular field of
sample `s`, the stored matrix is

```text
G_st = <a_s,a_t>/N_samples.
```

`G` and the full CSD have identical nonzero eigenvalues. Consequently
`Tr(G^2)/Tr(G)^2` is the global transverse degree of coherence, and
diagonalizing the much smaller sample matrix gives the coherent-mode weights.
This remains distinct from `coherent_fraction_spectrum`, which measures only
the ensemble-mean field and is sensitive to a common carrier-phase jitter.

`/temporal_coherence` stores
`two_frequency_cross_spectral_density_J_per_eV_sr[reference,f1,f2,2,2]`
and its normalized degree of coherence. This is the first-order temporal
coherence information; a mutual-coherence function can be formed by the
corresponding two-frequency inverse transform. Use a dense, uniformly spaced
selection for that operation.

The selected-reference CSD gives spatial cuts through the coherence operator
at bounded cost. The Gram eigenvalues give global coherent-mode weights, but
reconstructing the corresponding angular eigenfunctions would additionally
require the per-window angular fields; those are intentionally not duplicated
in this first storage-bounded output.
