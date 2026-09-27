# Standalone SASE bunch generator

## Scope

`sase_bunch_generator` creates an APRL HDF5-v2 common-laboratory-time
particle snapshot for SASE startup studies. It is a preprocessing utility,
not part of the simulation core. The runtime kernel does not generate,
renormalize, or alter the loaded shot noise.

The longitudinal density is a flat top over an integer number of resonant
wavelengths. Transverse position is Gaussian. Optional normalized emittance,
Twiss alpha, and uncorrelated gamma spread provide Gaussian transverse and
energy phase space.

This utility does not model an injector, compressor, wakefield, coherent
microbunching instability, cathode nonuniformity, or quantum photon emission.
Use an Elegant input when those upstream correlations matter.

## Why fully random macroparticles are not SASE shot noise

For `N_e` independent physical electrons, the shot-noise bunching at one
harmonic obeys

```text
sqrt(<|b|^2>) = 1/sqrt(N_e).
```

If `N_macro << N_e` equal-weight macroparticles are placed at completely
random longitudinal phases, their numerical bunching is instead of order
`1/sqrt(N_macro)`. The resulting startup signal is too large and the predicted
distance to exponential gain and saturation is too short.

The generator therefore uses the multibeamlet, multiharmonic quiet-start
method described by Fawley. Each wavelength cell contains `B` independent
beamlets. All members of one beamlet have the same `(x,y,ux,uy,gamma)` and
differ only in longitudinal phase. For requested shot noise through harmonic
`H`, the utility uses `2(H+1)` members:

```text
u_quiet(j,k) = u0(k) + j*pi/(H+1),
delta_u(j,k) = sum[m=1..H] (
  a(m,k)*cos(m*u_quiet) + b(m,k)*sin(m*u_quiet)),
a(m,k), b(m,k) ~ Normal(0, 2/(m^2*N_beamlet)).
```

Using quiet-start order `H+1`, rather than placing the highest loaded noise
harmonic on the boundary of a `2H` phase set, gives the requested harmonics
the Fawley variance and suppresses the adjacent unphysical harmonic that can
couple back into harmonic `H`. Random coefficients are independent between
beamlets and wavelength cells.

The expected local statistic is

```text
<N_cell*|b_m,cell|^2> = 1,  m=1..H.
```

The HDF5 metadata stores the measured value for every requested harmonic.
Only harmonics `1..H` are deliberately given physical shot-noise statistics;
higher harmonics are quiet-start structure, not a physical prediction. Choose
`H` at least as large as the highest harmonic that the planned APRL run will
resolve or diagnose, and verify convergence in both `H` and `B`.

## Transverse and energy phase space

For each beamlet, independent standard normal variables generate a Gaussian
transverse phase space. With rms size `sigma_x`, normalized rms emittance
`epsilon_nx`, and Twiss alpha `alpha_x`,

```text
x  = sigma_x*r1,
ux = epsilon_nx/sigma_x*(-alpha_x*r1 + r2).
```

The same definition is used in y. It gives
`det(Cov[x,ux]) = epsilon_nx^2`. Zero emittance produces zero transverse
proper velocity. Gamma is sampled from a Gaussian with the requested rms
spread, and positive `uz` is reconstructed exactly from

```text
uz = sqrt(gamma^2 - 1 - ux^2 - uy^2).
```

## Build and run

The target is built by default and may be disabled without changing APRL:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DAPRL_BUILD_SASE_BUNCH_GENERATOR=ON
cmake --build build --target sase_bunch_generator -j
```

Example:

```bash
./build/sase_bunch_generator \
  --output particles-sase.h5 \
  --electrons 6.24e8 \
  --length-m 30e-6 \
  --resonant-wavelength-m 100e-9 \
  --gamma 1000 \
  --sigma-gamma 0.2 \
  --sigma-x-m 30e-6 \
  --sigma-y-m 30e-6 \
  --emittance-n-x-m 1.0e-6 \
  --emittance-n-y-m 1.0e-6 \
  --harmonics 3 \
  --beamlets-per-wavelength 64 \
  --seed 12345
```

The effective bunch length is the nearest integer number of resonant
wavelengths and is reported in both the terminal summary and HDF5 metadata.
The macroparticle count is derived, not independently chosen:

```text
N_macro = N_wavelength_cells * B * 2*(H+1).
```

Check this number before producing an x-ray input; short wavelength and a long
bunch can create a large file. Output is streamed in bounded HDF5 chunks.
Existing files are refused unless `--overwrite` is explicit.

The generated file uses unit relative macro weights. The APRL input card must
repeat exactly the same physical electron population:

```yaml
beam:
  input:
    type: hdf5
    file: particles-sase.h5
    electrons: 6.24e8
    position_offset: [0.0, 0.0, -0.05]
```

HDF5-v2 is a common-time relative snapshot, so
`beam.reference.input_plane_z` is not used for this file. The z component of
`position_offset` places the bunch in the APRL laboratory beamline.

## Diagnostics and validation

The `/particles` group records:

- requested and effective bunch length;
- represented electrons, mean current, electrons per cell and per beamlet;
- wavelength cells, beamlets, members and total macroparticles;
- phase-space parameters and random seed;
- observed global `|b_m|`, expected global rms, and mean
  `N_cell*|b_m,cell|^2` for each loaded harmonic;
- `complete=1` only after all records and diagnostics are flushed.

Run the built-in statistical and harmonic-cancellation test with:

```bash
./build/sase_bunch_generator --self-test
```

Production SASE results still require an ensemble of independent seeds and a
scan of beamlets per wavelength. A single generated bunch is one physical
shot, not an ensemble average.

## References

1. C. Penman and B. W. J. McNeil, “Simulation of input electron noise in the
   free-electron laser,” *Optics Communications* 90, 82-84 (1992),
   DOI `10.1016/0030-4018(92)90333-M`.
2. W. M. Fawley, “Algorithm for loading shot noise microbunching in
   multidimensional, free-electron laser simulation codes,” *Phys. Rev. ST
   Accel. Beams* 5, 070701 (2002), DOI
   `10.1103/PhysRevSTAB.5.070701`.
