# Trajectory-to-far-field radiation tool

## Scope

`postprocess/trajectory_radiation` is a standalone C++/MPI program. It does not
link the simulation core and does not read the Maxwell grid. Its input is one
or more sets of laboratory-frame trajectory HDF5 files; its output is a single
HDF5 file containing the complex polarized far-field spectrum, spectral and
angular energy density, Stokes parameters, angle-integrated energy spectrum,
and optional first-order coherence matrices.

The separation is intentional: the expensive simulation stores reusable
particle histories once, while frequency ranges, angular grids, observation
distance, and ensemble statistics can be changed on a smaller analysis server
without rerunning particle dynamics.

This tool computes radiation fields only. It does not include the Coulomb or
velocity near field, so a particle moving uniformly produces zero output.

## Radiation model

For observation direction `n`, define

```text
F(beta) = n cross (n cross beta) / (1 - n dot beta).
```

The sampled path is treated as piecewise linear in laboratory spacetime. At
each internal trajectory record, the dimensionless chord velocity changes
from `beta_before` to `beta_after`. The charge-weighted spectral radiation
amplitude is accumulated coherently as

```text
A(omega,n) = sum_particles q sum_internal_knots
             [F(beta_after) - F(beta_before)]
             exp{i omega [t - n dot r/c]}.
```

This is the endpoint representation of the far-zone Lienard--Wiechert
acceleration field for the piecewise-linear path. The first and last artificial
start/stop endpoints are suppressed: the recorded path is assumed to continue
inertially outside its time range. Consequently, at least three records are
required per particle. CPML terminal events are valid final points, but the
non-physical CPML carriers never enter the calculation.

The program uses `long double` for retarded time, phase reduction, local
coherent accumulation, and MPI reduction. Output is float64. The trajectory
sampling must still resolve changes in velocity and the prescribed magnetic
field. Frequency-grid convergence and trajectory-cadence convergence are both
required for production results.

The stored complex electric spectrum at distance `R` is

```text
E_tilde(omega,n;R) = exp(i omega R/c) A / (4 pi epsilon_0 c R),
E_tilde = integral E(t) exp(+i omega t) dt.
```

The corresponding positive-frequency spectral-angular energy is

```text
d^2 W / (d omega d Omega)
  = |A|^2 / (16 pi^3 epsilon_0 c).
```

The implementation follows the arbitrary-trajectory spectral method described
by A. G. R. Thomas, *Phys. Rev. ST Accel. Beams* **13**, 020702 (2010),
<https://doi.org/10.1103/PhysRevSTAB.13.020702>. The present first kernel uses
the piecewise-linear endpoint form; a later higher-order path interpolant can
be added without changing the output contract.

## Build and run

Build it independently of the simulator:

```bash
cmake -S postprocess/trajectory_radiation \
      -B build-radiation \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-radiation -j
```

Run locally or with MPI:

```bash
./build-radiation/trajectory_radiation radiation.yaml
mpirun -n 16 ./build-radiation/trajectory_radiation radiation.yaml
```

Every input file is opened only by one analysis rank. Records are redistributed
by `particle_id`, so a trajectory that migrated between simulation ranks is
reassembled before radiation is calculated. Output is written only by analysis
rank zero. Frequency and theta-y blocking bound temporary field memory.

Trajectory records owned by an analysis rank remain in memory for the current
shot; the blocking applies to the radiation array, not to the trajectory set.
The one-rank path has no MPI message-size restriction. With multiple ranks, a
single peer-to-peer redistribution bucket must fit the MPI signed-`int` byte
count (about 2 GiB); increase the analysis rank count before that limit is
reached.

## YAML input

The complete example is
[example.yaml](../postprocess/trajectory_radiation/example.yaml).

```yaml
input:
  require_complete: true
  read_chunk_records: 65536
  shots:
    - name: shot-000
      files:
        - particles-rank-00000.h5
        - particles-rank-00001.h5

observation:
  axis: [0.0, 0.0, 1.0]
  horizontal: [1.0, 0.0, 0.0]
  distance_m: 100.0
  theta_x_rad: {min: -0.001, max: 0.001, count: 101}
  theta_y_rad: {min: -0.001, max: 0.001, count: 101}

spectrum:
  photon_energy_eV:
    min: 100.0
    max: 1000.0
    count: 181
    spacing: linear

calculation:
  frequency_block: 16
  theta_y_block: 4
  minimum_records_per_particle: 3

output:
  file: radiation-output.h5
  compression: 0
```

Paths are relative to the radiation YAML file. A shot is one statistically
independent realization and can contain any number of per-rank trajectory
files. Trajectory formats 1 and 2 are accepted. Setting `require_complete` to
false permits an intentionally interrupted but committed trajectory prefix.
All shots used in a field average or two-frequency CSD must share a meaningful
laboratory time origin. A deliberate random arrival-time distribution is
physical; accidental per-file time resets will instead erase phase-sensitive
ensemble quantities.

The central `axis` and projected `horizontal` vector define an orthonormal
basis. Grid directions use

```text
n = normalize(axis + tan(theta_x) horizontal + tan(theta_y) vertical).
```

Photon energy can use `linear` or `log` spacing. Angular axes are linear.

## HDF5 output

The `/far_field` group has `format_version=1` and a scalar `complete` marker.
Major axes are `photon_energy_eV`, `omega_rad_per_s`, `frequency_Hz`,
`wavelength_m`, `theta_x_rad`, and `theta_y_rad`. `observation_basis` records
the central, horizontal, and vertical basis vectors.

`complete` remains zero until all per-shot fields, ensemble summaries, spectra,
and coherence products have been written and flushed. An interrupted file may
therefore contain useful diagnostic blocks, but must not be treated as a
complete radiation result or as a restart file.

The principal datasets are:

| dataset | shape | meaning |
| --- | --- | --- |
| `electric_field_spectral` | `[shot,f,y,x,2]` complex | per-shot horizontal/vertical `E_tilde`, V s/m |
| `spectral_energy_density` | `[shot,f,y,x]` | per-shot `d²W/(dω dΩ)`, J s/sr |
| `stokes_spectral_energy_density` | `[shot,f,y,x,4]` | per-shot I,Q,U,V in the same energy unit |
| `mean_electric_field_spectral` | `[f,y,x,2]` complex | ensemble mean field |
| `mean_spectral_energy_density` | `[f,y,x]` | ensemble mean intensity |
| `coherent_spectral_energy_density` | `[f,y,x]` | intensity of the ensemble-mean field |
| `fluctuation_spectral_energy_density` | `[f,y,x]` | mean intensity minus coherent intensity |
| `energy_spectrum` | `[shot,f]` | angularly integrated `dW/dω`, J s |
| `mean_energy_spectrum` | `[f]` | ensemble mean `dW/dω`, J s |
| `mean_energy_per_log_frequency` | `[f]` | `dW/dln(ω)`, J |

The complex datatype has `real` and `imag` float64 members. Polarization index
0 is the projected horizontal basis and index 1 completes the right-handed
transverse basis. The stated Stokes convention is stored as an attribute.
Angular integration uses the exact Jacobian of the two tangent-angle
coordinates and trapezoidal weights. It is unavailable when either angular
axis contains only one point.

## Temporal and spatial coherence

First-order coherence is an ensemble property:

```text
W_ab(p,q) = mean_shots[conj(E_a(p)) E_b(q)].
```

For a single deterministic shot this outer product is necessarily rank one;
it cannot by itself establish partial coherence. The full per-shot complex
field is therefore always retained. With several statistically independent
shots, selected cross-spectral-density matrices can also be generated:

```yaml
coherence:
  spatial_photon_energy_eV: [500.0]
  spatial_reference_angles_rad:
    - [0.0, 0.0]
  temporal_photon_energy_eV: [450.0, 475.0, 500.0, 525.0, 550.0]
  temporal_reference_angles_rad:
    - [0.0, 0.0]
```

`/coherence/spatial_cross_spectral_density` stores the 2x2 polarization CSD
between each reference angle and every angular-grid point at selected photon
energies. `/coherence/temporal_cross_spectral_density` stores the 2x2 CSD
between every selected frequency pair at each reference angle. Requested
coordinates are mapped to the nearest actual grid point and the selected
coordinates are written alongside the matrices.

Cross-spectral density and mutual intensity require an ensemble average; this
definition and its connection to spectral density are reviewed in
“Coherence properties of the high-energy fourth-generation X-ray synchrotron
sources,” *J. Synchrotron Radiat.* **26** (2019),
<https://doi.org/10.1107/S1600577519013079>.

## Macroparticle and convergence cautions

- `charge_C` already contains the charge represented by a macroparticle and is
  summed coherently. The diagnostic `weight` is not multiplied again.
- Ensemble fluctuation terms describe physical shot noise only if the input
  realizations and macroparticle model represent that noise. Resampling a
  coarse smooth distribution mainly measures numerical macroparticle noise.
- The first/last endpoint suppression assumes inertial continuation. Record
  enough drift before and after the radiating element so the retained internal
  velocity changes cover the full acceleration region.
- A zero far-field result from a two-record particle is deliberate: one chord
  contains no resolved internal acceleration.
- Compare successively finer trajectory cadence, angular spacing, angular
  aperture, and photon-energy spacing. Coherent FEL phase is particularly
  sensitive to trajectory cadence and retarded-time accuracy.
- This is a far-zone tool. Near fields, finite-aperture Fresnel propagation,
  optical elements, and detector response require separate propagation tools.
