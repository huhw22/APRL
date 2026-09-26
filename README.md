# APRL — Accelerator Particles and Radiation in Lorentz Frames

[中文说明](README.zh-CN.md)

APRL is a C++/MPI program for accelerator-particle and radiation simulations
in Lorentz-boosted frames. The primary executable is `aprl`. Electromagnetic
evolution uses SI electric and magnetic fields directly; scalar and vector
potentials are not part of the time-advance state.

## Scope

Version 1 provides the following computation paths:

- self-consistent relativistic particle and Maxwell evolution;
- Cowan/CKC controlled-dispersion propagation along the configured z axis;
- compact CFS-CPML field boundaries and charge-continuous particle escape;
- analytical laboratory magnetic elements transformed into the simulation
  frame;
- HDF5 particle input, deterministic Gaussian test input, and native Elegant
  SDDS5 conversion through the official SDDS library;
- fixed laboratory field and particle detector planes;
- field-plane radiation, spectrum, polarization, coherence, and energy
  analysis on a separate analysis host;
- optional small-particle trajectory radiation and global energy accounting.

The general Cowan/CPML TF/SF laser or seed-field injection path is not
implemented in version 1. Seeded-FEL and laser-modulation calculations that
require that boundary source are outside the validated scope.

## Numerical model

- E and B are advanced on a staggered Yee lattice in SI units.
- The default `cowan-z` stencil removes vacuum phase error on the axial
  propagation direction; the standard Yee update is retained for regression.
- Relativistic particles use a Boris pusher with optional subcycling for
  analytical magnetic devices.
- Initial particle self-fields are obtained from a distributed CIC
  relativistic-Poisson solve that enforces the discrete Gauss constraint.
- Integer cell counts and physical cell sizes define the mesh. Physical
  extents and MPI slab offsets are constructed without floating-point cell
  counting.
- The inner CPML surface is the physical particle boundary. Escaping physical
  trajectories terminate there; compact output-free carriers preserve current
  continuity through the absorbing layer.

## Radiation workflow

The main radiation observable is a fixed laboratory field plane. A production
run stores the raw self-consistent Maxwell field; prescribed magnetic-device
fields are not copied into detector output.

Two field-analysis routes are available:

1. Direct analysis of the raw field plane, including any bound particle field.
2. Reconstruction and subtraction of a uniform-velocity particle background,
   using either a detector-owned ballistic-reference file or a colocated
   laboratory particle plane.

Background subtraction is performed on E/B amplitudes before Poynting power
or spectral energy is evaluated. The ballistic-reference route is valid only
when its left reference interval does not overlap a magnetic interaction
region. Full trajectories are an optional validation output and are not the
primary production record for large macro-particle counts.

The analysis programs are independent executables:

| Program | Function |
|---|---|
| `field_reconstruction` | Reconstruct and subtract the uniform-velocity particle field |
| `field_plane_analysis` | Compute angular spectra, energy spectra, Stokes quantities, and coherence |
| `trajectory_radiation` | Compute small-particle far-field radiation from laboratory trajectories |
| `field_power_compare` | Perform signal/baseline amplitude subtraction for the experimental retirement route |
| `energy_closure` | Compare particle-plane kinetic-energy change with field-plane radiation |
| `energy_ledger_report` | Summarize the runtime particle/field energy ledger |

## Input and output

- Configuration uses validated YAML mappings; unknown keys are rejected.
- Particle input uses versioned HDF5 records with SI coordinates, normalized
  proper velocity, stable particle identifiers, and positive relative macro
  weights.
- Each run creates a YAML manifest and embeds one `run_id`, configuration
  digest, source revision, and manifest path in all scientific HDF5 products.
- Existing scientific outputs are rejected unless `output.overwrite: true` is
  set explicitly.
- Detector and ledger files are written only by MPI rank 0. Optional trajectory
  output uses one file per rank.

## Runtime profiles

`runtime.mode: interactive` enables coordinated SIGINT/SIGTERM handling and
readable incomplete outputs for local testing. `runtime.mode: throughput`
removes signal polling and periodic durability flushes for scheduled cluster
runs. Disabled diagnostics allocate no buffers, open no files, and enter no
diagnostic communication paths.

The optional resource monitor reports estimated and measured resident memory,
output bounds, calibrated seconds per step, progress, and wall time in a format
suitable for batch-log capture.

## Validation requirements

The required regression suite verifies configuration loading, Lorentz
transforms, Elegant event reconstruction, Boris convergence, discrete charge
continuity, Cowan dispersion, CPML reflection, MPI consistency, detector and
stop behavior, HDF5 compatibility, output safety, both field-analysis routes,
and the runtime energy-ledger data path.

Passing the regression suite establishes software consistency, not convergence
of a physical problem. Production results require parameter studies covering
the field grid and time step, macro-particle count, particle substeps, CPML
thickness and reflection, transverse aperture, detector window, initial-field
boundary distance, and energy-ledger residual. The experimental particle
retirement route is disabled by default and requires a matched zero-radiation
baseline.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ -DBUILD_TESTING=ON
cmake --build build -j
cmake --build build --target verify_required
./build/aprl config/example.yaml
```

Parallel HDF5 is required by default for multi-rank particle input. Complete
dependency, installation, and cluster procedures are specified in the linked
manuals.

## Documentation

### Operation and file formats

- [Build, conversion, and execution](docs/BUILD_AND_RUN.md)
- [Installation, uninstall, and dependency maintenance](docs/INSTALLATION.md)
- [YAML input-card specification](docs/YAML_CONFIGURATION.md)
- [Particle HDF5 input specification](docs/PARTICLE_INPUT_HDF5.md)
- [Run identity and overwrite protection](docs/RUN_PROVENANCE.md)
- [Laboratory detector HDF5 output](docs/DETECTOR_OUTPUT_HDF5.md)
- [Laboratory trajectory HDF5 output](docs/TRAJECTORY_OUTPUT_HDF5.md)

### Numerical methods

- [Cowan-z Maxwell kernel and CPML boundary](docs/MAXWELL_COWAN_CPML.md)
- [Gauss-consistent initial particle self-field](docs/INITIAL_SELF_FIELD.md)
- [Particle subcycling](docs/PARTICLE_SUBCYCLING.md)
- [CPML-aware particle boundary](docs/PARTICLE_OPEN_BOUNDARY.md)
- [Field-detector ballistic reference region](docs/FIELD_DETECTOR_REFERENCE.md)

### Radiation and energy analysis

- [Particle-background field reconstruction](docs/FIELD_RECONSTRUCTION.md)
- [Field-plane spectrum and coherence analysis](docs/FIELD_PLANE_ANALYSIS.md)
- [Trajectory-to-far-field radiation](docs/TRAJECTORY_RADIATION.md)
- [Particle/field energy closure](docs/ENERGY_CLOSURE.md)
- [Runtime particle/field energy ledger](docs/ENERGY_LEDGER.md)
- [Laboratory-frame energy diagnostics](docs/LAB_FRAME_ENERGY_DIAGNOSTICS.md)
- [Experimental particle retirement and power comparison](docs/FIELD_POWER_COMPARISON.md)

### Verification and release status

- [Test framework](docs/TESTING.md)
- [Numerical validation record](docs/VALIDATION.md)
- [Release readiness and production requirements](docs/RELEASE_AUDIT.md)

Example inputs are provided in [`config/`](config/) and
[`postprocess/`](postprocess/).
