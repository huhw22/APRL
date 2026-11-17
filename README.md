# APRL — Accelerator Particles and Radiation in Lorentz Frames

[中文说明](README.zh-CN.md)

The release is **APRL 2.0.0**, tagged `v2.0.0`. The default and preferred baseline is
**energy-conserving coupling with the longitudinal three-point kernel**
(`conservative-source + linear-binomial3`). See the [release record](docs/releases/2.0.0.md)
and [explicit configuration fragment](config/profiles/v2.0-three-point.yaml).

APRL is a C++/MPI program for accelerator-particle and radiation simulations
in Lorentz-boosted frames. The primary executable is `aprl`. Electromagnetic
evolution uses SI electric and magnetic fields directly; scalar and vector
potentials are not part of the time-advance state.

## Algorithms in 2.0

The installed solver accepts exactly two particle–field coupling algorithms:

| Choice | `mesh.longitudinal_shape` | Use |
|---|---|---|
| Energy-conserving three-point | `linear-binomial3` | Default; preferred baseline |
| Energy-conserving quintic compensation (historically “momentum matching”) | `quintic-compensated` | Alternative comparison; no guarantee of strict total momentum conservation |

Both use `particle_coupling: conservative-source` and `current_filter: none`.
Omitting the algorithm keys selects the three-point method. Explicit Boris,
plain linear, seven-point/notch, unknown methods and the separate experimental
current filter fail during configuration loading; there is no legacy-mode CLI
switch. Low-level reference kernels remain in the verification code.

```yaml
mesh:
  field_solver: cowan-z
  particle_coupling: conservative-source
  longitudinal_shape: linear-binomial3  # or quintic-compensated
  current_filter: none
  self_field_feedback: true
```

Both algorithms require the Cowan moving window, its matched reference supply,
held fixed-plane HDF5 input and self-field feedback. Legacy fixed-grid/Boris
cards and checkpoint/retirement routes are not supported production inputs.

These are algorithm settings, not a complete beamline card. Start from
[config/example.yaml](config/example.yaml), set the actual beam, mesh, magnets
and diagnostics, and use `aprl --version` to confirm the installed version.

## Validated applicability and limits

The release separates executable/software acceptance from physical accuracy.
The following conclusions summarize the accepted tests; the README deliberately
does not embed the complete figure archive.

| Test or use | Supported conclusion | Limit |
|---|---|---|
| Vacuum propagation, CPML, short drift and cold-charge benchmarks | The numerical and diagnostic paths have bounded reference/regression checks; cold spherical expansion uses the same relativistic shell theory sampled at each particle’s crossing event. | Accuracy still depends on grid, time step, domain, boundary and particle sampling. |
| Self-consistent bending radiation, two independent 100k unit-electron inputs | Three-point and quintic results agree closely with the prescribed-trajectory Liénard–Wiechert reference when the same finite window, coordinates, aperture and field definition are used. Three-point complete fixed-lab-tube residuals are about −0.262% and −0.313% in those tests. | Full near-field-inclusive intensity, radiation-only aperture intensity, axis spectrum and total-frequency energy balance are different observables; the reference is not the exact self-consistent trajectory. |
| Six-segment FEL with independent random input, w≈355 | The corrected loading removes the previously identified artificial startup enhancement from longitudinal phase clones. Both algorithms support this workflow and comparison with Puffin. | The complete fixed-lab energy-budget residual remains **26.51% (three-point)** and **37.56% (quintic)**. This release does not establish quantitative lab-energy closure or universal saturated-FEL accuracy. |
| Radiation diagnostics | Raw E/B can be compared with a reference including particle velocity fields; away from the beam, radiation-only agreement can be assessed separately. | There is no generally validated automatic near-field subtraction. A forward spectrum alone is not a complete energy ledger. |

Use independently sampled six-dimensional input for the random PIC baseline;
do not clone transverse phase-space points along z. Macro-particle weight has
no universal accepted upper bound: the w≈355 campaign is a specific resource
choice, not proof of preserved one-electron shot noise. Converge the quantity
you intend to report. Retain signed field ports, side flux, time caps, particle
transport and recorded source/numerical channels when checking energy; a closed
boost-frame discrete sum is not proof of laboratory physical closure.

See [release notes](docs/releases/2.0.0.md), [acceptance procedure](docs/RELEASE_ACCEPTANCE_2.0.md)
and [release audit](docs/RELEASE_AUDIT.md) for scope and reproduction.

## Scope

Version 2.0 retains and provides the following computation paths:

- self-consistent relativistic particle and Maxwell evolution;
- Cowan/CKC controlled-dispersion propagation along the configured z axis;
- compact CFS-CPML field boundaries and charge-continuous particle escape;
- analytical laboratory magnetic elements transformed into the simulation
  frame;
- HDF5 particle input, deterministic Gaussian test input, and native Elegant
  SDDS5 conversion through the official SDDS library;
- a standalone multiharmonic quiet-start/shot-noise bunch generator for SASE
  input, without adding a generator path to the simulation kernel;
- fixed laboratory field and particle detector planes;
- field-plane radiation, spectrum, polarization, coherence, and energy
  analysis on a separate analysis host;
- optional small-particle trajectory radiation and global energy accounting.

The general Cowan/CPML TF/SF laser or seed-field injection path is not
implemented in this release. Seeded-FEL and laser-modulation calculations that
require that boundary source are outside the validated scope.

## Numerical model

- E and B are advanced on a staggered Yee lattice in SI units.
- The default `cowan-z` stencil removes vacuum phase error on the axial
  propagation direction; the standard Yee update is retained for regression.
- The preferred coupling is `conservative-source` with the `linear-binomial3`
  longitudinal kernel. The alternative is `quintic-compensated`; retired Boris,
  linear and seven-point input choices are rejected. Analytical magnetic devices
  retain their prescribed-field substeps.
- HDF5-v4 fixed-lab-plane input offers legacy synchronization, frozen
  three-dimensional ballistic startup, and frozen held-transverse startup.
  For high-boost input where back-projection puts many particles outside a
  finite transverse mesh, held is the recommended tested starting choice,
  not a claim of radiation convergence; see the
  [startup configuration](docs/YAML_CONFIGURATION.md#beam-reference-and-input).
- Initial particle self-fields are obtained from a distributed CIC
  relativistic-Poisson solve that enforces the discrete Gauss constraint.
- Integer cell counts and physical cell sizes define the mesh. Physical
  extents and MPI slab offsets are constructed without floating-point cell
  counting.
- The inner CPML surface is the physical particle boundary. Escaping physical
  trajectories terminate there; compact output-free carriers preserve current
  continuity through the absorbing layer.
- Experimental `mesh.ballistic_outside_transverse: true` retains particles
  outside the transverse aperture as field-free ballistic particles and
  injects their clipped current trajectories on re-entry. The initial Poisson
  source normally uses particles inside the physical region; an optional
  one-time solve also includes particles inside the transverse CPML.
  This finite-window approximation does not guarantee global Gauss/energy
  closure; see [boundary-injection limits](docs/PARTICLE_OPEN_BOUNDARY.md).

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

On a networked host, APRL, parallel HDF5 with zlib, and yaml-cpp can be
installed into one user-owned prefix. A managed cluster reuses its site MPI
by default:

```bash
./install_online.sh --profile full
```

Large online-install archives use the official CMake, zlib, and HDF Group
download sites instead of GitHub releases. Every archive is pinned by SHA-256
and cached for reuse.

The installer verifies pinned source archives, builds parallel HDF5 with the
selected MPI, and gives installed executables relative runtime search paths.
Only one `aprl-runtime` tree is needed on a shared filesystem.  A single-node
workstation without a platform MPI can use the lower-level
`./install.sh --bootstrap-deps --bundle-mpi` path.

The direct build remains available when all dependencies already exist:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ -DBUILD_TESTING=ON
cmake --build build -j
cmake --build build --target verify_required
mkdir demo-2.0
cp config/example.yaml demo-2.0/run.yaml
./build/uniform_lab_plane_bunch_generator --output demo-2.0/particles.h5 \
  --length-m 1e-6 --current-A 0.001 --central-momentum-MeV-c 2 \
  --emit-x-m 1e-12 --emit-y-m 1e-12 --beta-x-m 1 --beta-y-m 1 --sigma-dp 0
(cd demo-2.0 && ../build/aprl run.yaml)
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
- [Standalone SASE bunch generator](docs/SASE_BUNCH_GENERATOR.md)
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
