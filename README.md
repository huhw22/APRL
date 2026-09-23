# Unnamed boosted-frame FEL program

The program intentionally has no final public name yet. `simulator` is only a
temporary executable name.

## Basic idea

This is a lightweight C++/MPI code for free-electron-laser and related
relativistic beam simulations. Electric and magnetic fields are advanced
directly in SI units on a staggered Yee lattice. The primary Maxwell kernel is
a Cowan/CKC controlled-dispersion stencil with exact axial vacuum propagation
for the configured z direction; the standard Yee update remains as a
regression option. Particles are advanced self-consistently with those fields
while laboratory seed fields and magnetic devices are converted through an
explicit Lorentz boost. A/phi is not part of the evolution state.

The main reusable result is a compact, laboratory-frame particle history.
Optional fixed laboratory field and particle detector planes provide focused
measurements without forcing a large multidimensional field dump during every
supercomputer run.  The outputs are intended for radiation reconstruction and
other analysis on a separate, less expensive machine.

An independent C++/MPI post-processor reconstructs the complex polarized
far-field spectrum directly from those histories. It provides angular and
integrated energy spectra, Stokes data, and optional ensemble cross-spectral
density without linking or rerunning the simulation core. Stable parts of a
single pulse can alternatively be treated as a Hann-window ensemble in reduced
observer time for spatial and two-frequency coherence analysis.

## Intended applications

- seeded FEL and laser-modulation studies;
- boosted-frame electron motion through undulators and other magnetic devices;
- scalable trajectory production for radiation post-processing;
- numerical experiments on a direct SI E/B Maxwell-particle formulation.

Runs can end either after all still-valid particles pass the final element's
finite interaction region, or when the laboratory boost-reference centre
reaches a configured downstream z coordinate. Laboratory field and particle
detector planes are zero-length elements in this same ordering model. Only
MPI rank zero writes their HDF5 files; disabling detectors constructs no
detector object and enters no detector communication.

The runtime policy is selected independently of outputs: `interactive` is the
small-server test path with coordinated clean signal stopping, while
`throughput` removes signal polling and periodic durability flushes for
scheduled supercomputer runs.

The implementation is still a development solver. Compact CFS-CPML is
available for no-seed Cowan runs, while the generalized Cowan/CPML TF/SF
seed-wave correction, a Gauss-consistent initial particle self-field, and
particle subcycling remain planned work. Production radiation results still
require problem-scale convergence and reflection validation.

The inner CPML surface is also the physical particle boundary. A physical
trajectory ends exactly there and becomes a compact, output-free ballistic
carrier whose current is damped with the CPML conductivity profile. The outer
face performs charge-conserving residual cleanup. This avoids treating the
absorber as an observation region without creating the discontinuity caused by
immediate particle deletion. The production path uses no dense boundary arrays
or particle-boundary I/O.

## Documentation

- [Build, conversion, and execution](docs/BUILD_AND_RUN.md)
- [YAML input-card specification](docs/YAML_CONFIGURATION.md)
- [Cowan-z kernel and CPML boundary](docs/MAXWELL_COWAN_CPML.md)
- [CPML-aware particle boundary](docs/PARTICLE_OPEN_BOUNDARY.md)
- [Particle HDF5 file specification](docs/PARTICLE_INPUT_HDF5.md)
- [Laboratory trajectory HDF5 output](docs/TRAJECTORY_OUTPUT_HDF5.md)
- [Trajectory-to-far-field radiation tool](docs/TRAJECTORY_RADIATION.md)
- [Laboratory detector HDF5 output](docs/DETECTOR_OUTPUT_HDF5.md)
- [HDF5 input example](config/example.yaml)
- [Generated Gaussian test example](config/generated_gaussian.yaml)
