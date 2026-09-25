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
The mesh input treats integer cell counts and physical cell sizes as the only
authoritative geometry. Full extents are obtained by multiplication, and MPI z
slabs use integer quotient/remainder offsets, avoiding a floating length/spacing
division when constructing unequal-rank partitions.


The main radiation result is a fixed laboratory field plane. Each field plane
can own a diagnostic-only ballistic-reference region on its left: particles
continue their fully self-consistent push, while a single crossing record
supplies the straight-line charged-particle background reference for later
subtraction.
Small test runs can additionally compare that virtual straight line with the
same physical particle at the field plane, one paired record per particle,
without enabling full trajectory output.
Particle detector planes remain independent. Full laboratory particle
histories are retained only as an optional small-particle/debug path because
their production-scale storage cost is prohibitive. Outputs are intended for
analysis on a separate, less expensive machine.

For small validation runs, an independent C++/MPI post-processor reconstructs
the complex polarized far-field spectrum directly from trajectories. It
provides angular and integrated energy spectra, Stokes data, and optional
ensemble cross-spectral density without linking or rerunning the simulation
core. Stable parts of a single pulse can alternatively be treated as a
Hann-window ensemble in reduced observer time for spatial and two-frequency
coherence analysis.

A second standalone C++ post-processor combines a raw laboratory E/B plane
with either its ballistic-reference crossings or a colocated particle plane.
It deposits the crossings once, reconstructs the uniform-velocity bunch field
with a three-dimensional FFT, writes a background-subtracted field, and
reports the change in forward Poynting power and energy. This keeps the
production simulation light while retaining a reproducible field-based
radiation path on the analysis server.

The cleaned or raw laboratory plane can then be passed to a threaded C++/FFTW
analysis tool. It writes the forward propagating angular energy spectrum,
Stokes data, integrated spectrum, coherent/fluctuation split, selected spatial
cross-spectral densities, a global transverse-coherence/Gram-matrix diagnostic
and selected two-frequency coherence. Stable portions of one run can be
divided into energy-normalized Hann windows, keeping this large-run radiation
path independent of trajectory storage.

A lightweight energy-closure diagnostic joins two laboratory particle planes
by particle ID and compares their stable per-particle kinetic-energy loss with
the forward field-plane radiation energy. It supports a matched zero-radiation
run and subtracts that baseline per particle before summation, so high-gamma
roundoff is not amplified by subtracting two total beam energies.

An alternative, explicitly experimental radiation path ends physical particle
push/output after the magnetic system and tapers a compact ballistic current
carrier to zero before the field plane. It does not change the Maxwell kernel
and is disabled by default. Because removing net charge inside the domain is
not continuity exact, every such run must be paired with a zero-radiation
baseline. A third C++ post-processor subtracts that baseline at the E/B-amplitude
level, reports instantaneous forward power and band energy, and can compare
both peak power and energy with the trajectory far field for small tests.

The retirement route can optionally bind one terminal field detector to a
manual frequency-protection guard. Startup uses the actual maximum laboratory
particle gamma to verify the requested C2 transition cycles and reserves the
same conservative gamma-times-transverse-diagonal causal distance used by the
straight-line diagnostic. A small read-only utility estimates the
characteristic planar-undulator resonance before the user chooses the protected
band; it never changes the input card.

## Intended applications

- seeded FEL and laser-modulation studies;
- boosted-frame electron motion through undulators and other magnetic devices;
- scalable laboratory field-plane production for radiation analysis;
- numerical experiments on a direct SI E/B Maxwell-particle formulation.

Runs can end either after all still-valid particles pass the final element's
finite interaction region, or when the laboratory boost-reference centre
reaches a configured downstream z coordinate. Laboratory field and particle
detectors participate in this same ordering model. A field plane has a
left-only diagnostic interaction extent; a particle plane remains
geometrically zero-length. Only MPI rank zero writes their HDF5 files;
disabling detectors constructs no detector object and enters no detector
communication.

The runtime policy is selected independently of outputs: `interactive` is the
small-server test path with coordinated clean signal stopping, while
`throughput` removes signal polling and periodic durability flushes for
scheduled supercomputer runs.

Production HDF5 input supports a positive relative weight per particle and
normalizes those weights to the total electron count in the YAML card. The
generated Gaussian path instead uses the requested total charge and macro-
particle count with uniform weights. An optional root-only resource report
gives a pre-run memory/time estimate, periodic batch-log progress, and measured
peak resident memory and wall time without requiring a live terminal.

An independent optional runtime energy ledger closes the boosted-frame budget
across active and escaped particle kinetic energy, physical-interior E/B
energy, six inner-CPML Poynting fluxes and prescribed-device work. It also
tracks mean laboratory gamma, projected energy spread, linear chirp and the
linearly detrended spread. Disabled mode adds no loops, collectives or files;
enabled mode writes one small rank-zero HDF5 stream and has a C++ summary tool.


The implementation is still a development solver. Compact CFS-CPML is
available for no-seed Cowan runs, and a temporary distributed
CIC/relativistic-Poisson
solve now gives the bunch a relativistic-Poisson, Gauss-consistent initial E/B
field without retaining A/phi. Particle Boris substeps resolve analytical
laboratory devices without raising the Maxwell cadence; grid E/B sampling,
current deposition and field
detectors intentionally remain on the field step. The generalized Cowan/CPML
TF/SF seed-wave correction remains planned work. Production radiation results
still require problem-scale convergence, initial-field distance-to-CPML
convergence, reflection validation, and convergence of the runtime energy
ledger. The first dense-bunch ledger closes the driven boosted-frame exchange
to about four percent, which validates the accounting path but is not yet a
production tolerance.
A controlled one-electron-equivalent test gives the correct loss sign and the
same energy scale as its finite-band forward radiation, but the `10^6`-electron
test is dominated by changing collective/bound-field energy. One downstream
plane alone cannot close that ledger; stored E/B energy and all boundary fluxes
remain a production gate.

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
- [Gauss-consistent initial particle self-field](docs/INITIAL_SELF_FIELD.md)
- [Laboratory-frame energy diagnostics and 50 A scale estimate](docs/LAB_FRAME_ENERGY_DIAGNOSTICS.md)
- [Particle subcycling and its field-step limits](docs/PARTICLE_SUBCYCLING.md)
- [CPML-aware particle boundary](docs/PARTICLE_OPEN_BOUNDARY.md)
- [Particle HDF5 file specification](docs/PARTICLE_INPUT_HDF5.md)
- [Laboratory trajectory HDF5 output](docs/TRAJECTORY_OUTPUT_HDF5.md)
- [Trajectory-to-far-field radiation tool](docs/TRAJECTORY_RADIATION.md)
- [Laboratory detector HDF5 output](docs/DETECTOR_OUTPUT_HDF5.md)
- [Field-detector ballistic reference region](docs/FIELD_DETECTOR_REFERENCE.md)
- [Particle-background field reconstruction](docs/FIELD_RECONSTRUCTION.md)
- [Field-plane spectrum and coherence analysis](docs/FIELD_PLANE_ANALYSIS.md)
- [Particle/field energy-closure diagnostic](docs/ENERGY_CLOSURE.md)
- [Runtime particle/field energy ledger](docs/ENERGY_LEDGER.md)
- [Particle retirement and matched power comparison](docs/FIELD_POWER_COMPARISON.md)
- [Numerical validation status](docs/VALIDATION.md)
- [Current release audit and production gates](docs/RELEASE_AUDIT.md)
- [HDF5 input example](config/example.yaml)
- [Generated Gaussian test example](config/generated_gaussian.yaml)
