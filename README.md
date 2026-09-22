# Unnamed boosted-frame FEL program

The program intentionally has no final public name yet. `simulator` is only a
temporary executable name.

## Basic idea

This is a lightweight C++/MPI code for free-electron-laser and related
relativistic beam simulations. Electric and magnetic fields are advanced
directly in SI units on a Yee mesh; particles are advanced self-consistently
with those fields while laboratory seed fields and magnetic devices are
converted through an explicit Lorentz boost. A/phi is not part of the
evolution state.

The main reusable result is a compact, laboratory-frame particle history. It
is intended to retain enough information for radiation reconstruction and
other analysis on a separate, less expensive machine, rather than forcing a
large multidimensional field dump during every supercomputer run.

## Intended applications

- seeded FEL and laser-modulation studies;
- boosted-frame electron motion through undulators and other magnetic devices;
- scalable trajectory production for radiation post-processing;
- numerical experiments on a direct SI E/B Maxwell-particle formulation.

The implementation is still a development solver. CPML, a Gauss-consistent
initial particle self-field, and particle subcycling remain planned work, so
production radiation results require further validation.

## Documentation

- [Build, conversion, and execution](docs/BUILD_AND_RUN.md)
- [YAML input-card specification](docs/YAML_CONFIGURATION.md)
- [Particle HDF5 file specification](docs/PARTICLE_INPUT_HDF5.md)
- [HDF5 input example](config/example.yaml)
- [Generated Gaussian test example](config/generated_gaussian.yaml)
