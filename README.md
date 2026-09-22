# Unnamed direct E/B FEL simulator

The program deliberately has no final public name yet. `simulator` is only a
temporary executable name and may be replaced without defining the project's
identity.

This is a standalone C++/MPI program for boosted-frame free-electron-laser
simulation. It advances electric and magnetic fields directly on a Yee lattice
in SI units, pushes relativistic particles with the Boris method, deposits
charge-conserving trajectory current, and writes compact laboratory-frame
particle histories for radiation post-processing.

The repository is independent: it has its own Git history and contains no
A/phi solver, predecessor solver selector, predecessor input parser,
compatibility output path, or compatibility adapter. The retained physical
idea is the explicit Lorentz
boost between the laboratory frame and the moving computational box.

## Current scope

- direct SI E/B leapfrog Maxwell update;
- z-slab MPI field exchange and particle migration;
- relativistic Boris particle push;
- charge-conserving current deposition;
- TF/SF injection for seed, modulation, and other incident waves;
- prescribed laboratory-frame undulators and dipoles;
- explicit bunch placement followed by free-drift simultaneity conversion;
- per-rank buffered HDF5 particle trajectories;
- interactive recoverable output and low-overhead HPC throughput output.

The current boundary is PEC. A Gauss-consistent initial particle self-field,
particle subcycling, and CPML remain future work. Until those are implemented
and benchmarked, results should be treated as integration/validation results,
not final radiation-production data.

## Dependencies

- CMake 3.16 or newer
- C++11 compiler
- MPI
- HDF5 development files
- yaml-cpp development files

On Ubuntu/Debian, the YAML dependency is `libyaml-cpp-dev`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

## Run

Single process:

```bash
./build/simulator config/example.yaml
```

MPI:

```bash
mpirun -n 4 ./build/simulator config/example.yaml
```

The complete input schema and source roles are described in
`docs/YAML_CONFIGURATION.md`. Generated HDF5 files, build trees, scheduler
logs, development tests, probes, benchmarks, and helper tools are excluded by
`.gitignore` and are not part of the source history.
