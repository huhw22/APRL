# Build, particle conversion, and execution

## Dependencies

- CMake 3.16 or newer;
- a C++11 compiler;
- MPI;
- HDF5 development libraries, with parallel HDF5 required for multi-rank
  particle input;
- yaml-cpp development libraries.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

This produces `build/simulator` and `build/particle_text_to_hdf5`.

## Convert a legacy text particle file

```bash
./build/particle_text_to_hdf5 \
  --input particles.txt \
  --output particles.h5 \
  --length-unit micrometer
```

The converter reads the source twice and writes HDF5 in bounded chunks, so its
memory use does not grow with the full particle count. The six text columns are
`x y z ux uy uz`; blank lines and `#` comments are accepted. The utility is for
preparation outside the expensive simulation allocation. Production runs read
the resulting HDF5 file directly and collectively, one contiguous range per
MPI rank.

The committed example HDF5 can be regenerated with:

```bash
./build/particle_text_to_hdf5 \
  --input examples/particles/example_particles.txt \
  --output examples/particles/example_particles.h5 \
  --length-unit micrometer
```

## Run

```bash
./build/simulator config/example.yaml
mpirun -n 4 ./build/simulator config/example.yaml
```

Use `config/generated_gaussian.yaml` for the file-free test input. For a local
debug run choose global `runtime.mode: interactive`; for an uninterrupted
batch run choose `runtime.mode: throughput`. The choice is global, so local
signal stopping still works when trajectory output is disabled and only a
detector plane is active.

The production example also enables one laboratory field plane and one
particle plane. Their two HDF5 files appear under
`output/example/detectors/`; MPI rank zero is their only writer. Remove the
`detectors` block or set `detectors.enabled: false` for runs that need no
detector allocation or communication. Detector datasets are described in
[DETECTOR_OUTPUT_HDF5.md](DETECTOR_OUTPUT_HDF5.md).

The current solver stops with an explicit error if its boosted initial bunch
does not fit the longitudinal box or overlaps the finite interaction region of
the first element. Increase the box or move `initial_center_z` upstream only
after checking the reported physical and interaction boundaries.

The examples select `mesh.field_solver: cowan-z`. This requires `dx >= dz`
and `dy >= dz`; an invalid card exits before field allocation and particle
input and prints a transverse-grid recommendation. The startup log records the
Cowan coefficients and transverse vacuum dispersion diagnostics. CPML is not
yet connected. A nonempty TF/SF `incident_waves` list currently requires
`mesh.field_solver: yee` until the incident-boundary stencil is generalized.
