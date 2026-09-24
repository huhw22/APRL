# Build, particle conversion, and execution

## Dependencies

- CMake 3.16 or newer;
- a C++11 compiler;
- MPI;
- HDF5 development libraries, with parallel HDF5 required for multi-rank
  particle input;
- yaml-cpp development libraries;
- FFTW3, including its threads library, for field reconstruction only.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

This produces `build/simulator`, `build/particle_text_to_hdf5`, and
`build/undulator_resonance`.

Before choosing a retirement protection band, inspect the characteristic
planar-undulator resonance:

```bash
./build/undulator_resonance config/generated_gaussian.yaml
./build/undulator_resonance config/example.yaml 1000
```

The second argument supplies laboratory gamma when an HDF5 beam card has no
`beam.input.gamma`. With one planar undulator the tool selects it
automatically; with several, mark exactly one magnetic element
`characteristic: true`. It reports the on-axis cold-beam fundamental as a
one-dimensional gain-centre estimate. If a field plane enables
`retirement_frequency_protection`, the same command reads its minimum photon
energy and cycle count, reports the required manually configured retirement
length, and prints a preflight pass/fail result. It never edits the card.

The trajectory radiation post-processor is deliberately a separate build:

```bash
cmake -S postprocess/trajectory_radiation -B build-radiation \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-radiation -j
```

It reads completed per-rank trajectory files, reassembles migrated particle
histories, and writes one rank-zero HDF5 result. See
[TRAJECTORY_RADIATION.md](TRAJECTORY_RADIATION.md).

The raw-field/particle-plane reconstruction tool is also a separate build:

```bash
cmake -S postprocess/field_reconstruction -B build-field-reconstruction \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-reconstruction -j
./build-field-reconstruction/field_reconstruction \
  postprocess/field_reconstruction/example.yaml
```

It uses threaded FFTW on one analysis node and never modifies the simulation
outputs. Its input assumptions, output datasets, intensity diagnostics, and
convergence requirements are documented in
[FIELD_RECONSTRUCTION.md](FIELD_RECONSTRUCTION.md).

The matched zero-radiation power comparison is a third independent build:

```bash
cmake -S postprocess/field_power_compare -B build-field-power-compare \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-power-compare -j
./build-field-power-compare/field_power_compare \
  postprocess/field_power_compare/example.yaml
```

It streams small transverse batches from two matching detector files, performs
field-amplitude subtraction and band filtering, and can compare peak power and
energy with a one-shot trajectory far field. See
[FIELD_POWER_COMPARISON.md](FIELD_POWER_COMPARISON.md).

## Convert a legacy text particle file

```bash
./build/particle_text_to_hdf5 \
  --input particles.txt \
  --output particles.h5 \
  --length-unit micrometer
```

The converter reads the source twice and writes HDF5 in bounded chunks, so its
memory use does not grow with the full particle count. The text columns are
`x y z ux uy uz [macro_weight]`; the optional seventh value is a positive
relative macro-particle weight and defaults to one. Blank lines and `#`
comments are accepted. The utility is for preparation outside the expensive
simulation allocation. Production runs read the resulting HDF5 file directly
and collectively, one contiguous range per MPI rank.

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

The committed cards enable the root-only resource report. Before the time
loop it prints an estimated per-rank peak resident set, aggregate modeled
memory, uncompressed output upper bounds, calibrated seconds per step, and a
wall-time upper estimate. During a long run it prints one compact progress
record at the configured interval, and shutdown prints measured current/peak
resident memory and loop/full wall time. Every record is a flushed single line
prefixed with `[resource]`, so ordinary batch redirection is sufficient:

```bash
sbatch --output=run-%j.log run.sh
# run.sh ultimately executes, for example:
srun ./build/simulator config/production.yaml
```

The startup timing is a short local microbenchmark, not a scheduler guarantee:
parallel-file-system contention, later particle migration and early physical
stopping are not predicted. Calibrate the memory/time safety factors on the
intended machine and MPI decomposition before requesting a large allocation.


The compact I/O example enables one laboratory field plane and one particle
plane, but disables the ballistic background reference because its tiny box
has no room for the required free drift. Their two HDF5 files appear under
`output/example/detectors/`; MPI rank zero is their only writer. Remove the
`detectors` block or set `detectors.enabled: false` for runs that need no
detector allocation or communication. Detector datasets are described in
[DETECTOR_OUTPUT_HDF5.md](DETECTOR_OUTPUT_HDF5.md).

Production radiation cards should enable the field plane's
`particle_background` block. Startup computes its left reference boundary
from the transverse mesh diagonal and maximum input-particle laboratory gamma,
then rejects overlap with magnetic interaction regions and prints a corrected
detector position. See
[FIELD_DETECTOR_REFERENCE.md](FIELD_DETECTOR_REFERENCE.md).

For a small validation card, enable
`particle_background.validation.enabled`. The run then records only the two
crossings needed to compare the physical orbit with the virtual straight line;
it does not require full trajectory output. A mandatory particle-count guard
prevents this diagnostic from being enabled accidentally on a production
bunch.

The current solver stops with an explicit error if its boosted initial bunch
does not fit the longitudinal box or overlaps the finite interaction region of
the first element. Increase the box or move `initial_center_z` upstream only
after checking the reported physical and interaction boundaries.

The examples select `mesh.field_solver: cowan-z`. This requires `dx >= dz`
and `dy >= dz`; an invalid card exits before field allocation and particle
input and prints a transverse-grid recommendation. The startup log records the
Cowan coefficients and transverse vacuum dispersion diagnostics. They also
select compact CFS-CPML and report its actual auxiliary memory. A nonempty
TF/SF `incident_waves` list is currently rejected with CPML, and Cowan TF/SF
remains unavailable until the injection surface is generalized. See
[MAXWELL_COWAN_CPML.md](MAXWELL_COWAN_CPML.md) for the numerical details and
current MPI restriction. Physical particle histories stop at the inner CPML
surface and continue only as compact, output-free carriers with matched current
damping. The outer face performs charge-conserving residual cleanup, with no
production face-array allocation or boundary I/O; see
[PARTICLE_OPEN_BOUNDARY.md](PARTICLE_OPEN_BOUNDARY.md). The trajectory v2
schema and exact CPML-entry events are documented in
[TRAJECTORY_OUTPUT_HDF5.md](TRAJECTORY_OUTPUT_HDF5.md).
