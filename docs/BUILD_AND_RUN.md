# Build, particle conversion, and execution

## Dependencies

- CMake 3.16 or newer;
- a C++11 compiler;
- MPI;
- HDF5 development libraries, with parallel HDF5 required for multi-rank
  particle input;
- yaml-cpp development libraries;
- the official SDDS C library/toolkit when building the direct Elegant SDDS
  converter;
- FFTW3, including its threads library, for field reconstruction and
  field-plane spectrum/coherence analysis.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

This produces `build/simulator`, `build/particle_text_to_hdf5`,
`build/elegant_sdds_to_hdf5` when SDDS is found,
`build/undulator_resonance`, `build/energy_ledger_report`, and
`build/lab_frame_energy_estimate`.

The SDDS converter is intentionally linked to the official implementation,
not a partial binary parser. If SDDS and Elegant were built in the recommended
sibling directory layout it is found automatically. Otherwise configure with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DSDDS_ROOT=/path/to/built/SDDS
```

Summarize a complete or cleanly interrupted runtime energy ledger without any
Python dependency:

```bash
./build/energy_ledger_report output/run-name/energy-ledger.h5
```

An optional second argument selects a zero-based committed record and reports
the initialization-to-that-sample change as well as the normal final summary.

It reads only the committed HDF5 prefix and prints the initial/final balance,
field fractions, mean gamma and energy-spread decomposition. See
[ENERGY_LEDGER.md](ENERGY_LEDGER.md).

Generate a separate laboratory-frame analytical scale estimate without
loading simulation particles or fields:

```bash
./build/lab_frame_energy_estimate config/lab_frame_energy_estimate.yaml
```

The example reports normalized 50 A bound-field and undulator-radiation
scales. See
[LAB_FRAME_ENERGY_DIAGNOSTICS.md](LAB_FRAME_ENERGY_DIAGNOSTICS.md).

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
one-dimensional gain-centre estimate. It also reports the planar-undulator
longitudinal estimate

```text
boost_gamma_recommended = gamma_lab/sqrt(1+K^2/2),
```

compares it with `mesh.boost_gamma`, and converts their velocity mismatch into
the predicted boosted-frame z drift across the characteristic undulator core,
in metres and z cells. The recommendation selects one constant inertial frame;
the simulator never changes boost during a run. Free-drift shift, fringes,
energy spread, emittance and collective evolution still require box margin.
If a field plane enables
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

The complete field-plane spectrum/coherence path is a fourth independent
build:

```bash
cmake -S postprocess/field_plane_analysis -B build-field-plane-analysis \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-plane-analysis -j
./build-field-plane-analysis/field_plane_analysis \
  postprocess/field_plane_analysis/example.yaml
```

It accepts raw detector fields, reconstructed particle-background-subtracted
fields, or a matched signal/baseline pair. Working-memory and output-size
guards are checked before output creation. See
[FIELD_PLANE_ANALYSIS.md](FIELD_PLANE_ANALYSIS.md).

The particle/radiation energy-closure report is another lightweight,
read-only build:

```bash
cmake -S postprocess/energy_closure -B build-energy-closure \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-energy-closure -j
./build-energy-closure/energy_closure \
  postprocess/energy_closure/example.yaml
```

It joins two particle planes by ID, performs stable high-gamma kinetic-energy
differencing, optionally subtracts a matched `K=0` result per particle, and
compares with either field-plane or trajectory far-field energy. See
[ENERGY_CLOSURE.md](ENERGY_CLOSURE.md).

## Convert particle input

For native Elegant output, prefer the direct converter:

```bash
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --page 1 \
  --output particles.h5
```

It reads `x,xp,y,yp,t,p` from one SDDS page through the official library,
preserves `particleID` when present, and writes HDF5 v4 fixed-plane crossing
events. No intermediate text file is produced. Use `--weight-column NAME`
only for a genuine positive per-row relative-weight column; the usual Elegant
page parameter `Charge` is not such a column. See
[PARTICLE_INPUT_HDF5.md](PARTICLE_INPUT_HDF5.md) for the exact mapping and
multi-page policy.

For an older six-column text export:

```bash
./build/particle_text_to_hdf5 \
  --input particles.txt \
  --output particles.h5 \
  --length-unit micrometer
```

The text converter reads the source twice and writes HDF5 in bounded chunks, so its
memory use does not grow with the full particle count. The HDF5-v3 text columns
are `x_plane y_plane zeta ux uy uz [macro_weight]`: x/y are sampled at the
fixed Elegant plane and `zeta` is the signed longitudinal offset used to
reconstruct the common-time bunch. The optional seventh value is a positive
relative macro-particle weight and defaults to one. Blank lines and `#`
comments are accepted. Set `beam.reference.input_plane_z` and
`initial_center_z` in the same laboratory coordinate system. The utility is
for preparation outside the expensive simulation allocation. Production runs
read the resulting HDF5 file directly and collectively, one contiguous range
per MPI rank.

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
does not fit the longitudinal box. For HDF5-v3/v4 input it also checks that every
record is projected forward from the configured Elegant plane and that the
reconstructed bunch fits between that plane and the first magnetic interaction
region. The head-anchored Lorentz synchronization span is reported separately
and is not treated as a physical entrance-drift requirement.

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
