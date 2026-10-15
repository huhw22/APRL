# Build, particle conversion, and execution

## Dependencies

- CMake 3.16 or newer for APRL itself; the online dependency bootstrap uses
  CMake 3.18 or supplies a private newer version for HDF5 1.14.6;
- a C++11 compiler;
- MPI;
- HDF5 1.10 or newer development libraries, with parallel HDF5 required for multi-rank
  particle input;
- yaml-cpp 0.6 or newer development libraries;
- the official SDDS C library/toolkit when building the direct Elegant SDDS
  converter;
- FFTW3, including its threads library, for field reconstruction and
  field-plane spectrum/coherence analysis.

Use the compiled dependency probe before building on a new machine or after
changing compiler/MPI modules. The complete installation, offline-cluster and
uninstall procedures are in [INSTALLATION.md](INSTALLATION.md).

## Build

On a new machine, prefer the unprivileged local dependency bootstrap.  It
reuses the selected MPI and places zlib-enabled parallel HDF5, yaml-cpp, and
APRL in one prefix:

```bash
./install_online.sh --profile core
```

The online installer rejects Conda MPI on a managed cluster, can load a site
MPI with `--mpi-module`, and skips the two-rank launch on a login node by
default. Repeat the MPI preflight inside an allocation before production.
The lower-level installer and `--bundle-mpi` single-node fallback remain
available; production multi-node runs should keep the site MPI. Cache,
offline, and runtime-library details are in [installation and dependency
maintenance](INSTALLATION.md).

When every system dependency is already present, configure directly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

The main build now prefers and, by default, requires parallel HDF5. A strictly
single-rank local build may opt out with
`-DAPRL_REQUIRE_PARALLEL_HDF5=OFF`; that binary cannot read particle HDF5 input
with multiple MPI ranks.

This produces `build/aprl`, `build/particle_text_to_hdf5`,
`build/sase_bunch_generator`,
`build/elegant_sdds_to_hdf5` when SDDS is found,
`build/undulator_resonance`, `build/energy_ledger_report`, and
`build/lab_frame_energy_estimate`.

The SASE input generator is a standalone preprocessor and may be omitted with
`-DAPRL_BUILD_SASE_BUNCH_GENERATOR=OFF`. Its Fawley quiet-start/shot-noise
model, command line, diagnostics, and applicability limits are documented in
[SASE_BUNCH_GENERATOR.md](SASE_BUNCH_GENERATOR.md).

Before a production run, execute the lightweight required regression suite:

```bash
cmake --build build --target verify_required
```

It covers startup plus the Lorentz/Elegant, Boris, charge-continuity, Cowan,
CPML, MPI, detector/stop, HDF5-compatibility and K=0 ledger paths. Optional
non-gating physical energy-budget examples are kept separate. See
[TESTING.md](TESTING.md) for the test matrix and interpretation.

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
APRL never changes boost during a run. Free-drift shift, fringes,
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

The reference example HDF5 can be regenerated with:

```bash
./build/particle_text_to_hdf5 \
  --input examples/particles/example_particles.txt \
  --output examples/particles/example_particles.h5 \
  --length-unit micrometer
```

## Run

For the recommended project-local installation:

```bash
module load <the same site MPI module used for installation>
source ./aprl-runtime/share/aprl/runtime-env.sh
aprl config/example.yaml
mpirun -n 4 aprl config/example.yaml
```

A developer who deliberately retained a direct CMake build can substitute
`./build/aprl` for `aprl`.

Use `config/generated_gaussian.yaml` for the file-free test input. For a local
debug run choose global `runtime.mode: interactive`; for an uninterrupted
batch run choose `runtime.mode: throughput`. The choice is global, so local
signal stopping still works when trajectory output is disabled and only a
detector plane is active.

The reference cards enable the root-only resource report. Before the time
loop it prints an estimated per-rank peak resident set, aggregate modeled
memory, uncompressed output upper bounds, calibrated seconds per step, and a
wall-time upper estimate. During a long run it prints one compact progress
record at the configured interval, and shutdown prints measured current/peak
resident memory and loop/full wall time. Every record is a flushed single line
prefixed with `[resource]`, so ordinary batch redirection is sufficient:

```bash
sbatch --output=run-%j.log run.sh
# run.sh ultimately executes, for example:
module load <the same site MPI module used for installation>
source /shared/path/APRL/aprl-runtime/share/aprl/runtime-env.sh
mpirun aprl config/production.yaml
```

The startup timing is a short local microbenchmark, not a scheduler guarantee:
parallel-file-system contention, later particle migration and early physical
stopping are not predicted. Calibrate the memory/time safety factors on the
intended machine and MPI decomposition before requesting a large allocation.

A successful process shutdown prints flushed `[teardown]` markers around the
final stop synchronization and energy sample, trajectory/detector/ledger
close, particle-boundary and resource reports, `Simulation` destruction, and
`MPI_Finalize`. Each line
contains the cached MPI rank, host name and process ID, a monotonic elapsed
time, and the phase. Rank zero prints them by default; for a diagnostic rerun,
export
`APRL_TEARDOWN_TRACE=all` to print the same markers from every rank. The
variable changes logging only: in particular, tracing adds no barrier or other
MPI synchronization. The marker after `MPI_Finalize` uses only cached process
identity and the local monotonic clock; it makes no MPI call.

Independently of tracing, the default shutdown watchdog gives each rank 600
seconds for the final MPI reductions, `Simulation` destruction, and
`MPI_Finalize`/`MPI_Abort`. It is deliberately inactive during physical
stepping and the normal explicit HDF5 close phase. Exceptional stack unwinding
is guarded as a whole, including any destructor fallback cleanup. APRL requests
`MPI_THREAD_FUNNELED`; the one sleeping
watchdog thread per rank never calls MPI or performs numerical work. A timeout
exits nonzero with status 124. With `APRL_TEARDOWN_TRACE=all`, each rank's
last flushed `[teardown]` enter marker identifies the guarded phase, and a
batch script cannot mistake a cleanup failure for success. Configure it under
`runtime.shutdown_watchdog`; the resolved setting is also stored in the run
manifest.

For an all-rank trace, interpret each rank's last phase as follows:

- the `*-enter` and `*-complete` pairs identify the exact final solve operation
  which did not return; `solve-complete` proves that `Simulation::solve()`
  reached its normal return path;
- `simulation-destructor-enter` means that rank did not finish destruction;
- `simulation-destructor-complete` followed by `mpi-finalize-enter`, but not
  `mpi-finalize-complete`, identifies a rank stuck in MPI finalization;
- `mpi-finalize-complete` from every rank places the stall after MPI. If no
  APRL processes remain while the scheduler step is still active, the launcher
  or scheduler cleanup is stuck. If APRL processes remain, inspect them for the
  much narrower post-finalize C/C++ process-exit path.


The compact I/O example enables one laboratory field plane and one particle
plane, but disables the ballistic background reference because its tiny box
has no room for the required free drift. Their two HDF5 files appear under
`output/example/detectors/`; MPI rank zero is their only writer. Remove the
`detectors` block or set `detectors.enabled: false` for runs that need no
detector allocation or communication. Detector datasets are described in
[DETECTOR_OUTPUT_HDF5.md](DETECTOR_OUTPUT_HDF5.md).

Production radiation cards that use straight-line particle-background
reconstruction should explicitly enable the field plane's
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
record is projected forward from the configured Elegant plane. When magnetic
devices exist, the reconstructed bunch must additionally fit before the first
magnetic interaction region; detector-only and element-free cards do not invent
that constraint. The head-anchored Lorentz synchronization span is reported
separately and is not treated as a physical entrance-drift requirement.

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
