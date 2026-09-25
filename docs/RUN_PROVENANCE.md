# Run identity, manifests and overwrite safety

Every APRL invocation creates one root-only YAML manifest before field
allocation and gives the run a unique `run_id`. The manifest records the exact
input card, its FNV-1a configuration digest, the source revision (including a
`-dirty` suffix when configured from a modified tree), compiler/build/HDF5
description, MPI implementation and rank count, and particle-input identity.
For a large HDF5 particle file the identity uses path, byte size and mtime; it
does not reread and hash the entire file at startup.

The same `run_id`, configuration digest, source revision and manifest path are
written into every APRL-produced field-plane, particle-plane, ballistic-
reference, trajectory and energy-ledger HDF5 group. These attributes are the
authoritative way to check that separately stored files came from one run.

```yaml
output:
  overwrite: false
  manifest: output/run-manifest.yaml
```

`overwrite` defaults to `false`. APRL checks every enabled output
path on all MPI ranks before starting and refuses an existing manifest, rank
trajectory, detector or ledger file. Set it to `true` only when intentional
replacement of the complete run is acceptable. It applies to all outputs from
that invocation, so mixed old/new files cannot be authorized accidentally one
at a time.

Every post-processing card has its own `output.overwrite`, also defaulting to
`false`. Post-processors reject unknown YAML keys at every mapping level and
use exclusive HDF5 creation when replacement is disabled. Use a new output
name for normal analysis iterations; explicit overwrite is intended mainly
for disposable tests.

The manifest is written only by rank zero. Detector and ledger HDF5 remain
root-only, while trajectory files remain one file per rank, so this mechanism
does not introduce shared MPI writes.
