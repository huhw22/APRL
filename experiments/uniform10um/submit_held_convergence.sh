#!/usr/bin/env bash
# Submit held-field convergence cards using an existing particle HDF5 file.
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./submit_held_convergence.sh --reference HELD_RUN_DIR \
         --input PARTICLES_H5 --output-dir NEW_SUITE_DIR --electrons N \
         [--variants wide,thick-cpml,fine-z] [--plan]

The input particle file is reused in place: it is never copied or regenerated.
HELD_RUN_DIR is the existing g929-1m-held-field run directory (with a ledger).
--plan only checks inputs and prints the intended Slurm resources; it writes nothing.

Environment: APRL_SRC_DIR, APRL_BIN, LEDGER_REPORT_BIN, PYTHON_BIN,
  PARTITION=ssmb, NODES=1, TASKS_PER_NODE=8, MPI_RANKS=NODES*TASKS_PER_NODE,
  MPI_LAUNCHER=mpiexec|srun, SRUN_MPI=pmix_v3, WALLTIME=144:0:0,
  MODULES, NODELIST, EXCLUDE_NODES, SBATCH_ACCOUNT, SBATCH_QOS.

IMPORTANT: wide and fine-z increase grid size and memory. The 1-node/8-rank
default is only a template, not a resource recommendation. Production runs
previously used 2 nodes / 128 MPI ranks; preflight memory and set NODES,
TASKS_PER_NODE, and MPI_RANKS for your cluster before submitting.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

reference_arg=
input_arg=
output_arg=
electrons=
variants=wide,thick-cpml,fine-z
plan_only=0
while (( $# )); do
  case "$1" in
    --reference|--input|--output-dir|--electrons|--variants)
      (( $# >= 2 )) || die "$1 requires a value"
      case "$1" in
        --reference) reference_arg=$2 ;;
        --input) input_arg=$2 ;;
        --output-dir) output_arg=$2 ;;
        --electrons) electrons=$2 ;;
        --variants) variants=$2 ;;
      esac
      shift 2 ;;
    --plan) plan_only=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$reference_arg" && -n "$input_arg" && -n "$output_arg" && -n "$electrons" ]] || {
  usage >&2
  die '--reference, --input, --output-dir, and --electrons are required'
}
[[ -d "$reference_arg" ]] || die "reference run directory not found: $reference_arg"
[[ -f "$input_arg" ]] || die "input particle file not found: $input_arg"
REFERENCE_DIR=$(readlink -f -- "$reference_arg")
INPUT_FILE=$(readlink -f -- "$input_arg")
[[ -f "$REFERENCE_DIR/run.log" ]] || die "reference run.log not found: $REFERENCE_DIR/run.log"
[[ -f "$REFERENCE_DIR/ledger/energy-ledger.h5" ]] ||
  die "reference energy ledger not found: $REFERENCE_DIR/ledger/energy-ledger.h5"

SUITE_DIR=$(realpath -m -- "$output_arg")
[[ ! -e "$SUITE_DIR" ]] || die "refusing to reuse existing suite: $SUITE_DIR"
[[ -d "$(dirname -- "$SUITE_DIR")" ]] ||
  die "suite parent directory does not exist: $(dirname -- "$SUITE_DIR")"
[[ "$variants" =~ ^(base|wide|thick-cpml|long-guard|fine-z)(,(base|wide|thick-cpml|long-guard|fine-z))*$ ]] ||
  die 'variants must be a comma-separated list of base,wide,thick-cpml,long-guard,fine-z'
IFS=, read -r -a selected_variants <<< "$variants"
seen=,
for variant in "${selected_variants[@]}"; do
  [[ "$seen" != *",${variant},"* ]] || die "duplicate variant: $variant"
  seen+="${variant},"
done

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
APRL_SRC_DIR=${APRL_SRC_DIR:-$(cd -- "$script_dir/../.." && pwd -P)}
[[ -d "$APRL_SRC_DIR" ]] || die "APRL source directory not found: $APRL_SRC_DIR"
APRL_SRC_DIR=$(cd -- "$APRL_SRC_DIR" && pwd -P)
[[ -f "$APRL_SRC_DIR/experiments/uniform10um/prepare_held_convergence.py" ]] ||
  die "missing preparer in APRL source: $APRL_SRC_DIR"
APRL_BIN=${APRL_BIN:-$APRL_SRC_DIR/build/aprl}
LEDGER_REPORT_BIN=${LEDGER_REPORT_BIN:-$APRL_SRC_DIR/build/energy_ledger_report}
PYTHON_BIN=${PYTHON_BIN:-python3}
[[ -x "$APRL_BIN" ]] || die "APRL executable not found: $APRL_BIN"
[[ -x "$LEDGER_REPORT_BIN" ]] || die "ledger reporter not found: $LEDGER_REPORT_BIN"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_BIN"
APRL_BIN=$(readlink -f -- "$APRL_BIN")
LEDGER_REPORT_BIN=$(readlink -f -- "$LEDGER_REPORT_BIN")
PYTHON_BIN=$(readlink -f -- "$(command -v "$PYTHON_BIN")")
"$PYTHON_BIN" -c 'import math,sys; v=float(sys.argv[1]); sys.exit(0 if math.isfinite(v) and v > 0 else 1)' \
  "$electrons" >/dev/null 2>&1 || die '--electrons must be positive and finite'
ELECTRONS=$electrons
VARIANTS=$variants

PARTITION=${PARTITION:-ssmb}
NODES=${NODES:-1}
TASKS_PER_NODE=${TASKS_PER_NODE:-8}
[[ "$NODES" =~ ^[1-9][0-9]*$ && "$TASKS_PER_NODE" =~ ^[1-9][0-9]*$ ]] ||
  die 'NODES and TASKS_PER_NODE must be positive integers'
MPI_RANKS=${MPI_RANKS:-$((NODES * TASKS_PER_NODE))}
[[ "$MPI_RANKS" =~ ^[1-9][0-9]*$ ]] || die 'MPI_RANKS must be a positive integer'
(( MPI_RANKS <= NODES * TASKS_PER_NODE )) ||
  die 'MPI_RANKS cannot exceed NODES * TASKS_PER_NODE'
MPI_LAUNCHER=${MPI_LAUNCHER:-mpiexec}
[[ "$MPI_LAUNCHER" == mpiexec || "$MPI_LAUNCHER" == srun ]] ||
  die 'MPI_LAUNCHER must be mpiexec or srun'
SRUN_MPI=${SRUN_MPI:-pmix_v3}
WALLTIME=${WALLTIME:-144:0:0}
NODELIST=${NODELIST:-}
EXCLUDE_NODES=${EXCLUDE_NODES:-}
[[ -z "$NODELIST" || -z "$EXCLUDE_NODES" ]] ||
  die 'NODELIST and EXCLUDE_NODES cannot both be set'
SBATCH_ACCOUNT_REQUEST=${SBATCH_ACCOUNT:-}
SBATCH_QOS_REQUEST=${SBATCH_QOS:-}

printf 'Suite: %s\nReference: %s\nParticles: %s\nRepresented electrons: %s\nVariants: %s\n' \
  "$SUITE_DIR" "$REFERENCE_DIR" "$INPUT_FILE" "$ELECTRONS" "$VARIANTS"
printf 'Slurm: partition=%s nodes=%s tasks/node=%s MPI ranks=%s launcher=%s time=%s\n' \
  "$PARTITION" "$NODES" "$TASKS_PER_NODE" "$MPI_RANKS" "$MPI_LAUNCHER" "$WALLTIME"
echo 'RESOURCE WARNING: wide and fine-z need more mesh memory; preflight resources before submission.'
if (( NODES == 1 && MPI_RANKS == 8 )); then
  echo 'RESOURCE WARNING: 1 node / 8 ranks is a template; prior production used 2 nodes / 128 ranks.'
fi
if (( plan_only )); then
  echo 'Plan only: no directory created, no particles generated, no Slurm job submitted.'
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || die 'sbatch not found; submit from a Slurm login node'
# Snapshot Git state before making a suite inside the source tree, so the new
# suite itself cannot make an otherwise clean checkout appear dirty.
source_revision=unknown
source_status=
source_status_known=0
if source_revision=$(git -C "$APRL_SRC_DIR" rev-parse HEAD 2>/dev/null); then
  if source_status=$(git -C "$APRL_SRC_DIR" status --porcelain 2>/dev/null); then
    source_status_known=1
  fi
fi
mkdir -- "$SUITE_DIR" || die "could not create suite: $SUITE_DIR"
printf 'reference=%s\ninput=%s\nelectrons=%s\nvariants=%s\nsource=%s\n' \
  "$REFERENCE_DIR" "$INPUT_FILE" "$ELECTRONS" "$VARIANTS" "$APRL_SRC_DIR" \
  > "$SUITE_DIR/submission.txt"
printf 'nodes=%s\ntasks_per_node=%s\nmpi_ranks=%s\nmpi_launcher=%s\nsrun_mpi=%s\n' \
  "$NODES" "$TASKS_PER_NODE" "$MPI_RANKS" "$MPI_LAUNCHER" "$SRUN_MPI" \
  >> "$SUITE_DIR/submission.txt"
printf '%s\n' "$source_revision" > "$SUITE_DIR/source-revision.txt"
if (( source_status_known )); then
  if [[ -n "$source_status" ]]; then
    printf '%s\n' "$source_status" > "$SUITE_DIR/source-status.txt"
    printf 'source_worktree_dirty=1\n' >> "$SUITE_DIR/submission.txt"
    echo 'SOURCE WARNING: APRL source has uncommitted changes; inspect source-status.txt.'
  else
    printf '' > "$SUITE_DIR/source-status.txt"
    printf 'source_worktree_dirty=0\n' >> "$SUITE_DIR/submission.txt"
  fi
else
  printf 'source_worktree_dirty=unknown\n' >> "$SUITE_DIR/submission.txt"
fi
export APRL_SRC_DIR APRL_BIN LEDGER_REPORT_BIN PYTHON_BIN
export SUITE_DIR REFERENCE_DIR INPUT_FILE ELECTRONS VARIANTS
export MPI_LAUNCHER MPI_RANKS SRUN_MPI
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_QOS
sbatch_args=(
  --parsable --job-name=aprl-held-convergence --partition="$PARTITION"
  --nodes="$NODES" --ntasks="$((NODES * TASKS_PER_NODE))"
  --ntasks-per-node="$TASKS_PER_NODE" --cpus-per-task=1 --time="$WALLTIME"
  --exclusive --mem=0
  --output="$SUITE_DIR/slurm-%j.out" --error="$SUITE_DIR/slurm-%j.err"
  --export=ALL
)
[[ -z "$SBATCH_ACCOUNT_REQUEST" ]] || sbatch_args+=(--account="$SBATCH_ACCOUNT_REQUEST")
[[ -z "$SBATCH_QOS_REQUEST" ]] || sbatch_args+=(--qos="$SBATCH_QOS_REQUEST")
[[ -z "$NODELIST" ]] || sbatch_args+=(--nodelist="$NODELIST")
[[ -z "$EXCLUDE_NODES" ]] || sbatch_args+=(--exclude="$EXCLUDE_NODES")
if job_id=$(sbatch "${sbatch_args[@]}" "$script_dir/run_held_convergence.sbatch"); then
  printf 'job_id=%s\n' "$job_id" >> "$SUITE_DIR/submission.txt"
  echo "Submitted Slurm job $job_id; suite: $SUITE_DIR"
else
  printf 'submission_failed=%s\n' "$(date --iso-8601=seconds)" >> "$SUITE_DIR/submission.txt"
  die "sbatch failed; the suite directory is retained for inspection: $SUITE_DIR"
fi
