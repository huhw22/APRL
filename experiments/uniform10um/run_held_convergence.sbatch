#!/usr/bin/env bash
# Slurm entry point for submit_held_convergence.sh. Paths arrive as absolute exports.
#SBATCH -J aprl-held-convergence
#SBATCH -p ssmb
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 144:0:0
#SBATCH --exclusive
#SBATCH --mem=0

set -Eeuo pipefail

die() { echo "ERROR: $*" >&2; exit 2; }
for required in APRL_SRC_DIR APRL_BIN LEDGER_REPORT_BIN PYTHON_BIN SUITE_DIR \
                REFERENCE_DIR INPUT_FILE ELECTRONS VARIANTS MPI_LAUNCHER MPI_RANKS; do
  [[ -n "${!required:-}" ]] || die "missing job argument: $required"
done
[[ -d "$SUITE_DIR" && -d "$REFERENCE_DIR" && -f "$INPUT_FILE" ]] ||
  die 'suite, reference, or input particle file is unavailable on the compute node'
[[ -f "$APRL_SRC_DIR/experiments/uniform10um/prepare_held_convergence.py" ]] ||
  die 'prepare_held_convergence.py is unavailable on the compute node'

if [[ -n "${MODULES:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
  fi
  type module >/dev/null 2>&1 || die 'module command is unavailable'
  module purge
  read -r -a selected_modules <<< "$MODULES"
  module load "${selected_modules[@]}"
  hash -r
fi
[[ -x "$APRL_BIN" && -x "$LEDGER_REPORT_BIN" ]] ||
  die 'APRL or energy_ledger_report executable is unavailable on the compute node'
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python is unavailable: $PYTHON_BIN"
[[ "$MPI_RANKS" =~ ^[1-9][0-9]*$ ]] || die 'MPI_RANKS must be a positive integer'
(( MPI_RANKS <= ${SLURM_NTASKS:-MPI_RANKS} )) ||
  die 'MPI_RANKS exceeds allocated Slurm tasks'

case "$MPI_LAUNCHER" in
  mpiexec)
    command -v mpiexec >/dev/null 2>&1 || die 'mpiexec is unavailable'
    launch=(mpiexec -n "$MPI_RANKS") ;;
  srun)
    command -v srun >/dev/null 2>&1 || die 'srun is unavailable'
    srun --mpi=list 2>&1 | grep -Fq -- "${SRUN_MPI:-pmix_v3}" ||
      die "Slurm does not advertise MPI plugin ${SRUN_MPI:-pmix_v3}"
    launch=(srun --mpi="${SRUN_MPI:-pmix_v3}" --ntasks="$MPI_RANKS" --kill-on-bad-exit=1) ;;
  *) die "invalid MPI_LAUNCHER: $MPI_LAUNCHER" ;;
esac

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
cd -- "$SUITE_DIR"
printf 'preparing start=%s host=%s job=%s\n' \
  "$(date --iso-8601=seconds)" "$(hostname)" "${SLURM_JOB_ID:-manual}"
printf 'preparing\n' > "$SUITE_DIR/job-status.txt"
if ! "$PYTHON_BIN" "$APRL_SRC_DIR/experiments/uniform10um/prepare_held_convergence.py" \
  --input "$INPUT_FILE" --output-dir "$SUITE_DIR" \
  --electrons "$ELECTRONS" --variants "$VARIANTS" \
  > "$SUITE_DIR/prepare.log" 2>&1; then
  printf 'prepare-failed\n' > "$SUITE_DIR/job-status.txt"
  die "card preparation failed; see $SUITE_DIR/prepare.log"
fi

printf 'variant\texit_code\tnormal_stop\tledger_complete\tdetectors_present\tstatus\trun_dir\n' \
  > "$SUITE_DIR/status.tsv"
IFS=, read -r -a selected_variants <<< "$VARIANTS"
failures=0
for variant in "${selected_variants[@]}"; do
  card="$SUITE_DIR/$variant/input.yaml"
  run_dir="$SUITE_DIR/$variant/runs/g929-1m-held-field"
  ledger="$run_dir/ledger/energy-ledger.h5"
  mkdir -p -- "$run_dir"
  printf 'running %s\n' "$variant" > "$SUITE_DIR/job-status.txt"
  printf 'starting=%s at=%s ranks=%s card=%s log=%s\n' \
    "$variant" "$(date --iso-8601=seconds)" "$MPI_RANKS" "$card" "$run_dir/run.log"
  exit_code=0
  if [[ ! -f "$card" ]]; then
    printf 'ERROR: missing card: %s\n' "$card" > "$run_dir/run.log"
    exit_code=2
  elif "${launch[@]}" "$APRL_BIN" "$card" > "$run_dir/run.log" 2>&1; then
    :
  else
    exit_code=$?
  fi
  normal_stop=0
  if grep -Fq 'Configured stop reached: all' "$run_dir/run.log"; then
    normal_stop=1
  fi
  ledger_complete=0
  if [[ -f "$ledger" ]] && \
     "$LEDGER_REPORT_BIN" "$ledger" > "$run_dir/ledger-report.txt" 2> "$run_dir/ledger-report.err" && \
     grep -Eq '(^|[[:space:]])complete=1([[:space:]]|$)' "$run_dir/ledger-report.txt"; then
    ledger_complete=1
  fi
  detectors_present=0
  if [[ -s "$run_dir/detectors/near-field.h5" &&
        -s "$run_dir/detectors/after-1m.h5" ]]; then
    detectors_present=1
  fi
  status=complete
  if (( exit_code != 0 || normal_stop != 1 || ledger_complete != 1 ||
        detectors_present != 1 )); then
    status=failed-or-incomplete
    (( failures += 1 ))
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$variant" "$exit_code" "$normal_stop" "$ledger_complete" \
    "$detectors_present" "$status" "$run_dir" \
    >> "$SUITE_DIR/status.tsv"
  printf 'finished=%s at=%s exit=%s normal_stop=%s ledger_complete=%s detectors_present=%s\n' \
    "$variant" "$(date --iso-8601=seconds)" "$exit_code" "$normal_stop" \
    "$ledger_complete" "$detectors_present"
done
if (( failures )); then
  printf 'failed-or-incomplete variants=%s\n' "$failures" > "$SUITE_DIR/job-status.txt"
  echo "Finished with $failures failed/incomplete variant(s); see $SUITE_DIR/status.tsv" >&2
  exit 1
fi
printf 'complete\n' > "$SUITE_DIR/job-status.txt"
echo "All held convergence variants completed: $SUITE_DIR/status.tsv"
