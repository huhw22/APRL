#!/usr/bin/env bash
#SBATCH --job-name=aprl-smoke
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --exclusive
# Reserve the whole 128-CPU node; exact job steps retain the small test layouts.
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=aprl-smoke-%j.out
#SBATCH --error=aprl-smoke-%j.err
set -Eeuo pipefail

source "${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}/tests/cluster_acceptance/common.sh"
init_run smoke
[[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} == 1 ]] || {
  echo 'ERROR: smoke requires exactly one allocated node.' >&2; exit 2;
}
(( ${SLURM_NTASKS:-0} >= 2 && ${SLURM_CPUS_PER_TASK:-0} >= 2 )) || {
  echo 'ERROR: smoke requires at least 2 tasks and 2 CPUs per task.' >&2; exit 2;
}

for ranks in 1 2; do
  for threads in 1 2; do
    export OMP_NUM_THREADS="$threads"
    srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks="$ranks" \
      --ntasks-per-node="$ranks" --cpus-per-task="$threads" \
      --cpu-bind=verbose,cores --kill-on-bad-exit=1 \
      "$workspace/build/cluster-tools/binding_probe" "$ranks" "$threads" \
      2>&1 | tee "$run_dir/binding-r${ranks}-t${threads}.log"
  done
done

"$PYTHON_BIN" "$APRL_SOURCE_DIR/verification/performance/run.py" \
  --baseline "$APRL_BASELINE" --candidate "$APRL_CANDIDATE" \
  --generator "$SASE_GENERATOR" --output-dir "$run_dir/harness" \
  --ranks 1,2 --threads 1,2 --repeats 1 --profile-interval 16 \
  --launcher srun --srun-mpi "$SRUN_MPI" \
  --launcher-arg=--exact --launcher-arg=--cpu-bind=verbose,cores \
  2>&1 | tee "$run_dir/harness.log"

"$PYTHON_BIN" - "$run_dir/harness/report.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as stream:
    report = json.load(stream)
if report.get("pass") is not True or report.get("implementation_equivalence_pass") is not True:
    raise SystemExit("ERROR: smoke implementation equivalence failed")
if report.get("input_unchanged") is not True:
    raise SystemExit("ERROR: generated smoke input changed during testing")
print("Smoke passed; inherited cross-rank consistency is reported separately:",
      report.get("cross_rank_consistency_pass"))
PY
printf 'smoke passed\n' > "$run_dir/PASS"
printf 'Results: %s\n' "$run_dir"
