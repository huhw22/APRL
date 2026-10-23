#!/usr/bin/env bash
#SBATCH --job-name=aprl-fixed4
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=aprl-fixed4-%j.out
#SBATCH --error=aprl-fixed4-%j.err
# Submit with explicit --ntasks=R --cpus-per-task=T, where R*T=4.
# This is the small integration harness, not the production SASE/ELEGANT case.
set -Eeuo pipefail

source "${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}/tests/cluster_acceptance/common.sh"
init_run fixed-4cpu
require_pass "${APRL_SMOKE_RESULT:?Set APRL_SMOKE_RESULT to the accepted smoke result directory}" smoke
require_pass "${APRL_REQUIRED_RESULT:?Set APRL_REQUIRED_RESULT to the accepted required result directory}" required

ranks="${SLURM_NTASKS:-0}"
threads="${SLURM_CPUS_PER_TASK:-0}"
[[ "$ranks" =~ ^[1-9][0-9]*$ && "$threads" =~ ^[1-9][0-9]*$ ]] || {
  echo 'ERROR: submit with explicit --ntasks and --cpus-per-task.' >&2; exit 2;
}
[[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} == 1 ]] || {
  echo 'ERROR: fixed-4cpu requires exactly one allocated node.' >&2; exit 2;
}
(( ranks * threads == 4 )) || {
  echo 'ERROR: fixed-4cpu accepts only 4x1, 2x2, or 1x4.' >&2; exit 2;
}
printf 'candidate_layout=%sx%s\ncandidate_cpu_budget=4\nbaseline_control_layout=%sx1\n' \
  "$ranks" "$threads" "$ranks" > "$run_dir/layout.txt"

# The same-rank one-thread baseline is a correctness control. Only baseline
# 4x1 is a four-CPU performance reference across the three candidate layouts.
probe_threads=(1)
if (( threads > 1 )); then probe_threads+=("$threads"); fi
for probe_thread_count in "${probe_threads[@]}"; do
  export OMP_NUM_THREADS="$probe_thread_count"
  srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks="$ranks" \
    --ntasks-per-node="$ranks" --cpus-per-task="$probe_thread_count" \
    --cpu-bind=verbose,cores --kill-on-bad-exit=1 \
    "$workspace/build/cluster-tools/binding_probe" "$ranks" "$probe_thread_count" \
    2>&1 | tee "$run_dir/binding-r${ranks}-t${probe_thread_count}.log"
done

"$PYTHON_BIN" "$APRL_SOURCE_DIR/verification/performance/run.py" \
  --baseline "$APRL_BASELINE" --candidate "$APRL_CANDIDATE" \
  --generator "$SASE_GENERATOR" --output-dir "$run_dir/harness" \
  --ranks "$ranks" --threads "$threads" --repeats 3 --profile-interval 0 \
  --launcher srun --srun-mpi "$SRUN_MPI" \
  --launcher-arg=--exact --launcher-arg=--cpu-bind=verbose,cores \
  2>&1 | tee "$run_dir/harness.log"

"$PYTHON_BIN" - "$run_dir/harness/report.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as stream:
    report = json.load(stream)
if report.get("pass") is not True or report.get("implementation_equivalence_pass") is not True:
    raise SystemExit("ERROR: fixed-CPU implementation equivalence failed")
if report.get("input_unchanged") is not True:
    raise SystemExit("ERROR: generated input changed during testing")
print("Fixed-CPU small harness passed; production performance remains unmeasured.")
PY
printf 'fixed-4cpu passed\n' > "$run_dir/PASS"
printf 'Results: %s\n' "$run_dir"
