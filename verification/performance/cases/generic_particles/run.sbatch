#!/usr/bin/env bash
#SBATCH --job-name=aprl-generic-particles
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=aprl-generic-particles-%j.out
#SBATCH --error=aprl-generic-particles-%j.err
set -Eeuo pipefail
source "${APRL_WORKSPACE:-/public/home/wfollowm/APRL_Workspace}/tests/cluster_acceptance/common.sh"
trap 'generic_exit_code=$?; if [[ -n "${run_dir:-}" ]]; then printf "%s\n" "$generic_exit_code" > "$run_dir/EXIT_CODE"; fi' EXIT
init_run generic-particles
require_pass "${APRL_SMOKE_RESULT:?Set accepted smoke result directory}" smoke
require_pass "${APRL_REQUIRED_RESULT:?Set accepted required result directory}" required
prepared="$workspace/tests/generic_particles"
cp "$prepared/prepared-manifest.sha256" "$run_dir/generic-prepared-manifest.sha256"
sha256sum --check --quiet "$run_dir/generic-prepared-manifest.sha256"
[[ "${GENERIC_CASE:?Set the absolute frozen generated-particle case path}" == /* ]] || fail 'GENERIC_CASE must be absolute.'
arguments=(--workspace "$workspace" --case "$GENERIC_CASE" --output-dir "$run_dir/scientific"
           --mode "${GENERIC_MODE:-calibration}" --srun-mpi "$SRUN_MPI")
if [[ -n "${GENERIC_REPEATS:-}" ]]; then arguments+=(--repeats "$GENERIC_REPEATS"); fi
if [[ -n "${GENERIC_LAYOUTS:-}" ]]; then arguments+=(--layouts "$GENERIC_LAYOUTS"); fi
"$PYTHON_BIN" "$prepared/run.py" "${arguments[@]}"
printf 'Generated-particle binding, completed physical endpoint, and same-rank equivalence passed\n' > "$run_dir/PASS"
printf 'Results: %s\n' "$run_dir"
