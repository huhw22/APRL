#!/usr/bin/env bash
#SBATCH --job-name=aprl-performance
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=aprl-performance-%j.out
#SBATCH --error=aprl-performance-%j.err
set -Eeuo pipefail
# Submit from the APRL source root, or set APRL_SOURCE_DIR explicitly.
source_dir="${APRL_SOURCE_DIR:-${SLURM_SUBMIT_DIR:-${PWD}}}"
python_bin="${PYTHON_BIN:-python3}"
: "${APRL_BASELINE:?Set APRL_BASELINE to the unchanged baseline executable}"
: "${APRL_CANDIDATE:?Set APRL_CANDIDATE to the candidate executable}"
output_dir="${PERFORMANCE_OUTPUT_DIR:-${source_dir}/work/performance/job-${SLURM_JOB_ID:-manual}}"
ranks="${RANK_COUNTS:-1,2}"
threads="${THREAD_COUNTS:-1,2}"
# No MPI/UCX/network settings are injected. Load the same runtime as both binaries.
"${python_bin}" "${source_dir}/verification/performance/run.py" \
  --baseline "${APRL_BASELINE}" --candidate "${APRL_CANDIDATE}" \
  --generator "${SASE_GENERATOR:-${source_dir}/build/sase_bunch_generator}" \
  --output-dir "${output_dir}" --ranks "${ranks}" --threads "${threads}" \
  --repeats "${REPEATS:-3}" --beamlets "${BEAMLETS:-16}" \
  --launcher srun --srun-mpi "${SRUN_MPI:-pmix_v3}" \
  --profile-interval "${APRL_PROFILE_INTERVAL:-16}" "$@"
