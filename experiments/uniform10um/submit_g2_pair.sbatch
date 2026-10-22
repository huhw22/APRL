#!/usr/bin/env bash
# Example only: submit from the extracted APRL source root after generating w1.
# `sbatch experiments/uniform10um/submit_g2_pair.sbatch`
# Slurm copies this script into its spool, so locate the project through
# SLURM_SUBMIT_DIR rather than through BASH_SOURCE or the current spool path.
#SBATCH -J APRL_10UM_G2
#SBATCH -p ssmb
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -t 144:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

set -euo pipefail

repo_root="${APRL_REPO_ROOT:-${SLURM_SUBMIT_DIR:?Submit from the APRL source root}}"
cd -- "${repo_root}"
if [[ ! -f CMakeLists.txt || ! -x build/aprl ]]; then
  echo "APRL root/build not found at ${repo_root}" >&2
  exit 2
fi
input_file="${repo_root}/work/uniform10um/particles-10um-50A-w1.h5"
if [[ ! -f "${input_file}" ]]; then
  echo "Generate and verify the w1 input before submitting: ${input_file}" >&2
  exit 2
fi
if ! command -v mpiexec >/dev/null 2>&1; then
  echo "mpiexec not found; load the MPI environment used to build APRL" >&2
  exit 2
fi
ranks="${APRL_MPI_RANKS:-${SLURM_NTASKS:-8}}"
if [[ ! "${ranks}" =~ ^[1-9][0-9]*$ ]] ||
   (( ranks > ${SLURM_NTASKS:-8} )); then
  echo "APRL_MPI_RANKS must be positive and <= allocated SLURM_NTASKS" >&2
  exit 2
fi

case_root="${repo_root}/work/uniform10um/g2-short-w1-job-${SLURM_JOB_ID:?}"
if [[ -e "${case_root}" ]]; then
  echo "Refusing to reuse existing output directory: ${case_root}" >&2
  exit 2
fi
python3 experiments/uniform10um/prepare.py \
  --input "${input_file}" --output-dir "${case_root}" \
  --profile g2-short --electrons 10409720 \
  --cases held-field,ballistic-field

for case_name in held-field ballistic-field; do
  run_dir="${case_root}/runs/g2-short-${case_name}"
  mkdir -p -- "${run_dir}"
  echo "Starting ${case_name} with ${ranks} MPI ranks; log=${run_dir}/run.log"
  mpiexec -n "${ranks}" "${repo_root}/build/aprl" \
    "${case_root}/cards/g2-short-${case_name}.yaml" \
    > "${run_dir}/run.log" 2>&1
done

python3 experiments/uniform10um/analyze.py \
  "${case_root}/runs/g2-short-held-field" \
  "${case_root}/runs/g2-short-ballistic-field" \
  --json-out "${case_root}/analysis.json" \
  --csv-out "${case_root}/energy-traces.csv"
echo "Completed paired run: ${case_root}"
