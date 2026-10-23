#!/usr/bin/env bash
#SBATCH -J aprl-sase4m-w1
#SBATCH -p ssmb
#SBATCH -N 2
#SBATCH -n 256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH -t 48:0:0
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -o logs/%x-%j.out
#SBATCH -e logs/%x-%j.err
#SBATCH --signal=B:TERM@120

set -Eeuo pipefail

CAMPAIGN_DIR="${CAMPAIGN_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)}"
export CAMPAIGN_DIR
set -a
# shellcheck source=/dev/null
source "${CAMPAIGN_DIR}/campaign.env"
set +a
RUN_MODE="${RUN_MODE:-full}"
JOB_TAG="${SLURM_JOB_ID:-manual}"
STATE_DIR="${CAMPAIGN_WORK_DIR}/metadata"
mkdir -p "${STATE_DIR}"

on_term() {
  printf 'STOPPED time=%s job=%s\n' "$(date --iso-8601=seconds)" "${JOB_TAG}" \
    > "${STATE_DIR}/STOPPED-${JOB_TAG}.txt"
  exit 143
}
on_error() {
  local status=$?
  printf 'FAILED status=%s line=%s time=%s\n' \
    "${status}" "${BASH_LINENO[0]:-unknown}" "$(date --iso-8601=seconds)" \
    > "${STATE_DIR}/FAILED-${JOB_TAG}.txt"
  exit "${status}"
}
trap on_term TERM INT
trap on_error ERR

if [[ "${RUN_MODE}" == "full" ]]; then
  command -v flock >/dev/null 2>&1 || {
    echo "ERROR: 找不到 flock，无法防止重复写入。" >&2
    exit 1
  }
  exec 9> "${STATE_DIR}/campaign-full.lock"
  if ! flock -n 9; then
    echo "ERROR: 另一个 full 作业正在使用 ${CAMPAIGN_WORK_DIR}。" >&2
    exit 1
  fi
fi

if [[ -n "${MODULES}" ]]; then
  if ! type module >/dev/null 2>&1; then
    [[ -r /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  fi
  type module >/dev/null 2>&1 || {
    echo "ERROR: MODULES 已设置，但计算节点没有 module 命令。" >&2
    exit 1
  }
  module purge
  read -r -a module_list <<< "${MODULES}"
  module load "${module_list[@]}"
  hash -r
fi

APRL="${APRL_PROGRAM:-}"
GENERATOR="${SASE_GENERATOR:-}"
[[ -x "${APRL}" ]] || { echo "ERROR: APRL_PROGRAM 不可执行: ${APRL}" >&2; exit 1; }
[[ -x "${GENERATOR}" ]] || { echo "ERROR: SASE_GENERATOR 不可执行: ${GENERATOR}" >&2; exit 1; }

APRL_REAL="$(readlink -f -- "${APRL}" 2>/dev/null || printf '%s\n' "${APRL}")"
runtime_candidates=(
  "${APRL_RUNTIME_ENV:-}"
  "$(dirname -- "${APRL_REAL}")/share/aprl/runtime-env.sh"
  "$(dirname -- "$(dirname -- "${APRL_REAL}")")/share/aprl/runtime-env.sh"
)
for runtime_env in "${runtime_candidates[@]}"; do
  if [[ -n "${runtime_env}" && -f "${runtime_env}" ]]; then
    # shellcheck source=/dev/null
    source "${runtime_env}"
    echo "Loaded APRL runtime environment: ${runtime_env}"
    break
  fi
done
export APRL GENERATOR

command -v srun >/dev/null 2>&1 || { echo "ERROR: 找不到 srun。" >&2; exit 1; }
if [[ ! "${SRUN_TASK_EXIT_WAIT_SECONDS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SRUN_TASK_EXIT_WAIT_SECONDS 必须是正整数，当前值: ${SRUN_TASK_EXIT_WAIT_SECONDS}" >&2
  exit 1
fi
mpi_plugin_list="$(srun --mpi=list 2>&1)"
if ! grep -q "${SRUN_MPI}" <<< "${mpi_plugin_list}"; then
  echo "ERROR: 当前 Slurm 不提供 --mpi=${SRUN_MPI}。" >&2
  printf '%s\n' "${mpi_plugin_list}" >&2
  exit 1
fi

expected_tasks=$((NODES * TASKS_PER_NODE))
[[ "${SLURM_JOB_NUM_NODES:-0}" -eq "${NODES}" ]] || {
  echo "ERROR: 实际节点数 ${SLURM_JOB_NUM_NODES:-0}，期望 ${NODES}。" >&2
  exit 1
}
[[ "${SLURM_NTASKS:-0}" -eq "${expected_tasks}" ]] || {
  echo "ERROR: 实际 rank 数 ${SLURM_NTASKS:-0}，期望 ${expected_tasks}。" >&2
  exit 1
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_btl OMPI_MCA_btl_openib_warn_no_device_params_found
ulimit -s unlimited

if [[ "${RUN_MODE}" == "full" ]]; then
  # The full run owns the exclusive campaign lock here, so this snapshot of
  # the submitted environment stays immutable for every later case.
  python3 "${CAMPAIGN_DIR}/scripts/generate_campaign.py"
fi
"${GENERATOR}" --self-test

# The preflight verdict reads the detector completion marker, so even a
# preflight-only job needs the campaign's h5py environment.  Full analysis
# checks its additional NumPy/Matplotlib imports after this gate.
PYTHON="${CAMPAIGN_DIR}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]] || ! "${PYTHON}" -c 'import h5py' >/dev/null 2>&1; then
  echo "ERROR: 计算节点不能使用 campaign .venv/h5py；请先运行 prepare_python_env.sh。" >&2
  exit 1
fi

environment_log="${STATE_DIR}/environment-${JOB_TAG}.txt"
{
  echo "start_time=$(date --iso-8601=seconds)"
  echo "job_id=${JOB_TAG}"
  echo "run_mode=${RUN_MODE}"
  echo "nodes=${SLURM_JOB_NUM_NODES}"
  echo "tasks=${SLURM_NTASKS}"
  echo "node_list=${SLURM_JOB_NODELIST}"
  echo "aprl=${APRL}"
  echo "generator=${GENERATOR}"
  echo "srun_mpi_plugin=${SRUN_MPI}"
  echo "srun_task_exit_wait_seconds=${SRUN_TASK_EXIT_WAIT_SECONDS}"
  echo "aprl_teardown_trace=${APRL_TEARDOWN_TRACE:-disabled}"
  echo "preflight_aprl_teardown_trace=all"
  echo "ompi_btl=${OMPI_MCA_btl}"
  if type module >/dev/null 2>&1; then module list 2>&1 || true; fi
  srun --mpi=list 2>&1
  ldd "${APRL}"
  ldd "${GENERATOR}"
  scontrol show job "${JOB_TAG}" 2>&1 || true
} > "${environment_log}"

host_raw="${STATE_DIR}/hosts-raw-${JOB_TAG}.txt"
host_unique="${STATE_DIR}/hosts-${JOB_TAG}.txt"
srun --nodes="${NODES}" --ntasks="${NODES}" --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --distribution=block:block --mpi=none \
  hostname > "${host_raw}"
sort -u "${host_raw}" > "${host_unique}"
actual_hosts="$(wc -l < "${host_unique}")"
[[ "${actual_hosts}" -eq "${NODES}" ]] || {
  echo "ERROR: 节点探针只看到 ${actual_hosts} 个不同主机。" >&2
  exit 1
}

ldd_log="${STATE_DIR}/ldd-per-node-${JOB_TAG}.txt"
srun --nodes="${NODES}" --ntasks="${NODES}" --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --distribution=block:block --mpi=none \
  bash -c 'echo "===== host=$(hostname) ====="; ldd "$APRL"; ldd "$GENERATOR"' \
  > "${ldd_log}" 2>&1
if grep -q 'not found' "${ldd_log}"; then
  echo "ERROR: 至少一个计算节点缺动态库；见 ${ldd_log}。" >&2
  grep -B2 -A2 'not found' "${ldd_log}" >&2 || true
  exit 1
fi
if grep -Ei 'libmpi.*(/\.conda/|/conda/)' "${ldd_log}" >/dev/null; then
  echo "ERROR: APRL 在计算节点误链接到 Conda MPI。" >&2
  exit 1
fi

mpi_runtime_log="${STATE_DIR}/mpi-runtime-per-node-${JOB_TAG}.txt"
srun --nodes="${NODES}" --ntasks="${NODES}" --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --distribution=block:block --mpi=none \
  bash -c '
    echo "===== host=$(hostname) ====="
    echo "--- MPI/PMIx command identities ---"
    for command_name in mpirun mpiexec ompi_info pmix_info; do
      command_path=$(command -v "${command_name}" 2>/dev/null || true)
      echo "${command_name}=${command_path:-MISSING}"
    done
    mpirun --version 2>&1 | head -n 4 || true
    ompi_info --version 2>&1 | head -n 8 || true
    pmix_path=$(command -v pmix_info 2>/dev/null || true)
    if [[ -n "${pmix_path}" ]]; then
      "${pmix_path}" --version 2>&1 | head -n 8 || true
      echo "--- pmix_info dynamic libraries ---"
      ldd "${pmix_path}" 2>&1 || true
    fi
    echo "--- visible PMIx libraries ---"
    /sbin/ldconfig -p 2>/dev/null | grep -Ei "libpmix|libpmi" || true
  ' > "${mpi_runtime_log}" 2>&1

# This smoke card uses cells:[...,auto], so an old APRL binary fails here
# before any large allocation.  It uses the full production rank fan-out while
# keeping the global mesh tiny, and checks two-node MPI/PMIx plus parallel
# HDF5-v2 particle input and detector output together.
smoke_dir="${CAMPAIGN_WORK_DIR}/preflight/job-${JOB_TAG}"
mkdir -p "${smoke_dir}"
"${GENERATOR}" \
  --output "${smoke_dir}/input.h5" \
  --electrons 64 \
  --length-m 1e-6 \
  --resonant-wavelength-m 1e-6 \
  --gamma 4 \
  --sigma-x-m 2e-8 \
  --sigma-y-m 2e-8 \
  --harmonics 1 \
  --beamlets-per-wavelength 8 \
  --seed 15 \
  --compression 0 \
  --overwrite
cp -- "${CAMPAIGN_DIR}/configs/mpi-smoke.yaml" "${smoke_dir}/mpi-smoke.yaml"
preflight_finalize_watchdog_raw="$(awk \
  '$1 == "finalize_timeout_seconds:" {print $2; exit}' \
  "${smoke_dir}/mpi-smoke.yaml")"
if [[ ! "${preflight_finalize_watchdog_raw}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: tiny preflight 未配置正数 finalize watchdog。" >&2
  exit 1
fi
preflight_finalize_watchdog_seconds=$((10#${preflight_finalize_watchdog_raw}))
# Keep the containment layers ordered: APRL's MPI_Finalize watchdog first, srun's
# first-rank-exit wait second, and GNU timeout last.  A user may lengthen the
# normal case wait, but cannot accidentally shorten the preflight below the
# finalize watchdog window.  The generic shutdown watchdog is deliberately
# not used here: once all ranks enter MPI_Finalize, the dedicated and shorter
# finalize timer is the one that actually contains this failure mode.
preflight_srun_wait_seconds="${SRUN_TASK_EXIT_WAIT_SECONDS}"
minimum_preflight_srun_wait=$((preflight_finalize_watchdog_seconds + 60))
if [[ "${preflight_srun_wait_seconds}" -lt "${minimum_preflight_srun_wait}" ]]; then
  preflight_srun_wait_seconds="${minimum_preflight_srun_wait}"
fi
# GNU timeout starts with srun, whereas srun --wait starts only after the
# first task exits.  Cover both serial windows and retain a separate cleanup
# margin so the two containment layers cannot expire simultaneously.
preflight_outer_timeout_seconds=$((
  preflight_finalize_watchdog_seconds + preflight_srun_wait_seconds + 120
))
preflight_log="${STATE_DIR}/mpi-hdf5-preflight-${JOB_TAG}.txt"
trace_preflight_report="${STATE_DIR}/teardown-preflight-${JOB_TAG}.txt"
h5_preflight_report="${STATE_DIR}/hdf5-preflight-${JOB_TAG}.txt"
preflight_decision_file="${STATE_DIR}/preflight-decision-${JOB_TAG}.txt"
if (
  cd -- "${smoke_dir}"
  APRL_TEARDOWN_TRACE=all timeout --kill-after=30 \
    "${preflight_outer_timeout_seconds}" \
    srun --kill-on-bad-exit=1 \
    --wait="${preflight_srun_wait_seconds}" \
    --nodes="${NODES}" --ntasks="${expected_tasks}" \
    --ntasks-per-node="${TASKS_PER_NODE}" \
    --cpus-per-task=1 --cpu-bind=cores --distribution=block:block \
    --mpi="${SRUN_MPI}" "${APRL}" "${smoke_dir}/mpi-smoke.yaml"
) > "${preflight_log}" 2>&1; then
  preflight_status=0
else
  preflight_status=$?
fi

trace_parsed=0
if python3 "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
    --expected-ranks "${expected_tasks}" "${preflight_log}" \
    > "${trace_preflight_report}" 2>&1; then
  trace_parsed=1
fi
hdf5_complete=0
if "${PYTHON}" "${CAMPAIGN_DIR}/scripts/h5_status.py" \
    "${smoke_dir}/detectors/particles.h5" /particle_plane \
    > "${h5_preflight_report}" 2>&1; then
  hdf5_complete=1
fi

preflight_outcome="rejected"
preflight_reason=""
if [[ "${trace_parsed}" -ne 1 ]]; then
  preflight_reason="teardown-trace-invalid"
elif [[ "${hdf5_complete}" -ne 1 ]]; then
  preflight_reason="particle-plane-hdf5-incomplete"
elif [[ "${preflight_status}" -eq 0 ]]; then
  if python3 "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
      --expected-ranks "${expected_tasks}" --require-clean-exit \
      "${preflight_log}" >/dev/null 2>&1; then
    preflight_outcome="clean"
  else
    preflight_reason="zero-status-with-incomplete-rank-trace"
  fi
elif python3 "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
    --expected-ranks "${expected_tasks}" --require-safe-finalize-hang \
    "${preflight_log}" >/dev/null 2>&1; then
  preflight_outcome="degraded-mpi-finalize-or-pmix"
else
  preflight_reason="nonzero-status-outside-safe-finalize-policy"
fi

if [[ "${preflight_outcome}" == "degraded-mpi-finalize-or-pmix" ]]; then
  # srun has returned, but do not start the expensive case if a watchdog-killed
  # preflight left one of its APRL processes behind in this exclusive job.
  residual_log="${STATE_DIR}/preflight-residual-processes-${JOB_TAG}.txt"
  if ! timeout --kill-after=10 60 srun --overlap \
      --nodes="${NODES}" --ntasks="${NODES}" --ntasks-per-node=1 \
      --cpus-per-task=1 --cpu-bind=cores --distribution=block:block --mpi=none \
      bash -c '
        echo "===== host=$(hostname) ====="
        command -v pgrep >/dev/null 2>&1 || exit 43
        live="$(pgrep -a -x aprl || true)"
        printf "%s\n" "${live}"
        [[ -z "${live}" ]]
      ' \
      > "${residual_log}" 2>&1; then
    preflight_outcome="rejected"
    preflight_reason="preflight-process-cleanup-not-proven"
  fi
fi

{
  echo "schema=aprl-preflight-decision-v1"
  echo "time=$(date --iso-8601=seconds)"
  echo "job=${JOB_TAG}"
  echo "outcome=${preflight_outcome}"
  echo "reason=${preflight_reason:-none}"
  echo "srun_status=${preflight_status}"
  echo "expected_ranks=${expected_tasks}"
  echo "all_rank_trace_forced=true"
  echo "trace_parsed=${trace_parsed}"
  echo "particle_plane_hdf5_complete=${hdf5_complete}"
  echo "finalize_watchdog_seconds=${preflight_finalize_watchdog_seconds}"
  echo "srun_wait_seconds=${preflight_srun_wait_seconds}"
  echo "outer_timeout_seconds=${preflight_outer_timeout_seconds}"
} > "${preflight_decision_file}"

if [[ "${preflight_outcome}" == "degraded-mpi-finalize-or-pmix" ]]; then
  {
    echo "status=PREFLIGHT_DEGRADED_FINALIZE_ONLY"
    echo "time=$(date --iso-8601=seconds)"
    echo "job=${JOB_TAG}"
    echo "srun_status=${preflight_status}"
    echo "trace_report=${trace_preflight_report}"
    echo "hdf5_report=${h5_preflight_report}"
    echo "residual_process_report=${residual_log}"
  } > "${STATE_DIR}/PREFLIGHT_DEGRADED_FINALIZE_ONLY-${JOB_TAG}.txt"
fi

if [[ "${preflight_outcome}" == "rejected" ]]; then
  echo "ERROR: 两节点 MPI/HDF5/auto-z 预检失败: ${preflight_reason}。" >&2
  cat "${trace_preflight_report}" >&2 || true
  cat "${h5_preflight_report}" >&2 || true
  tail -n 120 "${preflight_log}" >&2 || true
  exit 1
fi
cat "${trace_preflight_report}"
if [[ "${preflight_outcome}" == "degraded-mpi-finalize-or-pmix" ]]; then
  echo "WARNING: tiny preflight 在全部科学收尾和 particle-plane HDF5 完整关闭后，" >&2
  echo "WARNING: 仅 MPI_Finalize/PMIx 路径被 watchdog 终止；按受控降级策略继续。" >&2
fi
echo "Two-node ${expected_tasks}-rank APRL MPI/HDF5/auto-z preflight outcome=${preflight_outcome} on: $(paste -sd, "${host_unique}")"

if [[ "${RUN_MODE}" == "preflight" ]]; then
  printf 'PREFLIGHT_COMPLETE time=%s job=%s hosts=%s outcome=%s\n' \
    "$(date --iso-8601=seconds)" "${JOB_TAG}" \
    "$(paste -sd, "${host_unique}")" "${preflight_outcome}" \
    > "${STATE_DIR}/PREFLIGHT_COMPLETE-${JOB_TAG}.txt"
  echo "RUN_MODE=preflight：未生成大束团，未运行 4 m 单位权重工况。"
  exit 0
fi

if ! "${PYTHON}" -c 'import h5py, numpy, matplotlib' >/dev/null 2>&1; then
  echo "ERROR: 计算节点不能使用 campaign .venv；请先运行 prepare_python_env.sh。" >&2
  exit 1
fi
command -v sha256sum >/dev/null 2>&1 || {
  echo "ERROR: 找不到 sha256sum，无法安全判断续跑身份。" >&2
  exit 1
}

h5_complete() {
  "${PYTHON}" "${CAMPAIGN_DIR}/scripts/h5_status.py" "$1" "$2" >/dev/null 2>&1
}
h5_complete_or_explain() {
  local path=$1
  local group=$2
  local label=$3
  if ! "${PYTHON}" "${CAMPAIGN_DIR}/scripts/h5_status.py" \
      "${path}" "${group}"; then
    echo "ERROR: ${label} HDF5 output is missing or incomplete: ${path}:${group}" >&2
    return 1
  fi
}
hash_file() {
  sha256sum -- "$1" | awk '{print $1}'
}
case_provenance_matches() {
  local run_dir=$1
  local config=$2
  local plan=$3
  local particle_file=$4
  local case_fingerprint=$5
  local input_hash=$6
  "${PYTHON}" "${CAMPAIGN_DIR}/scripts/case_provenance_matches.py" \
    "${run_dir}" \
    --config "${config}" \
    --plan "${plan}" \
    --input "${particle_file}" \
    --fingerprint "${case_fingerprint}" \
    --input-sha256 "${input_hash}" \
    --generator-sha256 "${GENERATOR_HASH}"
}
write_ledger_report() {
  local ledger=$1
  local output=$2
  if [[ -x "${LEDGER_REPORT}" ]]; then
    local temporary="${output}.tmp-${JOB_TAG}-$$"
    if ! "${LEDGER_REPORT}" "${ledger}" > "${temporary}"; then
      rm -f -- "${temporary}"
      echo "ERROR: energy_ledger_report failed for ${ledger}." >&2
      return 1
    fi
    mv -- "${temporary}" "${output}"
  fi
}
write_complete_marker() {
  local run_dir=$1
  local case_fingerprint=$2
  local input_hash=$3
  local config_hash=$4
  local completion_source=${5:-launcher-returned}
  local recorded_srun_status=${6:-unknown}
  local cleanup_nonzero_accepted=${7:-unknown}
  local teardown_trace_verdict=${8:-unknown}
  local marker="${run_dir}/COMPLETE"
  local temporary="${marker}.tmp-${JOB_TAG}-$$"
  {
    echo "status=COMPLETE"
    echo "time=$(date --iso-8601=seconds)"
    echo "job=${JOB_TAG}"
    echo "hosts=$(paste -sd, "${host_unique}")"
    echo "fingerprint=${case_fingerprint}"
    echo "campaign_fingerprint=${campaign_fingerprint}"
    echo "aprl_sha256=${APRL_HASH}"
    echo "generator_sha256=${GENERATOR_HASH}"
    echo "input_sha256=${input_hash}"
    echo "config_sha256=${config_hash}"
    echo "completion_source=${completion_source}"
    echo "srun_status=${recorded_srun_status}"
    echo "srun_task_exit_wait_seconds=${SRUN_TASK_EXIT_WAIT_SECONDS}"
    echo "srun_cleanup_nonzero_accepted=${cleanup_nonzero_accepted}"
    echo "teardown_trace_verdict=${teardown_trace_verdict}"
  } > "${temporary}"
  mv -- "${temporary}" "${marker}"
}
write_input_provenance() {
  local path=$1
  local temporary="${path}.tmp-${JOB_TAG}-$$"
  {
    echo "schema=aprl-sase-input-provenance-v1"
    echo "generator_sha256=${GENERATOR_HASH}"
    echo "centering_script_sha256=${CENTER_SCRIPT_HASH}"
    echo "shot_noise_model=Fawley-2002-multibeamlet-multiharmonic-quiet-start"
  } > "${temporary}"
  mv -- "${temporary}" "${path}"
}
archive_path() {
  local path=$1
  local reason=${2:-incomplete}
  if [[ -e "${path}" ]]; then
    local stamp
    stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    mv -- "${path}" "${path}.${reason}-${stamp}"
  fi
}

plan="${CAMPAIGN_WORK_DIR}/metadata/campaign-plan.json"
lambda_r="$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["theory"]["resonant_wavelength_m"])' "${plan}")"
boost_gamma="$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["theory"]["boost_gamma"])' "${plan}")"
campaign_fingerprint="$("${PYTHON}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["campaign_fingerprint"])' "${plan}")"
APRL_HASH="$(hash_file "${APRL}")"
GENERATOR_HASH="$(hash_file "${GENERATOR}")"
CENTER_SCRIPT_HASH="$(hash_file "${CAMPAIGN_DIR}/scripts/center_sase_input.py")"

for beamlets in "${BEAMLETS_PER_WAVELENGTH}"; do
  particle_file="${CAMPAIGN_WORK_DIR}/inputs/sase-B$(printf '%04d' "${beamlets}").h5"
  provenance_file="${particle_file}.provenance"
  mkdir -p "$(dirname -- "${particle_file}")"
  match_args=(
    "${particle_file}"
    --seed "${RANDOM_SEED}"
    --electrons "${ELECTRONS}"
    --length-m "${BUNCH_LENGTH_M}"
    --resonant-wavelength-m "${lambda_r}"
    --gamma "${GAMMA_LAB}"
    --sigma-gamma "${SIGMA_GAMMA}"
    --sigma-x-m "${SIGMA_X_M}"
    --sigma-y-m "${SIGMA_Y_M}"
    --emittance-n-x-m "${EMITTANCE_N_X_M}"
    --emittance-n-y-m "${EMITTANCE_N_Y_M}"
    --harmonics "${SASE_HARMONICS}"
    --beamlets-per-wavelength "${beamlets}"
    --boost-gamma "${boost_gamma}"
  )
  if h5_complete "${particle_file}" /particles && \
     grep -Fxq "generator_sha256=${GENERATOR_HASH}" "${provenance_file}" 2>/dev/null && \
     grep -Fxq "centering_script_sha256=${CENTER_SCRIPT_HASH}" "${provenance_file}" 2>/dev/null && \
     "${PYTHON}" "${CAMPAIGN_DIR}/scripts/sase_input_matches.py" "${match_args[@]}"; then
    echo "[skip] matching SASE input ${particle_file}"
  else
    archive_path "${particle_file}" superseded
    archive_path "${provenance_file}" superseded
    "${GENERATOR}" \
      --output "${particle_file}" \
      --electrons "${ELECTRONS}" \
      --length-m "${BUNCH_LENGTH_M}" \
      --resonant-wavelength-m "${lambda_r}" \
      --gamma "${GAMMA_LAB}" \
      --sigma-gamma "${SIGMA_GAMMA}" \
      --sigma-x-m "${SIGMA_X_M}" \
      --sigma-y-m "${SIGMA_Y_M}" \
      --emittance-n-x-m "${EMITTANCE_N_X_M}" \
      --emittance-n-y-m "${EMITTANCE_N_Y_M}" \
      --harmonics "${SASE_HARMONICS}" \
      --beamlets-per-wavelength "${beamlets}" \
      --seed "${RANDOM_SEED}" \
      --compression "${HDF5_COMPRESSION}"
    "${PYTHON}" "${CAMPAIGN_DIR}/scripts/center_sase_input.py" \
      "${particle_file}" --boost-gamma "${boost_gamma}"
    h5_complete "${particle_file}" /particles
    "${PYTHON}" "${CAMPAIGN_DIR}/scripts/sase_input_matches.py" "${match_args[@]}"
    write_input_provenance "${provenance_file}"
  fi
done

selected_cases="${CASE_SELECTION:-${CASES}}"
if [[ -z "${selected_cases// }" ]]; then
  selected_cases="${CASES}"
fi
is_known_case() {
  local requested=$1
  local known
  for known in ${CASES}; do
    [[ "${requested}" == "${known}" ]] && return 0
  done
  return 1
}

APRL_BIN_DIR="$(cd -- "$(dirname -- "${APRL}")" && pwd -P)"
LEDGER_REPORT="${APRL_LEDGER_REPORT:-${APRL_BIN_DIR}/energy_ledger_report}"
declare -A INPUT_HASHES=()
for case_tag in ${selected_cases}; do
  is_known_case "${case_tag}" || {
    echo "ERROR: 未知 case '${case_tag}'；允许值: ${CASES}" >&2
    exit 1
  }
  config="${CAMPAIGN_WORK_DIR}/configs/${case_tag}.yaml"
  run_dir="${CAMPAIGN_WORK_DIR}/runs/${case_tag}"
  [[ -s "${config}" ]] || { echo "ERROR: 缺少 ${config}" >&2; exit 1; }
  case_beamlets="${case_tag##*-B}"
  particle_file="${CAMPAIGN_WORK_DIR}/inputs/sase-B${case_beamlets}.h5"
  [[ -s "${particle_file}" ]] || {
    echo "ERROR: 缺少 ${particle_file}" >&2
    exit 1
  }
  if [[ -z "${INPUT_HASHES[${particle_file}]:-}" ]]; then
    INPUT_HASHES[${particle_file}]="$(hash_file "${particle_file}")"
  fi
  config_hash="$(hash_file "${config}")"
  case_fingerprint="$(
    printf '%s\0' \
      'aprl-sase4m-case-v1' "${config_hash}" \
      "${INPUT_HASHES[${particle_file}]}" "${APRL_HASH}" "${GENERATOR_HASH}" \
      "${campaign_fingerprint}" |
      sha256sum | awk '{print $1}'
  )"

  archive_reason="incomplete"
  if [[ -f "${run_dir}/COMPLETE" && "${FORCE_RERUN}" != "1" ]]; then
    saved_fingerprint="$(awk -F= '$1 == "fingerprint" {print $2}' \
      "${run_dir}/COMPLETE" | tail -n 1)"
    if [[ "${saved_fingerprint}" == "${case_fingerprint}" ]] &&
       [[ -s "${run_dir}/run-manifest.yaml" ]] &&
       h5_complete "${run_dir}/ledger/energy-ledger.h5" /energy_ledger &&
       h5_complete "${run_dir}/detectors/entry-particles.h5" /particle_plane &&
       h5_complete "${run_dir}/detectors/exit-particles.h5" /particle_plane; then
      echo "[skip] complete matching case ${case_tag} fingerprint=${case_fingerprint}"
      continue
    fi
    if [[ "${saved_fingerprint}" != "${case_fingerprint}" ]]; then
      archive_reason="superseded"
      echo "[rerun] ${case_tag}: configuration/input/binary fingerprint changed"
    else
      echo "[rerun] ${case_tag}: COMPLETE exists but an output is incomplete"
    fi
  elif [[ "${FORCE_RERUN}" == "1" ]]; then
    archive_reason="forced"
  elif [[ -d "${run_dir}" ]]; then
    # APRL can close every required HDF5 writer before a small subset of MPI
    # ranks hangs during launcher/finalization cleanup.  In that case the
    # expensive result is safe to adopt only when all pre-run identities, all
    # three independent completion markers, and the strict all-rank teardown
    # trace still match the rendered case.
    if case_provenance_matches \
        "${run_dir}" "${config}" "${plan}" "${particle_file}" \
        "${case_fingerprint}" "${INPUT_HASHES[${particle_file}]}"; then
      outputs_complete=1
      h5_complete_or_explain \
        "${run_dir}/ledger/energy-ledger.h5" /energy_ledger \
        "${case_tag} energy ledger" || outputs_complete=0
      h5_complete_or_explain \
        "${run_dir}/detectors/entry-particles.h5" /particle_plane \
        "${case_tag} entry plane" || outputs_complete=0
      h5_complete_or_explain \
        "${run_dir}/detectors/exit-particles.h5" /particle_plane \
        "${case_tag} exit plane" || outputs_complete=0
      recovery_trace_verdict="invalid"
      if [[ -s "${run_dir}/aprl.log" ]]; then
        if "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
            --expected-ranks "${expected_tasks}" --require-clean-exit \
            "${run_dir}/aprl.log" >/dev/null 2>&1; then
          recovery_trace_verdict="clean-post-mpi-finalize"
        elif "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
            --expected-ranks "${expected_tasks}" \
            --require-safe-finalize-hang "${run_dir}/aprl.log" \
            >/dev/null 2>&1; then
          recovery_trace_verdict="safe-mpi-finalize-or-pmix"
        fi
      fi
      if [[ "${outputs_complete}" -eq 1 &&
            ( "${recovery_trace_verdict}" == "clean-post-mpi-finalize" ||
              "${recovery_trace_verdict}" == "safe-mpi-finalize-or-pmix" ) ]]; then
        write_ledger_report \
          "${run_dir}/ledger/energy-ledger.h5" \
          "${run_dir}/ledger-report.txt"
        write_complete_marker \
          "${run_dir}" "${case_fingerprint}" \
          "${INPUT_HASHES[${particle_file}]}" "${config_hash}" \
          adopted-complete-hdf5 unknown unknown "${recovery_trace_verdict}"
        echo "[adopt] ${case_tag}: complete HDF5 outputs and strict teardown trace recovered after interrupted launcher finalization"
        echo "[skip] adopted matching case ${case_tag} fingerprint=${case_fingerprint}"
        continue
      fi
    fi
    echo "[rerun] ${case_tag}: missing COMPLETE and recovery validation failed"
  fi
  archive_path "${run_dir}" "${archive_reason}"
  mkdir -p "${run_dir}"
  cp -- "${config}" "${run_dir}/submitted-config.yaml"
  cp -- "${plan}" "${run_dir}/submitted-plan.json"
  "${PYTHON}" "${CAMPAIGN_DIR}/scripts/export_sase_metadata.py" \
    "${particle_file}" "${run_dir}/input-noise.json" \
    --input-sha256 "${INPUT_HASHES[${particle_file}]}" \
    --generator-sha256 "${GENERATOR_HASH}"
  printf '%s\n' "${case_fingerprint}" > "${run_dir}/RUN_FINGERPRINT"
  echo "[run] ${case_tag}: ${NODES} nodes x ${TASKS_PER_NODE} ranks = ${expected_tasks}"
  echo "[run] Slurm will wait at most ${SRUN_TASK_EXIT_WAIT_SECONDS}s after the first rank exits before terminating remaining ranks"
  teardown_watcher_pid=""
  primary_step_sentinel="${run_dir}/.primary-srun-active-${JOB_TAG}"
  printf 'active job=%s time=%s\n' \
    "${JOB_TAG}" "$(date --iso-8601=seconds)" > "${primary_step_sentinel}"
  if [[ "${APRL_TEARDOWN_TRACE:-}" == "all" ]]; then
    "${CAMPAIGN_DIR}/scripts/capture_teardown_snapshot.sh" \
      "${run_dir}/aprl.log" "${run_dir}/live-process-snapshot.txt" \
      "${NODES}" 30 43200 "${primary_step_sentinel}" &
    teardown_watcher_pid=$!
  fi
  if (
    cd -- "${run_dir}"
    srun --kill-on-bad-exit=1 \
      --wait="${SRUN_TASK_EXIT_WAIT_SECONDS}" \
      --nodes="${NODES}" --ntasks="${expected_tasks}" \
      --ntasks-per-node="${TASKS_PER_NODE}" --cpus-per-task=1 \
      --cpu-bind=cores --distribution=block:block --mpi="${SRUN_MPI}" \
      "${APRL}" "${config}"
  ) 2>&1 | tee "${run_dir}/aprl.log"; then
    case_pipeline_status=("${PIPESTATUS[@]}")
  else
    # A timed-out straggler makes srun non-zero.  Capture both sides of the
    # pipeline inside the conditional so pipefail/ERR cannot discard the
    # already closed, complete simulation products before we validate them.
    case_pipeline_status=("${PIPESTATUS[@]}")
  fi
  srun_status="${case_pipeline_status[0]:-1}"
  tee_status="${case_pipeline_status[1]:-1}"
  rm -f -- "${primary_step_sentinel}"

  if [[ -n "${teardown_watcher_pid}" ]]; then
    if kill -0 "${teardown_watcher_pid}" 2>/dev/null &&
       ! grep -q 'phase=final-stop-sync-enter' \
         "${run_dir}/aprl.log" 2>/dev/null; then
      # The APRL step ended before the solve-tail trigger.  Do not leave the
      # read-only watcher polling until its long safety deadline.
      kill "${teardown_watcher_pid}" 2>/dev/null || true
    fi
    wait "${teardown_watcher_pid}" || true
  fi

  case_trace_verdict="disabled"
  if [[ "${APRL_TEARDOWN_TRACE:-}" == "all" ]]; then
    teardown_report="${run_dir}/teardown-report.txt"
    if ! "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
        --expected-ranks "${expected_tasks}" "${run_dir}/aprl.log" \
        > "${teardown_report}" 2>&1; then
      case_trace_verdict="invalid"
    elif "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
        --expected-ranks "${expected_tasks}" --require-clean-exit \
        "${run_dir}/aprl.log" >/dev/null 2>&1; then
      case_trace_verdict="clean-post-mpi-finalize"
    elif "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_teardown_trace.py" \
        --expected-ranks "${expected_tasks}" --require-safe-finalize-hang \
        "${run_dir}/aprl.log" >/dev/null 2>&1; then
      case_trace_verdict="safe-mpi-finalize-or-pmix"
    else
      case_trace_verdict="invalid"
    fi
    cat "${teardown_report}" || true
    echo "[teardown-verdict] case=${case_tag} verdict=${case_trace_verdict}"
    if [[ -s "${run_dir}/live-process-snapshot.txt" ]]; then
      echo "Live teardown snapshot: ${run_dir}/live-process-snapshot.txt"
      grep -E '^(===== host=|aprl_pid_count=|snapshot_status=|primary_step_observed_active=|mpi_finalize_complete_rank_count=|snapshot_command_status=)' \
        "${run_dir}/live-process-snapshot.txt" || true
      if grep -q '^Diagnosis: POST_MPI_FINALIZE$' "${teardown_report}"; then
        if [[ "${srun_status}" -eq 0 ]]; then
          echo "Live refinement: all ranks completed MPI_Finalize and srun returned success; no teardown hang was observed." \
            | tee -a "${teardown_report}"
        elif grep -q "^mpi_finalize_complete_rank_count=${expected_tasks}$" \
               "${run_dir}/live-process-snapshot.txt" &&
             grep -Eq '^aprl_pid_count=[1-9][0-9]*$' \
            "${run_dir}/live-process-snapshot.txt"; then
          echo "Live refinement: APRL processes remained after all mpi-finalize-complete markers; the stall is in post-finalize process/runtime exit, not Slurm-only cleanup." \
            | tee -a "${teardown_report}"
        elif grep -q "^mpi_finalize_complete_rank_count=${expected_tasks}$" \
               "${run_dir}/live-process-snapshot.txt" &&
             grep -q '^primary_step_observed_active=1$' \
               "${run_dir}/live-process-snapshot.txt" &&
             [[ "$(grep -c '^aprl_pid_count=0$' "${run_dir}/live-process-snapshot.txt" || true)" -eq "${NODES}" ]]; then
          echo "Live refinement: no APRL process remained on either node while the primary step was still being observed; the stall is in srun/slurmstepd cleanup." \
            | tee -a "${teardown_report}"
        else
          echo "Live refinement: the snapshot lacks the complete time-consistent evidence required to distinguish post-finalize process exit from Slurm cleanup, so no narrower classification is made." \
            | tee -a "${teardown_report}"
        fi
      fi
    fi
  fi

  ledger="${run_dir}/ledger/energy-ledger.h5"
  entry="${run_dir}/detectors/entry-particles.h5"
  exit_plane="${run_dir}/detectors/exit-particles.h5"
  outputs_complete=1
  h5_complete_or_explain "${ledger}" /energy_ledger \
    "${case_tag} energy ledger" || outputs_complete=0
  h5_complete_or_explain "${entry}" /particle_plane \
    "${case_tag} entry plane" || outputs_complete=0
  h5_complete_or_explain "${exit_plane}" /particle_plane \
    "${case_tag} exit plane" || outputs_complete=0
  echo "[srun-status] case=${case_tag} srun=${srun_status} tee=${tee_status} required_hdf5_complete=${outputs_complete}"
  if [[ "${tee_status}" -ne 0 ]]; then
    echo "ERROR: tee failed for ${case_tag} with status ${tee_status}; aprl.log may be incomplete." >&2
    false
  fi
  if [[ "${srun_status}" -ne 0 ]]; then
    if [[ "${outputs_complete}" -ne 1 ]]; then
      echo "ERROR: srun for ${case_tag} returned status ${srun_status} and at least one required HDF5 product is incomplete." >&2
      false
    elif [[ "${case_trace_verdict}" == "safe-mpi-finalize-or-pmix" ||
            "${case_trace_verdict}" == "clean-post-mpi-finalize" ]]; then
      echo "WARNING: srun for ${case_tag} returned status ${srun_status}, but all required HDF5 products are complete and strict all-rank trace verdict=${case_trace_verdict}; accepting the case." >&2
    else
      echo "ERROR: srun for ${case_tag} returned status ${srun_status}; complete HDF5 alone is insufficient and strict trace verdict=${case_trace_verdict}." >&2
      false
    fi
  elif [[ "${outputs_complete}" -ne 1 ]]; then
    echo "ERROR: srun for ${case_tag} returned success, but at least one required HDF5 product is incomplete." >&2
    false
  elif [[ "${APRL_TEARDOWN_TRACE:-}" == "all" &&
          "${case_trace_verdict}" != "clean-post-mpi-finalize" ]]; then
    echo "ERROR: srun for ${case_tag} returned success, but the enabled all-rank trace is not clean (verdict=${case_trace_verdict})." >&2
    false
  fi
  write_ledger_report "${ledger}" "${run_dir}/ledger-report.txt"
  cleanup_nonzero_accepted=0
  if [[ "${srun_status}" -ne 0 ]]; then
    cleanup_nonzero_accepted=1
  fi
  write_complete_marker \
    "${run_dir}" "${case_fingerprint}" \
    "${INPUT_HASHES[${particle_file}]}" "${config_hash}" \
    launcher-returned "${srun_status}" "${cleanup_nonzero_accepted}" \
    "${case_trace_verdict}"
  echo "[complete] ${case_tag}"
done

if ! "${CAMPAIGN_DIR}/analyze.sh"; then
  printf 'ANALYSIS_FAILED time=%s job=%s\n' \
    "$(date --iso-8601=seconds)" "${JOB_TAG}" \
    > "${STATE_DIR}/ANALYSIS_FAILED-${JOB_TAG}.txt"
  echo "WARNING: 仿真输出已保存，自动分析失败；可稍后手动运行 ./analyze.sh。" >&2
fi

printf 'CAMPAIGN_COMPLETE time=%s job=%s selected_cases=%s\n' \
  "$(date --iso-8601=seconds)" "${JOB_TAG}" "${selected_cases}" \
  > "${STATE_DIR}/CAMPAIGN_COMPLETE-${JOB_TAG}.txt"
echo "Campaign run complete: ${CAMPAIGN_WORK_DIR}"
