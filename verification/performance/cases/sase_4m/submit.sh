#!/usr/bin/env bash
set -Eeuo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export CAMPAIGN_DIR
set -a
# shellcheck source=/dev/null
source "${CAMPAIGN_DIR}/campaign.env"
set +a

mode="full"
case "${1:-}" in
  --plan) mode="plan" ;;
  --preflight) mode="preflight" ;;
  "") ;;
  *) echo "Usage: ./submit.sh [--plan|--preflight]" >&2; exit 2 ;;
esac
if [[ ! "${SRUN_TASK_EXIT_WAIT_SECONDS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SRUN_TASK_EXIT_WAIT_SECONDS 必须是正整数，当前值: ${SRUN_TASK_EXIT_WAIT_SECONDS}" >&2
  exit 2
fi

resolve_executable() {
  local name=$1
  local explicit=${2:-}
  local candidate
  for candidate in \
    "${explicit}" \
    "${APRL_PROGRAM_DIR%/}/${name}" \
    "${APRL_PROGRAM_DIR%/}/bin/${name}" \
    "${APRL_PROGRAM_DIR%/}/build/${name}"; do
    if [[ -n "${candidate}" && -x "${candidate}" ]]; then
      readlink -f -- "${candidate}" 2>/dev/null || printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

APRL="$(resolve_executable aprl "${APRL_PROGRAM:-}")" || {
  echo "ERROR: 找不到预编译 APRL。" >&2
  echo "请在 campaign.env 设置 APRL_PROGRAM 或 APRL_PROGRAM_DIR。" >&2
  exit 1
}
generator_request="${SASE_GENERATOR:-}"
if [[ -z "${generator_request}" && -x "$(dirname -- "${APRL}")/sase_bunch_generator" ]]; then
  generator_request="$(dirname -- "${APRL}")/sase_bunch_generator"
fi
GENERATOR="$(resolve_executable sase_bunch_generator "${generator_request}")" || {
  echo "ERROR: 找不到同版本 sase_bunch_generator。" >&2
  echo "请设置 SASE_GENERATOR，或把它与 APRL 放在同一目录。" >&2
  exit 1
}

# A locally bundled runtime may be needed even for the login-node self-test.
APRL_REAL="$(readlink -f -- "${APRL}" 2>/dev/null || printf '%s\n' "${APRL}")"
for runtime_env in \
  "${APRL_RUNTIME_ENV:-}" \
  "$(dirname -- "${APRL_REAL}")/share/aprl/runtime-env.sh" \
  "$(dirname -- "$(dirname -- "${APRL_REAL}")")/share/aprl/runtime-env.sh"; do
  if [[ -n "${runtime_env}" && -f "${runtime_env}" ]]; then
    # shellcheck source=/dev/null
    source "${runtime_env}"
    break
  fi
done

command -v python3 >/dev/null 2>&1 || {
  echo "ERROR: 登录节点找不到 python3。" >&2
  exit 1
}
command -v flock >/dev/null 2>&1 || {
  echo "ERROR: 找不到 flock，无法保护共享配置。" >&2
  exit 1
}
mkdir -p "${CAMPAIGN_WORK_DIR}/metadata"
exec 8> "${CAMPAIGN_WORK_DIR}/metadata/campaign-full.lock"
if ! flock -n 8; then
  echo "ERROR: 正式 campaign 正在运行；当前操作不会改写其配置。" >&2
  exit 1
fi
python3 "${CAMPAIGN_DIR}/scripts/generate_campaign.py"
"${GENERATOR}" --self-test
flock -u 8
exec 8>&-

if [[ "${mode}" == "plan" ]]; then
  echo "PLAN_ONLY：单位权重、1.35 nm 单工况输入卡和静态约束已生成；未生成大束团，未提交作业。"
  echo "Plan: ${CAMPAIGN_WORK_DIR}/metadata/campaign-plan.md"
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || {
  echo "ERROR: 找不到 sbatch；请在 Slurm 登录节点提交。" >&2
  exit 1
}
[[ "${NODES}" -eq 2 ]] || {
  echo "ERROR: 这个对照脚本固定使用 2 nodes。" >&2
  exit 1
}
[[ -z "${NODELIST}" || -z "${EXCLUDE_NODES}" ]] || {
  echo "ERROR: NODELIST 与 EXCLUDE_NODES 不能同时设置。" >&2
  exit 1
}
if [[ "${mode}" != "plan" ]]; then
  PYTHON="${CAMPAIGN_DIR}/.venv/bin/python"
  if [[ ! -x "${PYTHON}" ]] || ! "${PYTHON}" -c 'import h5py' >/dev/null 2>&1; then
    echo "ERROR: 请先运行 ./prepare_python_env.sh（preflight 需要 h5py 核对完整标记）。" >&2
    exit 1
  fi
  if [[ "${mode}" == "full" ]] && \
     ! "${PYTHON}" -c 'import h5py, numpy, matplotlib' >/dev/null 2>&1; then
    echo "ERROR: 请先运行 ./prepare_python_env.sh。" >&2
    exit 1
  fi
fi

mkdir -p "${CAMPAIGN_DIR}/logs" "${CAMPAIGN_WORK_DIR}"
requested_partition="${SBATCH_PARTITION}"
requested_account="${SBATCH_ACCOUNT}"
requested_qos="${SBATCH_QOS}"
requested_nodelist="${NODELIST}"
requested_exclude="${EXCLUDE_NODES}"
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_QOS

if [[ "${mode}" == "preflight" ]]; then
  # Exercise the same cross-node rank fan-out as the production cases while
  # using a tiny mesh, so OpenMPI/PMIx failures appear before the expensive run.
  submit_tasks=$((NODES * TASKS_PER_NODE))
  submit_tasks_per_node="${TASKS_PER_NODE}"
  job_name="aprl-sase4m-w1-check"
else
  submit_tasks=$((NODES * TASKS_PER_NODE))
  submit_tasks_per_node="${TASKS_PER_NODE}"
  job_name="aprl-sase4m-w1"
fi

# Slurm gives variables already present in the submission environment
# precedence when --export contains ALL.  Publish the resolved executables and
# selected mode here, then pass the environment as a whole; otherwise an empty
# APRL_PROGRAM inherited from campaign.env can shadow the resolved path.
export APRL_PROGRAM="${APRL}"
export SASE_GENERATOR="${GENERATOR}"
export RUN_MODE="${mode}"
export CASE_SELECTION

args=(
  --parsable
  --job-name="${job_name}"
  --partition="${requested_partition}"
  --nodes="${NODES}"
  --ntasks="${submit_tasks}"
  --ntasks-per-node="${submit_tasks_per_node}"
  --cpus-per-task=1
  --time="${WALLTIME}"
  --exclusive
  --mem=0
  --output="${CAMPAIGN_DIR}/logs/%x-%j.out"
  --error="${CAMPAIGN_DIR}/logs/%x-%j.err"
  --export=ALL
)
[[ -n "${requested_account}" ]] && args+=(--account="${requested_account}")
[[ -n "${requested_qos}" ]] && args+=(--qos="${requested_qos}")
[[ -n "${requested_nodelist}" ]] && args+=(--nodelist="${requested_nodelist}")
[[ -n "${requested_exclude}" ]] && args+=(--exclude="${requested_exclude}")

echo "Slurm request: mode=${mode} partition=${requested_partition} nodes=${NODES} tasks-per-node=${submit_tasks_per_node} total-tasks=${submit_tasks} time=${WALLTIME}"
echo "Slurm first-rank-exit wait: ${SRUN_TASK_EXIT_WAIT_SECONDS}s"
echo "APRL: ${APRL}"
echo "SASE generator: ${GENERATOR}"
if [[ -n "${CASE_SELECTION}" ]]; then
  echo "Selected cases: ${CASE_SELECTION}"
else
  echo "Cases: ${CASES}"
fi
if [[ -n "${requested_nodelist}" ]]; then
  echo "Requested nodes: ${requested_nodelist}"
elif [[ -n "${requested_exclude}" ]]; then
  echo "Excluded nodes: ${requested_exclude}"
else
  echo "Nodes: 由 Slurm 选择。"
fi

job_id="$(sbatch "${args[@]}" "${CAMPAIGN_DIR}/run_campaign.sbatch")"
echo "Submitted job ${job_id}"
echo "Queue: squeue -j ${job_id}"
echo "Output: ${CAMPAIGN_DIR}/logs/${job_name}-${job_id}.out"
echo "Errors: ${CAMPAIGN_DIR}/logs/${job_name}-${job_id}.err"
echo "Accounting: sacct -j ${job_id} --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ExitCode"
