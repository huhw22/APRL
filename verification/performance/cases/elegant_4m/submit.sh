#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
APRL_BIN="${APRL_BIN:-/public/home/wfollowm/aprl_dev/aprl-runtime-lab/bin/aprl}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PARTITION="${PARTITION:-ssmb}"
NODES="${NODES:-2}"
TASKS_PER_NODE="${TASKS_PER_NODE:-128}"
THREADS_PER_RANK="${THREADS_PER_RANK:-1}"
CPU_BIND="${CPU_BIND:-cores}"
WALLTIME="${WALLTIME:-144:0:0}"
SRUN_MPI="${SRUN_MPI:-pmix_v3}"
NODELIST="${NODELIST:-}"
EXCLUDE_NODES="${EXCLUDE_NODES:-}"
SBATCH_ACCOUNT_REQUEST="${SBATCH_ACCOUNT:-}"
SBATCH_QOS_REQUEST="${SBATCH_QOS:-}"

if [[ "${1:-}" != "" && "${1:-}" != "--plan" ]] || (( $# > 1 )); then
  echo "Usage: ./submit.sh [--plan]" >&2
  exit 2
fi
[[ -f "${ROOT}/particles.h5" ]] || {
  echo "ERROR: 请把 Elegant 转换得到的 particles.h5 放在 ${ROOT}/" >&2
  exit 1
}
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
  echo "ERROR: 找不到 Python: ${PYTHON_BIN}" >&2; exit 1;
}
[[ "${NODES}" =~ ^[1-9][0-9]*$ && "${TASKS_PER_NODE}" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: NODES 和 TASKS_PER_NODE 必须是正整数" >&2; exit 2;
}
[[ "${THREADS_PER_RANK}" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: THREADS_PER_RANK 必须是正整数" >&2; exit 2;
}
[[ -z "${NODELIST}" || -z "${EXCLUDE_NODES}" ]] || {
  echo "ERROR: NODELIST 和 EXCLUDE_NODES 不能同时设置" >&2; exit 2;
}

run_dir="${ROOT}/runs/$(date -u +%Y%m%dT%H%M%SZ)-$$"
"${PYTHON_BIN}" "${ROOT}/prepare_case.py" --run-dir "${run_dir}"
printf 'nodes=%s\ntasks_per_node=%s\nthreads_per_rank=%s\ncpus_per_node_requested=%s\ncpu_bind=%s\nmpi_plugin=%s\n' \
  "${NODES}" "${TASKS_PER_NODE}" "${THREADS_PER_RANK}" \
  "$((TASKS_PER_NODE * THREADS_PER_RANK))" "${CPU_BIND}" "${SRUN_MPI}" \
  > "${run_dir}/execution-plan.txt"
if [[ "${1:-}" == "--plan" ]]; then
  echo "只生成并检查输入卡；未提交作业。"
  exit 0
fi
[[ -x "${APRL_BIN}" ]] || { echo "ERROR: APRL 不可执行: ${APRL_BIN}" >&2; exit 1; }
command -v sbatch >/dev/null 2>&1 || {
  echo "ERROR: 找不到 sbatch；请在 Slurm 登录节点运行。" >&2; exit 1;
}

export APRL_BIN CASE_RUN_DIR="${run_dir}" CASE_CONFIG="${run_dir}/run.yaml" SRUN_MPI THREADS_PER_RANK CPU_BIND PYTHON_BIN
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_QOS
args=(
  --parsable
  --job-name=aprl-elegant-4m
  --partition="${PARTITION}"
  --nodes="${NODES}"
  --ntasks="$((NODES * TASKS_PER_NODE))"
  --ntasks-per-node="${TASKS_PER_NODE}"
  --cpus-per-task="${THREADS_PER_RANK}"
  --time="${WALLTIME}"
  --exclusive --mem=0
  --output="${run_dir}/slurm-%j.out"
  --error="${run_dir}/slurm-%j.err"
  --export=ALL
)
[[ -n "${SBATCH_ACCOUNT_REQUEST}" ]] && args+=(--account="${SBATCH_ACCOUNT_REQUEST}")
[[ -n "${SBATCH_QOS_REQUEST}" ]] && args+=(--qos="${SBATCH_QOS_REQUEST}")
[[ -n "${NODELIST}" ]] && args+=(--nodelist="${NODELIST}")
[[ -n "${EXCLUDE_NODES}" ]] && args+=(--exclude="${EXCLUDE_NODES}")
job_id="$(sbatch "${args[@]}" "${ROOT}/run.sbatch")"
echo "已提交作业 ${job_id}；运行目录: ${run_dir}"
echo "日志: ${run_dir}/slurm-${job_id}.out"
