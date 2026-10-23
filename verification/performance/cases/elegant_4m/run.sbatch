#!/usr/bin/env bash
#SBATCH -J aprl-elegant-4m
#SBATCH -p ssmb
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 144:0:0
#SBATCH --exclusive
#SBATCH --mem=0

set -Eeuo pipefail

[[ -x "${APRL_BIN:-}" && -f "${CASE_CONFIG:-}" && -d "${CASE_RUN_DIR:-}" ]] || {
  echo "ERROR: APRL 或运行配置不可用" >&2; exit 1;
}
if [[ -n "${MODULES:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
  fi
  type module >/dev/null 2>&1 || { echo "ERROR: 找不到 module 命令" >&2; exit 1; }
  module purge
  read -r -a selected_modules <<< "${MODULES}"
  module load "${selected_modules[@]}"
  hash -r
fi
APRL_BIN="$(readlink -f -- "${APRL_BIN}")"
prefix="$(dirname -- "$(dirname -- "${APRL_BIN}")")"
runtime_env="${prefix}/share/aprl/runtime-env.sh"
if [[ -f "${runtime_env}" ]]; then
  # shellcheck source=/dev/null
  source "${runtime_env}"
fi
THREADS_PER_RANK="${THREADS_PER_RANK:-1}"
[[ "${THREADS_PER_RANK}" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: THREADS_PER_RANK 必须是正整数" >&2; exit 2;
}
if (( THREADS_PER_RANK > ${SLURM_CPUS_PER_TASK:-1} )); then
  echo "ERROR: THREADS_PER_RANK 超过 Slurm 分配的每任务 CPU 数" >&2; exit 2
fi
export OMP_NUM_THREADS="${THREADS_PER_RANK}"
export OMP_PLACES="${OMP_PLACES:-cores}" OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CPU_BIND="${CPU_BIND:-cores}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,vader,tcp}"
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
ulimit -s unlimited

"${PYTHON_BIN}" - "${CASE_RUN_DIR}/input-plan.json" <<'PY'
import json
import os
import sys
plan = json.load(open(sys.argv[1], encoding="utf-8"))
stat = os.stat(plan["input_file"])
if (stat.st_size, stat.st_mtime_ns) != (plan["input_size_bytes"], plan["input_mtime_ns"]):
    raise SystemExit("ERROR: particles.h5 在提交后发生变化；拒绝使用旧配置启动")
PY

"${PYTHON_BIN}" - "${CASE_RUN_DIR}/execution-environment.json" <<'PYENV'
import json, os, platform, sys
keys = ["APRL_PROFILE_INTERVAL", "OMP_NUM_THREADS", "OMP_PLACES", "OMP_PROC_BIND", "OMP_DYNAMIC",
        "SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK",
        "SLURM_TASKS_PER_NODE", "SRUN_MPI", "CPU_BIND", "OMPI_MCA_btl",
        "OMPI_MCA_pml", "OMPI_MCA_btl_tcp_if_include", "OMPI_MCA_btl_tcp_if_exclude",
        "UCX_TLS", "UCX_NET_DEVICES", "FI_PROVIDER", "HDF5_USE_FILE_LOCKING"]
with open(sys.argv[1], "w", encoding="utf-8") as output:
    json.dump({"host": platform.node(), "platform": platform.platform(),
               "environment": {key: os.environ.get(key) for key in keys},
               "launcher_cpu_affinity": sorted(os.sched_getaffinity(0))}, output, indent=2)
    output.write("\n")
PYENV
command -v lscpu >/dev/null 2>&1 && lscpu > "${CASE_RUN_DIR}/cpu-topology.txt"
echo "start=$(date --iso-8601=seconds) host=$(hostname) job=${SLURM_JOB_ID:-manual}"
echo "APRL=${APRL_BIN} input=${CASE_CONFIG}"
echo "MPI plugin=${SRUN_MPI} tasks=${SLURM_NTASKS:-unknown}"
if ! srun --mpi=list 2>&1 | grep -q "${SRUN_MPI}"; then
  echo "ERROR: Slurm 不支持 MPI 插件 ${SRUN_MPI}" >&2
  exit 1
fi
if ! srun --mpi="${SRUN_MPI}" --cpus-per-task="${THREADS_PER_RANK}" \
  --cpu-bind="${CPU_BIND}" --kill-on-bad-exit=1 --wait=300 \
  bash -c 'affinity=$(sed -n "s/^Cpus_allowed_list:[[:space:]]*//p" /proc/self/status); printf "[launch] rank=%s host=%s cpus=%s threads=%s\n" "${SLURM_PROCID:-unknown}" "$(hostname)" "$affinity" "$OMP_NUM_THREADS"; exec "$@"' \
  bash "${APRL_BIN}" "${CASE_CONFIG}" > "${CASE_RUN_DIR}/aprl.log" 2>&1; then
  echo "ERROR: APRL 运行失败；见 ${CASE_RUN_DIR}/aprl.log" >&2
  exit 1
fi
"${PYTHON_BIN}" - "${CASE_RUN_DIR}" <<'PYCHECK'
import json, sys
from pathlib import Path
import h5py
root = Path(sys.argv[1])
log = (root / "aprl.log").read_text(errors="replace")
if "Configured stop reached:" not in log:
    raise SystemExit("ERROR: APRL 未报告到达物理终点，不能只按退出码判定成功")
expected = ["detectors/entry-particles.h5", "detectors/after-1m-particles.h5",
            "detectors/after-2m-particles.h5", "detectors/after-3m-particles.h5",
            "detectors/exit-field-ballistic-reference.h5", "detectors/exit-field.h5",
            "ledger/energy-ledger.h5"]
products = [root / name for name in expected]
missing = [str(path) for path in products if not path.is_file()]
if missing:
    raise SystemExit("ERROR: 缺少诊断 HDF5 输出: " + ", ".join(missing))
checks = []
for path in sorted(products):
    with h5py.File(path, "r") as file:
        found = False
        for name in ("particle_plane", "field_plane", "ballistic_reference", "energy_ledger"):
            if name not in file:
                continue
            found = True
            group = file[name]
            if int(group["complete"][()]) != 1:
                raise SystemExit("ERROR: 未完成的诊断输出 " + str(path))
            dataset = "fields" if name == "field_plane" else "records"
            marker = "committed_samples" if name == "field_plane" else "committed_records"
            committed = int(group[marker][()])
            if committed != len(group[dataset]) or committed == 0 or (name == "field_plane" and len(group["time_s"]) != committed):
                raise SystemExit("ERROR: 诊断记录为空或未完全提交 " + str(path))
            checks.append({"file": str(path), "group": name, "records": committed, "complete": True})
        if not found:
            raise SystemExit("ERROR: 无法识别诊断输出结构 " + str(path))
(root / "completion-check.json").write_text(json.dumps(checks, indent=2) + "\n")
PYCHECK
echo "finished=$(date --iso-8601=seconds)"
