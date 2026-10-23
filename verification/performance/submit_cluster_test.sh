#!/usr/bin/env bash
# Build both revisions with one toolchain, then submit the small comparison.
set -Eeuo pipefail

usage() {
  cat <<'EOF'
用法: bash verification/performance/submit_cluster_test.sh [--prepare-only] 旧源码目录 依赖前缀

例如:
  bash verification/performance/submit_cluster_test.sh /public/home/wfollowm/aprl_dev /public/home/wfollowm/aprl_dev/aprl-runtime-lab

沿用当前 MPI/编译器环境，不自动加载模块或更改网络设置。
可选环境变量: JOBS=4 PARTITION=ssmb SRUN_MPI=pmix_v3 REPEATS=1 PROFILE=16
PYTHON_BIN 可选择已安装 numpy 和 h5py 的 Python。--prepare-only 只检查和编译，不提交。
APRL_BASELINE_BUILD_DIR / APRL_CANDIDATE_BUILD_DIR 可外置两个构建目录。
PERFORMANCE_SUBMIT_DIR 可外置提交工作目录和 Slurm 日志；PERFORMANCE_OUTPUT_DIR 可外置结果。
EOF
}

fail() { printf '错误: %s\n' "$*" >&2; exit 1; }
prepare_only=0
paths=()
for argument in "$@"; do
  case "$argument" in
    --prepare-only) prepare_only=1 ;;
    -h|--help) usage; exit 0 ;;
    --*) fail "未知选项: ${argument}" ;;
    *) paths+=("$argument") ;;
  esac
done
if [[ ${#paths[@]} != 2 ]]; then usage >&2; exit 2; fi

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
[[ -d "${paths[0]}" ]] || fail "旧源码目录不存在: ${paths[0]}"
[[ -d "${paths[1]}" ]] || fail "依赖前缀不存在: ${paths[1]}"
old_source="$(cd -- "${paths[0]}" && pwd -P)"
dependency_prefix="$(cd -- "${paths[1]}" && pwd -P)"
[[ "$old_source" != "$source_dir" ]] || fail '旧源码与新版源码指向同一目录，无法进行版本对照。'
for revision_source in "$old_source" "$source_dir"; do
  [[ -f "$revision_source/CMakeLists.txt" ]] || fail "不是 APRL 源码目录，缺少 CMakeLists.txt: $revision_source"
done
[[ -f "$source_dir/verification/performance/run.sbatch" ]] || fail '新版源码缺少 verification/performance/run.sbatch。'
[[ -f "$source_dir/verification/performance/run.py" ]] || fail '新版源码缺少 verification/performance/run.py。'

printf '[1/5] 核对源码和运行环境\n新版源码: %s\n旧版源码: %s\n依赖前缀: %s\n' "$source_dir" "$old_source" "$dependency_prefix"
printf '沿用当前环境；不会自动加载模块或选择 MPI 网络。\n'
runtime_env="$dependency_prefix/share/aprl/runtime-env.sh"
if [[ -f "$runtime_env" ]]; then
  printf '读取运行环境: %s\n' "$runtime_env"
  # shellcheck source=/dev/null
  source "$runtime_env"
fi

# Preserve the executable basename: MPI compiler wrappers can depend on argv[0].
absolute_command() {
  local selected
  selected="$(command -v -- "$1")" || fail "找不到命令: $1；请先准备现有编译器/MPI/Python 环境。"
  [[ -f "$selected" && -x "$selected" ]] || fail "命令不是可执行文件: $selected"
  printf '%s/%s\n' "$(cd -- "$(dirname -- "$selected")" && pwd -P)" "$(basename -- "$selected")"
}
cmake_bin="$(absolute_command cmake)"
mpi_cxx="$(absolute_command mpic++)"
python_bin="$(absolute_command "${PYTHON_BIN:-python3}")"
srun_bin="$(absolute_command srun)"
if (( ! prepare_only )); then sbatch_bin="$(absolute_command sbatch)"; fi
"$python_bin" -c 'import sys, numpy, h5py; print("Python:", sys.executable, "numpy:", numpy.__version__, "h5py:", h5py.__version__)' \
  || fail "Python 缺少 numpy/h5py 或无法导入: $python_bin；可通过 PYTHON_BIN 指定已有环境。"

jobs="${JOBS:-4}"
export REPEATS="${REPEATS:-1}"
export APRL_PROFILE_INTERVAL="${PROFILE:-${APRL_PROFILE_INTERVAL:-16}}"
export SRUN_MPI="${SRUN_MPI:-pmix_v3}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || fail 'JOBS 必须是正整数。'
[[ "$REPEATS" =~ ^[1-9][0-9]*$ ]] || fail 'REPEATS 必须是正整数。'
[[ "$APRL_PROFILE_INTERVAL" =~ ^[0-9]+$ ]] || fail 'PROFILE/APRL_PROFILE_INTERVAL 必须是非负整数。'
[[ "$SRUN_MPI" =~ ^[[:alnum:]_]+$ ]] || fail 'SRUN_MPI 必须是 srun --mpi=list 中的插件名称。'
if ! mpi_plugins="$("$srun_bin" --mpi=list 2>&1)"; then
  printf '%s\n' "$mpi_plugins" >&2
  fail '无法读取 Slurm MPI 插件列表；请先检查当前 srun/MPI 环境。'
fi
if ! grep -Eq "(^|[^[:alnum:]_])${SRUN_MPI}([^[:alnum:]_]|$)" <<< "$mpi_plugins"; then
  printf '当前 srun --mpi=list 输出:\n%s\n' "$mpi_plugins" >&2
  fail "未找到指定插件 ${SRUN_MPI}；请按上述实际名称设置 SRUN_MPI，例如 SRUN_MPI=pmi2，然后重试。"
fi
printf 'CMake: %s\nMPI C++: %s\nPython: %s\nSlurm MPI 插件: %s\n' "$cmake_bin" "$mpi_cxx" "$python_bin" "$SRUN_MPI"

before_build="${APRL_BASELINE_BUILD_DIR:-$source_dir/build-cluster-before}"
after_build="${APRL_CANDIDATE_BUILD_DIR:-$source_dir/build-cluster-after}"
submit_dir="${PERFORMANCE_SUBMIT_DIR:-$source_dir}"
mkdir -p -- "$before_build" "$after_build" "$submit_dir"
before_build="$(cd -- "$before_build" && pwd -P)"
after_build="$(cd -- "$after_build" && pwd -P)"
submit_dir="$(cd -- "$submit_dir" && pwd -P)"
[[ "$before_build" != "$after_build" ]] || fail '基线与新版构建目录必须不同。'
for build_dir in "$before_build" "$after_build"; do
  [[ "$build_dir" != "$source_dir" && "$build_dir" != "$old_source" ]] \
    || fail '构建目录不得等于新版或旧版源码根目录。'
done
common_options=(
  -DCMAKE_BUILD_TYPE=Release
  "-DCMAKE_CXX_COMPILER=$mpi_cxx"
  "-DCMAKE_PREFIX_PATH=$dependency_prefix"
  "-DHDF5_ROOT=$dependency_prefix"
  -DAPRL_REQUIRE_PARALLEL_HDF5=ON
  -DAPRL_BUILD_POSTPROCESSORS=OFF
  -DBUILD_TESTING=OFF
)
printf '[2/5] 编译旧版基线: %s\n' "$before_build"
"$cmake_bin" -S "$old_source" -B "$before_build" "${common_options[@]}"
"$cmake_bin" --build "$before_build" --target aprl -j "$jobs"
printf '[3/5] 编译新版及同一输入生成器: %s\n' "$after_build"
"$cmake_bin" -S "$source_dir" -B "$after_build" "${common_options[@]}" \
  -DAPRL_ENABLE_OPENMP=ON -DAPRL_BUILD_SASE_BUNCH_GENERATOR=ON
"$cmake_bin" --build "$after_build" --target aprl sase_bunch_generator -j "$jobs"

export APRL_SOURCE_DIR="$source_dir"
export APRL_BASELINE="$before_build/aprl"
export APRL_CANDIDATE="$after_build/aprl"
export SASE_GENERATOR="$after_build/sase_bunch_generator"
export PYTHON_BIN="$python_bin"
export RANK_COUNTS=1,2 THREAD_COUNTS=1,2
for executable in "$APRL_BASELINE" "$APRL_CANDIDATE" "$SASE_GENERATOR"; do
  [[ -x "$executable" ]] || fail "编译未产生可执行文件: $executable"
done
printf '[4/5] 对照已准备\n基线: %s\n新版: %s\n生成器: %s\n' "$APRL_BASELINE" "$APRL_CANDIDATE" "$SASE_GENERATOR"
printf 'MPI ranks=%s；每 rank 线程=%s；重复=%s；profile 间隔=%s\n' "$RANK_COUNTS" "$THREAD_COUNTS" "$REPEATS" "$APRL_PROFILE_INTERVAL"
printf '提交入口: %s\n提交工作目录及 Slurm 日志: %s\n结果目录: %s\n' \
  "$source_dir/verification/performance/run.sbatch" "$submit_dir" \
  "${PERFORMANCE_OUTPUT_DIR:-${source_dir}/work/performance/job-<作业号>}"
if (( prepare_only )); then
  printf '[5/5] --prepare-only 完成，未提交作业。审核后去掉该选项重新运行即可提交。\n'
  exit 0
fi

printf '[5/5] 提交小规模对照：分区=%s，1 节点，2 MPI slots，每 slot 2 CPU\n' "${PARTITION:-ssmb}"
cd -- "$submit_dir"
"$sbatch_bin" --export=ALL --partition="${PARTITION:-ssmb}" --nodes=1 --ntasks=2 \
  --cpus-per-task=2 --time=00:30:00 --chdir="$submit_dir" \
  --output="$submit_dir/aprl-performance-%j.out" --error="$submit_dir/aprl-performance-%j.err" \
  "$source_dir/verification/performance/run.sbatch"
