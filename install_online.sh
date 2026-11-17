#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./install_online.sh [options]

Inspect the current Linux environment, select one consistent platform MPI,
download verified third-party sources, and install APRL plus its private
runtime under the source tree. No root privilege or system package install
is used.

Options:
  --profile core|full    core, or core plus field post-processors (default: core)
  --prefix DIR           private install prefix (default: ./aprl-runtime)
  --jobs N               build parallelism (default: 8; at most 16 in a job)
  --mpi-module NAME      load this site MPI module before detection
  --mpi-run              perform the two-rank check during installation
  --skip-mpi-run         defer MPI checks even inside an allocation
  --skip-tests           skip APRL's lightweight required regressions
  --keep-build           retain third-party and APRL build trees
  --plan-only            inspect without downloading/building/installing
  -h, --help             show this help

Environment:
  APRL_MPI_MODULE        same as --mpi-module

On a managed multi-node cluster this script deliberately reuses the site MPI.
It locally builds CMake when needed, yaml-cpp, zlib, parallel HDF5, and FFTW
for the full profile. Source archives are SHA-256 verified and cached for
reuse.
EOF
}

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
profile="core"
prefix="$source_root/aprl-runtime"
jobs=""
mpi_module="${APRL_MPI_MODULE:-}"
mpi_check_mode="auto"
skip_tests=false
keep_build=false
plan_only=false
initial_conda_prefix="${CONDA_PREFIX:-}"

while (($#)); do
  case "$1" in
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --prefix)
      [[ $# -ge 2 ]] || { echo "--prefix needs a value" >&2; exit 2; }
      prefix="$2"; shift 2 ;;
    --jobs)
      [[ $# -ge 2 ]] || { echo "--jobs needs a value" >&2; exit 2; }
      jobs="$2"; shift 2 ;;
    --mpi-module)
      [[ $# -ge 2 ]] || { echo "--mpi-module needs a value" >&2; exit 2; }
      mpi_module="$2"; shift 2 ;;
    --mpi-run) mpi_check_mode="run"; shift ;;
    --skip-mpi-run) mpi_check_mode="skip"; shift ;;
    --skip-tests) skip_tests=true; shift ;;
    --keep-build) keep_build=true; shift ;;
    --plan-only) plan_only=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$profile" != "core" && "$profile" != "full" ]]; then
  echo "--profile must be core or full" >&2
  exit 2
fi
if [[ -n "$jobs" && ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "--jobs must be a positive integer" >&2
  exit 2
fi
if [[ -z "$jobs" ]]; then
  if [[ -n "${SLURM_JOB_ID:-}" &&
        "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
    detected_jobs="$SLURM_CPUS_PER_TASK"
    [[ "$detected_jobs" =~ ^[1-9][0-9]*$ ]] || detected_jobs=1
    if ((detected_jobs > 16)); then jobs=16; else jobs="$detected_jobs"; fi
  else
    jobs=8
  fi
fi

absolute_path() {
  local value="$1"
  if [[ "$value" == /* ]]; then
    realpath -m -- "$value"
  else
    realpath -m -- "$(pwd -P)/$value"
  fi
}

prefix="$(absolute_path "$prefix")"
if [[ -z "$prefix" || "$prefix" == "/" || "$prefix" == "$source_root" ]]; then
  echo "Refusing unsafe install prefix: $prefix" >&2
  exit 2
fi

initialize_modules() {
  if type module >/dev/null 2>&1; then return 0; fi
  local initialization=""
  for initialization in \
      /etc/profile.d/modules.sh \
      /usr/share/Modules/init/bash \
      /apps/support/tools/modules5.4.0/init/bash; do
    if [[ -r "$initialization" ]]; then
      # shellcheck source=/dev/null
      source "$initialization"
      type module >/dev/null 2>&1 && return 0
    fi
  done
  return 1
}

is_conda_path() {
  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    case "$1" in
      "$CONDA_PREFIX"|"$CONDA_PREFIX"/*) return 0 ;;
    esac
  fi
  case "$1" in
    *"/.conda/"*|*"/conda/"*|*"/miniconda/"*|*"/miniconda3/"*|\
    *"/anaconda/"*|*"/anaconda3/"*|*"/miniforge/"*|*"/miniforge3/"*|\
    *"/mambaforge/"*|*"/micromamba/"*) return 0 ;;
    *) return 1 ;;
  esac
}

# Prevent an activated Conda environment from leaking incompatible wrappers,
# headers or libraries into the platform-MPI build. Module paths and ordinary
# system paths are retained.
sanitize_colon_variable() {
  local variable_name="$1"
  local current_value="${!variable_name:-}"
  local clean_value=""
  local entry=""
  local entries=()
  IFS=':' read -r -a entries <<<"$current_value"
  for entry in "${entries[@]}"; do
    [[ -n "$entry" ]] || continue
    is_conda_path "$entry" && continue
    clean_value+="${clean_value:+:}$entry"
  done
  if [[ -n "$clean_value" ]]; then
    printf -v "$variable_name" '%s' "$clean_value"
    export "$variable_name"
  else
    unset "$variable_name"
  fi
}

for search_variable in PATH \
    LD_LIBRARY_PATH LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH \
    CMAKE_PREFIX_PATH PKG_CONFIG_PATH; do
  sanitize_colon_variable "$search_variable"
done
unset CC CXX CPP FC F77 F90 CFLAGS CXXFLAGS CPPFLAGS LDFLAGS
unset HDF5_ROOT HDF5_DIR YAML_CPP_DIR CMAKE_TOOLCHAIN_FILE
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL
unset _CE_CONDA _CE_M

valid_platform_mpi() {
  local cxx="$(command -v mpic++ 2>/dev/null || true)"
  local wrapper_dir=""
  local cc=""
  local launcher=""
  [[ -n "$cxx" && -x "$cxx" ]] || return 1
  wrapper_dir="$(dirname -- "$cxx")"
  cc="$wrapper_dir/mpicc"
  if [[ -x "$wrapper_dir/mpirun" ]]; then
    launcher="$wrapper_dir/mpirun"
  elif [[ -x "$wrapper_dir/mpiexec" ]]; then
    launcher="$wrapper_dir/mpiexec"
  fi
  [[ -x "$cc" && -n "$launcher" ]] || return 1
  is_conda_path "$cxx" && return 1
  is_conda_path "$cc" && return 1
  is_conda_path "$launcher" && return 1
  return 0
}

load_mpi_module() {
  local requested="$1"
  initialize_modules || {
    echo "Cannot load MPI module '$requested': module command unavailable." >&2
    return 1
  }
  module load "$requested" || return 1
  hash -r
  valid_platform_mpi || {
    echo "Module '$requested' did not provide non-Conda mpicc and mpic++." >&2
    return 1
  }
}

loaded_mpi_module() {
  local module_name=""
  local lower_name=""
  local loaded_modules=()
  IFS=':' read -r -a loaded_modules <<<"${LOADEDMODULES:-}"
  for module_name in "${loaded_modules[@]}"; do
    lower_name="${module_name,,}"
    case "$lower_name" in
      *openmpi*|*mpich*|*intelmpi*|*intel-mpi*|mpi/*|*/mpi/*)
        printf '%s\n' "$module_name"
        return 0 ;;
    esac
  done
  return 1
}

selected_module=""
if [[ -n "$mpi_module" ]]; then
  load_mpi_module "$mpi_module"
  selected_module="$mpi_module"
else
  detected_loaded_module="$(loaded_mpi_module || true)"
  if [[ -n "$detected_loaded_module" ]] && valid_platform_mpi; then
    selected_module="$detected_loaded_module"
  else
    # Known Yinxiu site profile. Other sites can pass --mpi-module.
    known_module="intel/openmpi/5.0.3/gcc8.5.0"
    if initialize_modules; then
      known_availability="$(module -t avail "$known_module" 2>&1 || true)"
      if grep -Fq "$known_module" <<<"$known_availability"; then
        if load_mpi_module "$known_module"; then
          selected_module="$known_module"
        fi
      fi
    fi
  fi
fi

if ! valid_platform_mpi; then
  current_mpi="$(command -v mpic++ 2>/dev/null || echo missing)"
  cat >&2 <<EOF
No usable platform MPI was selected (current mpic++: $current_mpi).
Load the cluster MPI module first, or run:
  ./install_online.sh --mpi-module MODULE_NAME
Conda MPI is intentionally rejected for a managed multi-node cluster.
EOF
  exit 1
fi

if [[ -z "$selected_module" ]]; then
  echo "WARNING: using the non-Conda MPI already on PATH; no MPI module was identified." >&2
  echo "Pass --mpi-module NAME on a module-managed cluster if this is not the site MPI." >&2
fi

mpi_cxx="$(command -v mpic++)"
mpi_wrapper_dir="$(dirname -- "$mpi_cxx")"
mpi_cc="$mpi_wrapper_dir/mpicc"
if [[ -x "$mpi_wrapper_dir/mpirun" ]]; then
  mpi_launcher="$mpi_wrapper_dir/mpirun"
else
  mpi_launcher="$mpi_wrapper_dir/mpiexec"
fi

if [[ "$mpi_check_mode" == "auto" ]]; then
  allocated_tasks="${SLURM_NTASKS:-${SLURM_NPROCS:-0}}"
  allocated_nodes="${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}}"
  if [[ -n "${SLURM_JOB_ID:-}" && "$allocated_tasks" =~ ^[0-9]+$ &&
        "$allocated_tasks" -ge 2 && "$allocated_nodes" =~ ^[0-9]+$ &&
        "$allocated_nodes" -eq 1 ]]; then
    run_mpi_check=true
  else
    run_mpi_check=false
  fi
elif [[ "$mpi_check_mode" == "run" ]]; then
  run_mpi_check=true
else
  run_mpi_check=false
fi

os_name="unknown"
os_version="unknown"
if [[ -r /etc/os-release ]]; then
  # shellcheck source=/dev/null
  source /etc/os-release
  os_name="${PRETTY_NAME:-${NAME:-unknown}}"
  os_version="${VERSION_ID:-unknown}"
fi

report_dir="$source_root/.cache/aprl-install"
report="$report_dir/environment-$(date -u +%Y%m%dT%H%M%SZ).txt"
mkdir -p -- "$report_dir"
if command -v flock >/dev/null 2>&1; then
  exec {online_lock_fd}>"$report_dir/install.lock"
  if ! flock -n "$online_lock_fd"; then
    echo "Another APRL online installation is using this source tree." >&2
    exit 1
  fi
fi
{
  echo "generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "source_root=$source_root"
  echo "host=$(hostname 2>/dev/null || echo unknown)"
  echo "os=$os_name"
  echo "os_version=$os_version"
  echo "architecture=$(uname -m)"
  echo "slurm_job_id=${SLURM_JOB_ID:-none}"
  echo "slurm_tasks=${SLURM_NTASKS:-${SLURM_NPROCS:-unknown}}"
  echo "slurm_nodes=${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-unknown}}"
  echo "profile=$profile"
  echo "prefix=$prefix"
  echo "jobs=$jobs"
  echo "conda_environment_removed=$([[ -n "$initial_conda_prefix" ]] && echo yes || echo no)"
  echo "selected_module=${selected_module:-preloaded-environment}"
  echo "mpi_cc=$mpi_cc"
  echo "mpi_cxx=$mpi_cxx"
  echo "mpi_launcher=$mpi_launcher"
  echo "mpi_checks=$($run_mpi_check && echo run || echo deferred)"
  echo "mpi_cxx_show=$($mpi_cxx --showme:command 2>/dev/null || $mpi_cxx -show 2>/dev/null || true)"
  echo "mpi_launcher_version_begin"
  "$mpi_launcher" --version 2>&1 | head -n 6 || true
  echo "mpi_launcher_version_end"
  if type module >/dev/null 2>&1; then
    echo "modules_begin"
    module list 2>&1 || true
    echo "modules_end"
  fi
} | tee "$report"

cat <<EOF
=== APRL online installation plan ===
The platform MPI above is retained for Slurm/network integration.
Verified sources will be downloaded for missing private dependencies.
Environment report: $report
EOF

if $plan_only; then
  "$source_root/bootstrap_dependencies.sh" \
    --prefix "$prefix" --profile "$profile" --jobs "$jobs" \
    --mpi-cc "$mpi_cc" --mpi-cxx "$mpi_cxx" \
    --mpi-launcher "$mpi_launcher" --plan-only
  echo "RESULT=PLAN: no download, build, or installation was performed."
  exit 0
fi

install_args=(
  --bootstrap-deps
  --profile "$profile"
  --prefix "$prefix"
  --jobs "$jobs"
  --mpi-cc "$mpi_cc"
  --mpi-cxx "$mpi_cxx"
  --mpi-launcher "$mpi_launcher"
)
$run_mpi_check || install_args+=(--skip-mpi-run)
$skip_tests && install_args+=(--skip-tests)
$keep_build && install_args+=(--keep-build)

"$source_root/install.sh" "${install_args[@]}"

installed_aprl="$prefix/bin/aprl"
[[ -x "$installed_aprl" ]] || {
  echo "Installation returned without the APRL executable: $installed_aprl" >&2
  exit 1
}
if command -v ldd >/dev/null 2>&1; then
  ldd_output="$(ldd "$installed_aprl" 2>&1)"
  missing="$(printf '%s\n' "$ldd_output" | grep 'not found' || true)"
  if [[ -n "$missing" ]]; then
    echo "Installed APRL still has unresolved runtime libraries:" >&2
    echo "$missing" >&2
    exit 1
  fi
  for private_library in libhdf5 libyaml-cpp; do
    resolved_library="$(printf '%s\n' "$ldd_output" |
      awk -v wanted="$private_library" \
        '$1 ~ wanted && $2 == "=>" {print $3; exit}')"
    if [[ -n "$resolved_library" ]]; then
      resolved_library="$(realpath -m -- "$resolved_library")"
    fi
    case "$resolved_library" in
      "$prefix/lib/"*|"$prefix/lib64/"*) ;;
      "")
        echo "Installed APRL does not resolve $private_library through ldd." >&2
        exit 1 ;;
      *)
        echo "Installed APRL resolved $private_library outside the private prefix:" >&2
        echo "  $resolved_library" >&2
        exit 1 ;;
    esac
  done
fi

compression_probe="$report_dir/hdf5-deflate-probe-$$.h5"
rm -f -- "$compression_probe"
echo "Verifying private HDF5 gzip/deflate support"
if ! "$prefix/bin/sase_bunch_generator" \
    --output "$compression_probe" \
    --electrons 100000 \
    --length-m 8.0e-6 \
    --resonant-wavelength-m 1.0e-6 \
    --gamma 50.0 \
    --sigma-x-m 2.0e-5 \
    --sigma-y-m 3.0e-5 \
    --beamlets-per-wavelength 4 \
    --compression 1 >/dev/null; then
  rm -f -- "$compression_probe"
  echo "The installed HDF5 cannot create a compressed APRL particle file." >&2
  exit 1
fi
rm -f -- "$compression_probe"

cat <<EOF
RESULT=PASS: APRL online installation is complete.
  executable: $installed_aprl
  runtime environment: $prefix/share/aprl/runtime-env.sh
  dependency lock: $prefix/share/aprl/dependencies.lock
  environment report: $report

For cluster jobs, load the same MPI module before sourcing runtime-env.sh.
EOF
if [[ -n "$selected_module" ]]; then
  echo "  module load $selected_module"
  echo "  source $prefix/share/aprl/runtime-env.sh"
fi
