#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./install.sh [options]

Validate dependencies, build, test and install. Successful installation
removes compilation products by default, but never removes source files.

Options:
  --prefix DIR             install prefix (/usr/local; with bootstrap: ./aprl-runtime)
  --build-dir DIR          disposable build root (default: build-install)
  --profile core|full      main program only, or include post-processors
  --jobs N                 parallel build jobs
  --bootstrap-deps         build verified dependencies into --prefix
  --bundle-mpi             with --bootstrap-deps, build a private MPICH
  --dependency-cache DIR   verified source archive cache
  --bootstrap-build-dir DIR disposable third-party build tree
  --offline-deps           bootstrap only from the verified local cache
  --mpi-cc COMMAND         MPI C wrapper used to build parallel HDF5
  --mpi-cxx COMMAND        MPI C++ wrapper (default: $MPICXX or mpic++)
  --mpi-launcher COMMAND   launcher accepting '-n 2' (default: mpiexec)
  --generator NAME         CMake generator (default: Ninja when available)
  --sdds-root DIR          use and validate an official built SDDS tree
  --require-sdds           fail if the native SDDS converter is unavailable
  --allow-serial-hdf5      single-rank/local build only
  --skip-mpi-run           skip dependency and regression MPI launches
  --skip-tests             do not run the required lightweight regressions
  --keep-build             retain build products after successful install
  --cmake-arg ARG          extra CMake argument; may be repeated
  -h, --help               show this help
EOF
}

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
prefix="/usr/local"
prefix_explicit=false
build_dir="$source_root/build-install"
profile="full"
jobs=""
bootstrap_deps=false
bundle_mpi=false
dependency_cache="$source_root/.cache/aprl-dependencies"
bootstrap_build_dir="$source_root/build-bootstrap"
offline_deps=false
mpi_cc="${MPICC:-}"
mpi_cxx="${MPICXX:-mpic++}"
mpi_launcher="${MPIEXEC:-mpiexec}"
generator=""
sdds_root="${SDDS_ROOT:-}"
require_sdds=false
allow_serial_hdf5=false
skip_mpi_run=false
skip_tests=false
keep_build=false
extra_cmake_args=()

while (($#)); do
  case "$1" in
    --prefix)
      [[ $# -ge 2 ]] || { echo "--prefix needs a value" >&2; exit 2; }
      prefix="$2"; prefix_explicit=true; shift 2 ;;
    --build-dir)
      [[ $# -ge 2 ]] || { echo "--build-dir needs a value" >&2; exit 2; }
      build_dir="$2"; shift 2 ;;
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --jobs)
      [[ $# -ge 2 ]] || { echo "--jobs needs a value" >&2; exit 2; }
      jobs="$2"; shift 2 ;;
    --bootstrap-deps) bootstrap_deps=true; shift ;;
    --bundle-mpi) bundle_mpi=true; shift ;;
    --dependency-cache)
      [[ $# -ge 2 ]] || { echo "--dependency-cache needs a value" >&2; exit 2; }
      dependency_cache="$2"; shift 2 ;;
    --bootstrap-build-dir)
      [[ $# -ge 2 ]] || { echo "--bootstrap-build-dir needs a value" >&2; exit 2; }
      bootstrap_build_dir="$2"; shift 2 ;;
    --offline-deps) offline_deps=true; shift ;;
    --mpi-cc)
      [[ $# -ge 2 ]] || { echo "--mpi-cc needs a value" >&2; exit 2; }
      mpi_cc="$2"; shift 2 ;;
    --mpi-cxx)
      [[ $# -ge 2 ]] || { echo "--mpi-cxx needs a value" >&2; exit 2; }
      mpi_cxx="$2"; shift 2 ;;
    --mpi-launcher)
      [[ $# -ge 2 ]] || { echo "--mpi-launcher needs a value" >&2; exit 2; }
      mpi_launcher="$2"; shift 2 ;;
    --generator)
      [[ $# -ge 2 ]] || { echo "--generator needs a value" >&2; exit 2; }
      generator="$2"; shift 2 ;;
    --sdds-root)
      [[ $# -ge 2 ]] || { echo "--sdds-root needs a value" >&2; exit 2; }
      sdds_root="$2"; shift 2 ;;
    --require-sdds) require_sdds=true; shift ;;
    --allow-serial-hdf5) allow_serial_hdf5=true; shift ;;
    --skip-mpi-run) skip_mpi_run=true; shift ;;
    --skip-tests) skip_tests=true; shift ;;
    --keep-build) keep_build=true; shift ;;
    --cmake-arg)
      [[ $# -ge 2 ]] || { echo "--cmake-arg needs a value" >&2; exit 2; }
      extra_cmake_args+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$profile" != "core" && "$profile" != "full" ]]; then
  echo "--profile must be core or full" >&2
  exit 2
fi
if $require_sdds && [[ -z "$sdds_root" ]]; then
  echo "--require-sdds needs --sdds-root DIR (or SDDS_ROOT)" >&2
  exit 2
fi
if $bundle_mpi && ! $bootstrap_deps; then
  echo "--bundle-mpi requires --bootstrap-deps" >&2
  exit 2
fi
if $offline_deps && ! $bootstrap_deps; then
  echo "--offline-deps requires --bootstrap-deps" >&2
  exit 2
fi
if $bootstrap_deps && ! $prefix_explicit; then
  prefix="$source_root/aprl-runtime"
fi
if [[ -n "$jobs" && ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "--jobs must be a positive integer" >&2
  exit 2
fi
if [[ -z "$jobs" ]]; then
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" &&
        "${SLURM_CPUS_PER_TASK}" =~ ^[1-9][0-9]*$ ]]; then
    jobs="$SLURM_CPUS_PER_TASK"
  elif command -v nproc >/dev/null 2>&1; then
    jobs="$(nproc)"
  else
    jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
  fi
  if ((jobs > 16)); then jobs=16; fi
fi

absolute_path() {
  local value="$1"
  local candidate=""
  if [[ "$value" == /* ]]; then
    candidate="$value"
  else
    candidate="$(pwd -P)/$value"
  fi
  realpath -m -- "$candidate"
}

prefix="$(absolute_path "$prefix")"
build_dir="$(absolute_path "$build_dir")"
dependency_cache="$(absolute_path "$dependency_cache")"
bootstrap_build_dir="$(absolute_path "$bootstrap_build_dir")"
if [[ -z "$prefix" || "$prefix" == "/" ]]; then
  echo "Refusing an empty or root install prefix" >&2
  exit 2
fi
if [[ -z "$build_dir" || "$build_dir" == "/" ||
      "$build_dir" == "$source_root" ]]; then
  echo "Refusing unsafe build directory: $build_dir" >&2
  exit 2
fi
case "$source_root/" in
  "$build_dir/"*)
    echo "Build directory may not contain the source tree: $build_dir" >&2
    exit 2 ;;
esac
case "$prefix/" in
  "$build_dir/"*)
    echo "Install prefix may not be inside the disposable build tree: $prefix" >&2
    exit 2 ;;
esac
case "$build_dir/" in
  "$prefix/"*)
    echo "Disposable build tree may not be inside the install prefix: $build_dir" >&2
    exit 2 ;;
esac

mkdir -p -- "$(dirname -- "$build_dir")"
if command -v flock >/dev/null 2>&1; then
  install_lock_dir="$source_root/.cache/aprl-locks"
  mkdir -p -- "$install_lock_dir"
  exec {install_lock_fd}>"$install_lock_dir/install.lock"
  if ! flock -n "$install_lock_fd"; then
    echo "Another APRL installation is using this source tree." >&2
    exit 1
  fi
fi

if $bootstrap_deps; then
  if ! $bundle_mpi; then
    if [[ "$mpi_cxx" == */* ]]; then
      [[ -x "$mpi_cxx" ]] || {
        echo "MPI C++ wrapper is not executable: $mpi_cxx" >&2
        exit 1
      }
      mpi_cxx="$(command -v -- "$mpi_cxx")"
    else
      mpi_cxx="$(command -v -- "$mpi_cxx" 2>/dev/null || true)"
    fi
    if [[ -z "$mpi_cxx" ]]; then
      echo "No platform MPI C++ wrapper found; load the site MPI or add --bundle-mpi." >&2
      exit 1
    fi
    if [[ -n "$mpi_cc" ]]; then
      if [[ "$mpi_cc" == */* ]]; then
        [[ -x "$mpi_cc" ]] || {
          echo "MPI C wrapper is not executable: $mpi_cc" >&2
          exit 1
        }
        mpi_cc="$(command -v -- "$mpi_cc")"
      else
        mpi_cc="$(command -v -- "$mpi_cc" 2>/dev/null || true)"
      fi
    else
      mpi_wrapper_dir="$(dirname -- "$mpi_cxx")"
      for mpi_cc_candidate in mpicc mpiicc mpcc; do
        if [[ -x "$mpi_wrapper_dir/$mpi_cc_candidate" ]]; then
          mpi_cc="$mpi_wrapper_dir/$mpi_cc_candidate"
          break
        fi
      done
    fi
    if [[ -z "$mpi_cc" ]]; then
      echo "Could not derive the MPI C wrapper next to $mpi_cxx; pass --mpi-cc." >&2
      exit 1
    fi
    if [[ "$mpi_launcher" == */* ]]; then
      [[ -x "$mpi_launcher" ]] || {
        echo "MPI launcher is not executable: $mpi_launcher" >&2
        exit 1
      }
      mpi_launcher="$(command -v -- "$mpi_launcher")"
    elif command -v -- "$mpi_launcher" >/dev/null 2>&1; then
      mpi_launcher="$(command -v -- "$mpi_launcher")"
    else
      echo "MPI launcher not found: $mpi_launcher" >&2
      exit 1
    fi
  fi

  bootstrap_args=(
    --prefix "$prefix"
    --build-dir "$bootstrap_build_dir"
    --cache-dir "$dependency_cache"
    --profile "$profile"
    --jobs "$jobs"
    --mpi-cxx "$mpi_cxx"
    --mpi-launcher "$mpi_launcher"
  )
  [[ -n "$mpi_cc" ]] && bootstrap_args+=(--mpi-cc "$mpi_cc")
  $bundle_mpi && bootstrap_args+=(--bundle-mpi)
  $offline_deps && bootstrap_args+=(--offline)
  $keep_build && bootstrap_args+=(--keep-build)
  "$source_root/bootstrap_dependencies.sh" "${bootstrap_args[@]}"

  if $bundle_mpi; then
    mpi_cc="$prefix/bin/mpicc"
    mpi_cxx="$prefix/bin/mpic++"
    mpi_launcher="$prefix/bin/mpiexec"
  fi
  # shellcheck source=/dev/null
  source "$prefix/share/aprl/runtime-env.sh"
  export HDF5_ROOT="$prefix"
fi

if [[ -z "$generator" ]] && command -v ninja >/dev/null 2>&1; then
  generator="Ninja"
fi

dependency_report="${build_dir}.dependency-report.log"
dependency_args=(--profile "$profile" --report "$dependency_report"
  --mpi-cxx "$mpi_cxx" --mpi-launcher "$mpi_launcher")
$allow_serial_hdf5 && dependency_args+=(--allow-serial-hdf5)
$skip_mpi_run && dependency_args+=(--skip-mpi-run)
[[ -n "$sdds_root" ]] && dependency_args+=(--sdds-root "$sdds_root")
$require_sdds && dependency_args+=(--require-sdds)
"$source_root/check_dependencies.sh" "${dependency_args[@]}"

echo "Removing any stale disposable build tree: $build_dir"
cmake -E remove_directory "$build_dir"
cmake -E make_directory "$build_dir"

generator_args=()
[[ -n "$generator" ]] && generator_args=(-G "$generator")
parallel_requirement=ON
$allow_serial_hdf5 && parallel_requirement=OFF
testing=ON
$skip_tests && testing=OFF

common_cmake_args=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_CXX_COMPILER="$mpi_cxx"
  -DMPIEXEC_EXECUTABLE="$mpi_launcher"
  -DCMAKE_INSTALL_PREFIX="$prefix"
  -DHDF5_PREFER_PARALLEL=TRUE
)
if $bootstrap_deps; then
  common_cmake_args+=(
    -DCMAKE_C_COMPILER="$mpi_cc"
    -DCMAKE_PREFIX_PATH="$prefix"
    -DHDF5_ROOT="$prefix"
    -DAPRL_PORTABLE_INSTALL=ON
  )
fi

root_configure=(cmake -S "$source_root" -B "$build_dir/core")
root_configure+=("${generator_args[@]}")
root_configure+=("${common_cmake_args[@]}")
root_configure+=(
  -DBUILD_TESTING="$testing"
  -DAPRL_REQUIRE_PARALLEL_HDF5="$parallel_requirement"
  -DAPRL_BUILD_POSTPROCESSORS=$([[ "$profile" == "full" ]] && echo ON || echo OFF)
)
[[ -n "$sdds_root" ]] && root_configure+=("-DSDDS_ROOT=$sdds_root")
root_configure+=("${extra_cmake_args[@]}")

echo "Configuring core build"
"${root_configure[@]}"
cmake --build "$build_dir/core" --parallel "$jobs"

if $require_sdds && [[ ! -x "$build_dir/core/elegant_sdds_to_hdf5" ]]; then
  echo "The native Elegant SDDS converter was required but was not built." >&2
  exit 1
elif [[ -n "$sdds_root" && ! -x "$build_dir/core/elegant_sdds_to_hdf5" ]]; then
  echo "WARNING: SDDS_ROOT was supplied but the converter was not built." >&2
fi

if ! $skip_tests; then
  if $allow_serial_hdf5; then
    echo "Running required tests except the parallel-HDF5 version test"
    test_args=(--test-dir "$build_dir/core" --output-on-failure -L required
      -E '^required\.hdf5-versions$')
    $skip_mpi_run && test_args+=(-LE mpi)
    ctest "${test_args[@]}"
  elif $skip_mpi_run; then
    echo "Running required non-MPI tests; MPI tests are deferred"
    cmake --build "$build_dir/core" --target verify_required_no_mpi
  else
    cmake --build "$build_dir/core" --target verify_required
  fi
fi

build_roots=("$build_dir/core")

receipt_dir="$prefix/share/aprl"
receipt="$receipt_dir/install-manifest.txt"
install_info="$receipt_dir/install-info.txt"
previous_manifest=""
if [[ -f "$receipt" ]]; then
  previous_manifest="$build_dir/previous-install-manifest.txt"
  cp -- "$receipt" "$previous_manifest"
  while IFS= read -r previous_path || [[ -n "$previous_path" ]]; do
    [[ -n "$previous_path" ]] || continue
    normalized_previous="$(realpath -m -- "$previous_path")"
    case "$normalized_previous" in
      "$prefix"/*) ;;
      *)
        echo "Existing receipt contains a path outside this prefix: $previous_path" >&2
        exit 1 ;;
    esac
  done < "$previous_manifest"
  echo "Existing installation receipt found; obsolete files will be removed"
  echo "only after the replacement installation and new receipt succeed."
fi

combined_manifest="$build_dir/install-manifest.txt"
: > "$combined_manifest"
for component_build in "${build_roots[@]}"; do
  cmake --install "$component_build"
  if [[ ! -s "$component_build/install_manifest.txt" ]]; then
    echo "Missing CMake install manifest: $component_build" >&2
    exit 1
  fi
  {
    cat "$component_build/install_manifest.txt"
    printf '\n'
  } >> "$combined_manifest"
done

cmake -E make_directory "$receipt_dir"
source_revision="unknown"
if command -v git >/dev/null 2>&1; then
  source_revision="$(git -C "$source_root" rev-parse HEAD 2>/dev/null || echo unknown)"
fi
cat > "$build_dir/install-info.txt" <<EOF
installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
source_root=$source_root
source_revision=$source_revision
profile=$profile
prefix=$prefix
mpi_cxx=$mpi_cxx
mpi_cc=${mpi_cc:-auto}
mpi_launcher=$mpi_launcher
mpi_tests_run=$(! $skip_mpi_run && echo yes || echo no)
bootstrapped_dependencies=$bootstrap_deps
bundled_mpi=$bundle_mpi
parallel_hdf5_required=$parallel_requirement
sdds_converter=$([[ -x "$build_dir/core/elegant_sdds_to_hdf5" ]] && echo yes || echo no)
EOF
cmake -E copy "$build_dir/install-info.txt" "$install_info"
printf '%s\n%s\n' "$install_info" "$receipt" >> "$combined_manifest"
sed '/^$/d' "$combined_manifest" |
  sort -u > "$build_dir/install-manifest.sorted.txt"
while IFS= read -r installed_path || [[ -n "$installed_path" ]]; do
  [[ -n "$installed_path" ]] || continue
  normalized_installed="$(realpath -m -- "$installed_path")"
  case "$normalized_installed" in
    "$prefix"/*) ;;
    *)
      echo "New install manifest escapes the prefix: $installed_path" >&2
      exit 1 ;;
  esac
done < "$build_dir/install-manifest.sorted.txt"
if $bootstrap_deps && command -v ldd >/dev/null 2>&1; then
  echo "Verifying installed runtime library resolution"
  runtime_missing=0
  runtime_ld_path=""
  IFS=':' read -r -a runtime_library_dirs <<<"${LD_LIBRARY_PATH:-}"
  for runtime_library_dir in "${runtime_library_dirs[@]}"; do
    [[ -n "$runtime_library_dir" ]] || continue
    case "$runtime_library_dir" in
      "$prefix/lib"|"$prefix/lib64") continue ;;
    esac
    runtime_ld_path+="${runtime_ld_path:+:}$runtime_library_dir"
  done
  while IFS= read -r installed_path || [[ -n "$installed_path" ]]; do
    case "$installed_path" in
      "$prefix/bin/"*) ;;
      *) continue ;;
    esac
    [[ -x "$installed_path" ]] || continue
    if [[ -n "$runtime_ld_path" ]]; then
      ldd_output="$(env LD_LIBRARY_PATH="$runtime_ld_path" \
        ldd "$installed_path" 2>&1 || true)"
    else
      ldd_output="$(env -u LD_LIBRARY_PATH ldd "$installed_path" 2>&1 || true)"
    fi
    if grep -q 'not found' <<<"$ldd_output"; then
      echo "Unresolved libraries for $installed_path:" >&2
      echo "$ldd_output" >&2
      runtime_missing=1
    fi
  done < "$build_dir/install-manifest.sorted.txt"
  if ((runtime_missing)); then
    echo "Installation completed, but its local runtime is incomplete." >&2
    exit 1
  fi
  echo "Installed runtime libraries resolved without 'not found' entries."
fi

if [[ -n "$previous_manifest" ]]; then
  while IFS= read -r previous_path || [[ -n "$previous_path" ]]; do
    [[ -n "$previous_path" ]] || continue
    normalized_previous="$(realpath -m -- "$previous_path")"
    if ! grep -Fxq -- "$previous_path" \
         "$build_dir/install-manifest.sorted.txt"; then
      if [[ -f "$normalized_previous" || -L "$normalized_previous" ]]; then
        rm -f -- "$normalized_previous"
        echo "Removed obsolete installed file: $normalized_previous"
      fi
    fi
  done < "$previous_manifest"
fi

# Publish the new receipt only after runtime verification and obsolete-file
# cleanup succeed.  Rename within the receipt directory so replacement is
# atomic and a failed update leaves the previous installation auditable.
receipt_temp="$receipt_dir/.install-manifest.txt.tmp.$$"
cp -- "$build_dir/install-manifest.sorted.txt" "$receipt_temp"
mv -f -- "$receipt_temp" "$receipt"

echo "Installation complete"
echo "  prefix: $prefix"
echo "  profile: $profile"
echo "  manifest: $receipt"
echo "Add $prefix/bin to PATH if the prefix is not already searched."
if $bootstrap_deps; then
  echo "Portable dependency environment: $prefix/share/aprl/runtime-env.sh"
  echo "Installed executables use relative RPATH for $prefix/lib and lib64."
fi

if $keep_build; then
  echo "Build products retained: $build_dir"
else
  cmake -E remove_directory "$build_dir"
  echo "Compilation products removed after successful installation."
fi
