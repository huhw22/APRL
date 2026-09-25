#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./install.sh [options]

Validate dependencies, build, test and install. Successful installation
removes compilation products by default, but never removes source files.

Options:
  --prefix DIR             install prefix (default: /usr/local)
  --build-dir DIR          disposable build root (default: build-install)
  --profile core|full      main program only, or include post-processors
  --jobs N                 parallel build jobs
  --mpi-cxx COMMAND        MPI C++ wrapper (default: $MPICXX or mpic++)
  --mpi-launcher COMMAND   launcher accepting '-n 2' (default: mpiexec)
  --generator NAME         CMake generator (default: Ninja when available)
  --sdds-root DIR          use and validate an official built SDDS tree
  --require-sdds           fail if the native SDDS converter is unavailable
  --allow-serial-hdf5      single-rank/local build only
  --skip-mpi-run           skip dependency check's two-rank launch
  --skip-tests             do not run the required lightweight regressions
  --keep-build             retain build products after successful install
  --cmake-arg ARG          extra CMake argument; may be repeated
  -h, --help               show this help
EOF
}

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
prefix="/usr/local"
build_dir="$source_root/build-install"
profile="full"
jobs=""
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
      prefix="$2"; shift 2 ;;
    --build-dir)
      [[ $# -ge 2 ]] || { echo "--build-dir needs a value" >&2; exit 2; }
      build_dir="$2"; shift 2 ;;
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --jobs)
      [[ $# -ge 2 ]] || { echo "--jobs needs a value" >&2; exit 2; }
      jobs="$2"; shift 2 ;;
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
if [[ -n "$jobs" && ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "--jobs must be a positive integer" >&2
  exit 2
fi
if [[ -z "$jobs" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    jobs="$(nproc)"
  else
    jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
  fi
fi
if [[ -z "$generator" ]] && command -v ninja >/dev/null 2>&1; then
  generator="Ninja"
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
  -DCMAKE_INSTALL_PREFIX="$prefix"
  -DHDF5_PREFER_PARALLEL=TRUE
)

root_configure=(cmake -S "$source_root" -B "$build_dir/core")
root_configure+=("${generator_args[@]}")
root_configure+=("${common_cmake_args[@]}")
root_configure+=(
  -DBUILD_TESTING="$testing"
  -DFEL_REQUIRE_PARALLEL_HDF5="$parallel_requirement"
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
    ctest --test-dir "$build_dir/core" --output-on-failure -L required \
      -E '^required\.hdf5-versions$'
  else
    cmake --build "$build_dir/core" --target verify_required
  fi
fi

build_roots=("$build_dir/core")
if [[ "$profile" == "full" ]]; then
  postprocessors=(trajectory_radiation field_reconstruction
    field_plane_analysis field_power_compare energy_closure)
  for component in "${postprocessors[@]}"; do
    component_build="$build_dir/postprocess-$component"
    echo "Configuring post-processor: $component"
    component_configure=(cmake
      -S "$source_root/postprocess/$component"
      -B "$component_build")
    component_configure+=("${generator_args[@]}")
    component_configure+=("${common_cmake_args[@]}")
    component_configure+=("${extra_cmake_args[@]}")
    "${component_configure[@]}"
    cmake --build "$component_build" --parallel "$jobs"
    build_roots+=("$component_build")
  done
fi

receipt_dir="$prefix/share/unnamed_fel_program"
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
cmake -E copy "$build_dir/install-manifest.sorted.txt" "$receipt"

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

echo "Installation complete"
echo "  prefix: $prefix"
echo "  profile: $profile"
echo "  manifest: $receipt"
echo "Add $prefix/bin to PATH if the prefix is not already searched."

if $keep_build; then
  echo "Build products retained: $build_dir"
else
  cmake -E remove_directory "$build_dir"
  echo "Compilation products removed after successful installation."
fi
