#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./bootstrap_dependencies.sh [options]

Download, verify, build and install APRL's third-party dependencies into a
user-owned Linux prefix.  Nothing is installed into the operating system.

By default the platform MPI is reused and parallel HDF5 is compiled against
that exact MPI wrapper.  --bundle-mpi provides a portable single-node/
workstation fallback; a managed multi-node cluster should normally use its
site MPI for scheduler and high-speed-network integration.

Options:
  --prefix DIR          dependency/runtime prefix (default: ./aprl-runtime)
  --build-dir DIR       disposable dependency build tree
  --cache-dir DIR       verified download cache
  --profile core|full   full additionally builds FFTW3 + threads
  --jobs N              parallel build jobs
  --mpi-cc COMMAND      platform MPI C wrapper (default: derive from mpic++)
  --mpi-cxx COMMAND     platform MPI C++ wrapper (default: $MPICXX or mpic++)
  --mpi-launcher CMD    recorded runtime launcher (default: $MPIEXEC or mpiexec)
  --bundle-mpi          build a private MPICH instead of using platform MPI
  --offline             use only already verified cache files
  --keep-build          retain dependency build products
  --plan-only           validate arguments and print the dependency plan
  -h, --help            show this help

Pinned sources may be overridden for controlled mirrors with matching
APRL_CMAKE_*, APRL_YAML_CPP_*, APRL_HDF5_*, APRL_FFTW_* or APRL_MPICH_*
URL/SHA256 environment variables.  The selected values are recorded in the
generated dependency lock file.
EOF
}

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
prefix="$source_root/aprl-runtime"
build_dir="$source_root/build-bootstrap"
cache_dir="$source_root/.cache/aprl-dependencies"
profile="full"
jobs=""
mpi_cc="${MPICC:-}"
mpi_cxx="${MPICXX:-mpic++}"
mpi_launcher="${MPIEXEC:-mpiexec}"
bundle_mpi=false
offline=false
keep_build=false
plan_only=false

# Versions and hashes are deliberately pinned.  A changed mirror is accepted
# only when it provides byte-identical content.
cmake_version="3.31.8"
yaml_cpp_version="0.8.0"
hdf5_version="1.14.6"
fftw_version="3.3.10"
mpich_version="5.0.2"

yaml_cpp_url="${APRL_YAML_CPP_URL:-https://github.com/jbeder/yaml-cpp/archive/refs/tags/${yaml_cpp_version}.tar.gz}"
yaml_cpp_sha256="${APRL_YAML_CPP_SHA256:-fbe74bbdcee21d656715688706da3c8becfd946d92cd44705cc6098bb23b3a16}"
hdf5_url="${APRL_HDF5_URL:-https://github.com/HDFGroup/hdf5/releases/download/hdf5_${hdf5_version}/hdf5-${hdf5_version}.tar.gz}"
hdf5_sha256="${APRL_HDF5_SHA256:-e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b}"
fftw_url="${APRL_FFTW_URL:-https://fftw.org/pub/fftw/fftw-${fftw_version}.tar.gz}"
fftw_sha256="${APRL_FFTW_SHA256:-56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467}"
mpich_url="${APRL_MPICH_URL:-https://www.mpich.org/static/downloads/${mpich_version}/mpich-${mpich_version}.tar.gz}"
mpich_sha256="${APRL_MPICH_SHA256:-928c2f18d350a91443fe8024ad01ce2c6009e9c551e0a73e797fc567f119137d}"

while (($#)); do
  case "$1" in
    --prefix)
      [[ $# -ge 2 ]] || { echo "--prefix needs a value" >&2; exit 2; }
      prefix="$2"; shift 2 ;;
    --build-dir)
      [[ $# -ge 2 ]] || { echo "--build-dir needs a value" >&2; exit 2; }
      build_dir="$2"; shift 2 ;;
    --cache-dir)
      [[ $# -ge 2 ]] || { echo "--cache-dir needs a value" >&2; exit 2; }
      cache_dir="$2"; shift 2 ;;
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --jobs)
      [[ $# -ge 2 ]] || { echo "--jobs needs a value" >&2; exit 2; }
      jobs="$2"; shift 2 ;;
    --mpi-cc)
      [[ $# -ge 2 ]] || { echo "--mpi-cc needs a value" >&2; exit 2; }
      mpi_cc="$2"; shift 2 ;;
    --mpi-cxx)
      [[ $# -ge 2 ]] || { echo "--mpi-cxx needs a value" >&2; exit 2; }
      mpi_cxx="$2"; shift 2 ;;
    --mpi-launcher)
      [[ $# -ge 2 ]] || { echo "--mpi-launcher needs a value" >&2; exit 2; }
      mpi_launcher="$2"; shift 2 ;;
    --bundle-mpi) bundle_mpi=true; shift ;;
    --offline) offline=true; shift ;;
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
  if [[ -n "${SLURM_CPUS_PER_TASK:-}" &&
        "${SLURM_CPUS_PER_TASK}" =~ ^[1-9][0-9]*$ ]]; then
    jobs="$SLURM_CPUS_PER_TASK"
  elif command -v nproc >/dev/null 2>&1; then
    jobs="$(nproc)"
  else
    jobs=1
  fi
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
cache_dir="$(absolute_path "$cache_dir")"
for safe_path in "$prefix" "$build_dir" "$cache_dir"; do
  if [[ -z "$safe_path" || "$safe_path" == "/" ||
        "$safe_path" == "$source_root" ]]; then
    echo "Refusing unsafe bootstrap path: $safe_path" >&2
    exit 2
  fi
done
for retained_path in "$source_root" "$prefix" "$cache_dir"; do
  case "$retained_path/" in
    "$build_dir/"*)
      echo "Build directory may not contain retained path: $retained_path" >&2
      exit 2 ;;
  esac
done

host_os="$(uname -s)"
if [[ "$host_os" != "Linux" ]]; then
  echo "The local dependency bootstrap currently supports Linux only." >&2
  exit 1
fi
host_arch="$(uname -m)"
case "$host_arch" in
  x86_64)
    cmake_platform="linux-x86_64"
    cmake_sha256="${APRL_CMAKE_SHA256:-630615d8e98ac33eba7fbe472626dff5c899c85af3c024585ae109166a6909d0}"
    ;;
  aarch64|arm64)
    cmake_platform="linux-aarch64"
    cmake_sha256="${APRL_CMAKE_SHA256:-609735983e3bdf24b6ab379d918458d64196fe72b98226f62dd5e9fe7b2997cc}"
    ;;
  *)
    cmake_platform=""
    cmake_sha256="${APRL_CMAKE_SHA256:-}"
    ;;
esac
cmake_url="${APRL_CMAKE_URL:-https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}-${cmake_platform}.tar.gz}"

cat <<EOF
=== APRL local dependency bootstrap ===
source_root=$source_root
prefix=$prefix
build_dir=$build_dir
cache_dir=$cache_dir
profile=$profile
jobs=$jobs
mpi_provider=$($bundle_mpi && echo bundled-mpich || echo platform)
mpi_cc_request=${mpi_cc:-auto}
mpi_cxx_request=$mpi_cxx
mpi_launcher_request=$mpi_launcher
offline=$offline
cmake=$cmake_version
yaml_cpp=$yaml_cpp_version
hdf5=$hdf5_version
fftw=$([[ "$profile" == "full" ]] && echo "$fftw_version" || echo disabled)
mpich=$($bundle_mpi && echo "$mpich_version" || echo disabled)
EOF

if $plan_only; then
  echo "RESULT=PLAN: no download, build or installation was performed."
  exit 0
fi

required_commands=(tar make cc c++ realpath)
for required_command in "${required_commands[@]}"; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    echo "Missing prerequisite command: $required_command" >&2
    echo "A basic C/C++ toolchain, make and tar must be supplied by the host." >&2
    exit 1
  fi
done
if ! command -v sha256sum >/dev/null 2>&1 &&
   ! command -v shasum >/dev/null 2>&1; then
  echo "Missing SHA-256 tool: install sha256sum or shasum." >&2
  exit 1
fi
if ! $offline &&
   ! command -v curl >/dev/null 2>&1 &&
   ! command -v wget >/dev/null 2>&1; then
  echo "Missing downloader: install curl or wget." >&2
  exit 1
fi

mkdir -p -- "$prefix" "$build_dir" "$cache_dir"
state_dir="$prefix/share/aprl/bootstrap"
mkdir -p -- "$state_dir"
if ! $bundle_mpi && compgen -G "$state_dir/mpich-*.installed" >/dev/null; then
  echo "This prefix already contains a bootstrapped MPICH." >&2
  echo "Use --bundle-mpi again or choose a fresh prefix for the platform MPI." >&2
  echo "Mixing the retained private MPI runtime with a new site MPI is unsafe." >&2
  exit 1
fi

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

sha256_text() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

download_verified() {
  local file_name="$1"
  local url="$2"
  local expected="$3"
  local destination="$cache_dir/$file_name"
  local temporary="$destination.partial.$$"
  local actual=""

  if [[ -f "$destination" ]]; then
    actual="$(sha256_file "$destination")"
    if [[ "$actual" == "$expected" ]]; then
      echo "Using verified cache: $destination" >&2
      printf '%s\n' "$destination"
      return 0
    fi
    echo "Discarding cache entry with wrong SHA-256: $destination" >&2
    rm -f -- "$destination"
  fi
  if $offline; then
    echo "Offline mode: missing verified cache entry $destination" >&2
    exit 1
  fi

  echo "Downloading $url" >&2
  rm -f -- "$temporary"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --connect-timeout 20 -o "$temporary" "$url"
  else
    wget -O "$temporary" "$url"
  fi
  actual="$(sha256_file "$temporary")"
  if [[ "$actual" != "$expected" ]]; then
    rm -f -- "$temporary"
    echo "SHA-256 mismatch for $file_name" >&2
    echo "  expected: $expected" >&2
    echo "  actual:   $actual" >&2
    exit 1
  fi
  mv -- "$temporary" "$destination"
  printf '%s\n' "$destination"
}

remove_tree() {
  local target="$1"
  case "$target" in
    "$build_dir"|"$build_dir"/*) ;;
    *) echo "Refusing to remove path outside build root: $target" >&2; exit 1 ;;
  esac
  if [[ -x "$prefix/bin/cmake" ]]; then
    "$prefix/bin/cmake" -E remove_directory "$target"
  elif command -v cmake >/dev/null 2>&1; then
    cmake -E remove_directory "$target"
  else
    rm -rf -- "$target"
  fi
}

extract_source() {
  local archive="$1"
  local destination="$2"
  local entries=()
  local nested=""
  remove_tree "$destination"
  mkdir -p -- "$destination"
  tar -xzf "$archive" -C "$destination"

  # Release archives disagree on whether their member names begin with "./".
  # Flatten a single top-level source directory after extraction instead of
  # relying on a fragile --strip-components count.
  shopt -s dotglob nullglob
  entries=("$destination"/*)
  shopt -u dotglob nullglob
  if ((${#entries[@]} == 1)) && [[ -d "${entries[0]}" ]]; then
    nested="${destination}.nested.$$"
    mv -- "${entries[0]}" "$nested"
    rmdir -- "$destination"
    mv -- "$nested" "$destination"
  fi
}

version_at_least() {
  local actual="$1"
  local required="$2"
  [[ "$(printf '%s\n%s\n' "$required" "$actual" | sort -V | head -n 1)" == "$required" ]]
}

cmake_command=""
if [[ -x "$prefix/bin/cmake" ]]; then
  candidate_cmake="$prefix/bin/cmake"
  candidate_version="$("$candidate_cmake" --version | awk 'NR==1 {print $3}')"
  if version_at_least "$candidate_version" "3.16"; then
    cmake_command="$candidate_cmake"
  fi
elif command -v cmake >/dev/null 2>&1; then
  candidate_cmake="$(command -v cmake)"
  candidate_version="$("$candidate_cmake" --version | awk 'NR==1 {print $3}')"
  if version_at_least "$candidate_version" "3.16"; then
    cmake_command="$candidate_cmake"
  fi
fi
if [[ -z "$cmake_command" ]]; then
  if [[ -z "$cmake_platform" ]]; then
    echo "No CMake >= 3.16 was found and no prebuilt CMake is pinned for $host_arch." >&2
    echo "Supply CMake on PATH or set APRL_CMAKE_URL and APRL_CMAKE_SHA256." >&2
    exit 1
  fi
  cmake_archive="$(download_verified \
    "cmake-${cmake_version}-${cmake_platform}.tar.gz" \
    "$cmake_url" "$cmake_sha256")"
  cmake_marker="$state_dir/cmake-${cmake_version}.installed"
  if [[ ! -x "$prefix/bin/cmake" || ! -f "$cmake_marker" ]]; then
    echo "Installing private CMake $cmake_version"
    tar -xzf "$cmake_archive" --strip-components=1 -C "$prefix"
    : > "$cmake_marker"
  fi
  cmake_command="$prefix/bin/cmake"
fi
echo "cmake_command=$cmake_command"

resolve_command() {
  local request="$1"
  # Preserve the public MPI wrapper name instead of resolving a wrapper
  # symlink to an implementation detail such as Open MPI's opal_wrapper.
  # CMake and the dependency lock should both record the command the user
  # selected.
  command -v -- "$request" 2>/dev/null || return 1
}

derive_mpi_cc() {
  local cxx_path="$1"
  local wrapper_dir=""
  local candidate=""
  wrapper_dir="$(dirname -- "$cxx_path")"
  for candidate in mpicc mpiicc mpcc; do
    if [[ -x "$wrapper_dir/$candidate" ]]; then
      printf '%s\n' "$wrapper_dir/$candidate"
      return 0
    fi
  done
  command -v mpicc 2>/dev/null || return 1
}

if $bundle_mpi; then
  mpich_marker="$state_dir/mpich-${mpich_version}-ch3-sock.installed"
  if [[ ! -x "$prefix/bin/mpic++" || ! -x "$prefix/bin/mpicc" ||
        ! -f "$mpich_marker" ]]; then
    mpich_archive="$(download_verified "mpich-${mpich_version}.tar.gz" \
      "$mpich_url" "$mpich_sha256")"
    echo "Building private MPICH $mpich_version"
    mpich_source="$build_dir/mpich-source"
    mpich_build="$build_dir/mpich-build"
    extract_source "$mpich_archive" "$mpich_source"
    remove_tree "$mpich_build"
    mkdir -p -- "$mpich_build"
    (
      cd "$mpich_build"
      "$mpich_source/configure" \
        --prefix="$prefix" \
        --disable-fortran \
        --enable-shared \
        --disable-static \
        --with-device=ch3:sock \
        --with-hwloc=no \
        --with-wrapper-dl-type=none
      make -j "$jobs"
      make install
    )
    find "$state_dir" -maxdepth 1 -type f -name 'mpich-*.installed' \
      ! -name "$(basename -- "$mpich_marker")" -delete
    : > "$mpich_marker"
  fi
  mpi_cc="$prefix/bin/mpicc"
  mpi_cxx="$prefix/bin/mpic++"
  mpi_launcher="$prefix/bin/mpiexec"
else
  mpi_cxx="$(resolve_command "$mpi_cxx" || true)"
  if [[ -z "$mpi_cxx" ]]; then
    echo "Platform MPI C++ wrapper not found; load the site MPI module or use --bundle-mpi." >&2
    exit 1
  fi
  if [[ -n "$mpi_cc" ]]; then
    mpi_cc="$(resolve_command "$mpi_cc" || true)"
  else
    mpi_cc="$(derive_mpi_cc "$mpi_cxx" || true)"
  fi
  if [[ -z "$mpi_cc" ]]; then
    echo "Could not derive the MPI C wrapper next to $mpi_cxx; pass --mpi-cc." >&2
    exit 1
  fi
fi

mpi_show="$($mpi_cc -show 2>/dev/null || $mpi_cc --showme 2>/dev/null || true)"
mpi_identity="$(printf '%s\n%s\n' "$mpi_cc" "$mpi_show" | sha256_text)"
echo "mpi_cc=$mpi_cc"
echo "mpi_cxx=$mpi_cxx"
echo "mpi_launcher=$mpi_launcher"
echo "mpi_identity=$mpi_identity"

yaml_marker="$state_dir/yaml-cpp-${yaml_cpp_version}.installed"
if [[ ! -f "$yaml_marker" || ! -e "$prefix/lib/libyaml-cpp.so" ]]; then
  yaml_archive="$(download_verified "yaml-cpp-${yaml_cpp_version}.tar.gz" \
    "$yaml_cpp_url" "$yaml_cpp_sha256")"
  echo "Building yaml-cpp $yaml_cpp_version"
  yaml_source="$build_dir/yaml-cpp-source"
  yaml_build="$build_dir/yaml-cpp-build"
  extract_source "$yaml_archive" "$yaml_source"
  remove_tree "$yaml_build"
  "$cmake_command" -S "$yaml_source" -B "$yaml_build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_INSTALL_RPATH='\$ORIGIN' \
    -DBUILD_SHARED_LIBS=ON \
    -DYAML_BUILD_SHARED_LIBS=ON \
    -DYAML_CPP_BUILD_TESTS=OFF \
    -DYAML_CPP_BUILD_TOOLS=OFF \
    -DYAML_CPP_BUILD_CONTRIB=OFF
  "$cmake_command" --build "$yaml_build" --parallel "$jobs"
  "$cmake_command" --install "$yaml_build"
  : > "$yaml_marker"
fi

hdf5_marker="$state_dir/hdf5-${hdf5_version}-shared-${mpi_identity}.installed"
if [[ ! -f "$hdf5_marker" || ! -e "$prefix/lib/libhdf5.so" ||
      ! -x "$prefix/bin/h5pcc" ]]; then
  hdf5_archive="$(download_verified "hdf5-${hdf5_version}.tar.gz" \
    "$hdf5_url" "$hdf5_sha256")"
  echo "Building parallel HDF5 $hdf5_version"
  hdf5_source="$build_dir/hdf5-source"
  hdf5_build="$build_dir/hdf5-build"
  extract_source "$hdf5_archive" "$hdf5_source"
  remove_tree "$hdf5_build"
  "$cmake_command" -S "$hdf5_source" -B "$hdf5_build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_C_COMPILER="$mpi_cc" \
    -DCMAKE_INSTALL_RPATH='\$ORIGIN' \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF \
    -DONLY_SHARED_LIBS=ON \
    -DBUILD_TESTING=OFF \
    -DHDF5_ENABLE_PARALLEL=ON \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_FORTRAN=OFF \
    -DHDF5_BUILD_JAVA=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF \
    -DHDF5_ENABLE_SZIP_SUPPORT=OFF
  "$cmake_command" --build "$hdf5_build" --parallel "$jobs"
  "$cmake_command" --install "$hdf5_build"
  find "$state_dir" -maxdepth 1 -type f -name 'hdf5-*.installed' \
    ! -name "$(basename -- "$hdf5_marker")" -delete
  : > "$hdf5_marker"
fi

if [[ "$profile" == "full" ]]; then
  fftw_marker="$state_dir/fftw-${fftw_version}-portable.installed"
  if [[ ! -f "$fftw_marker" || ! -e "$prefix/lib/libfftw3.so" ||
        ! -e "$prefix/lib/libfftw3_threads.so" ]]; then
    fftw_archive="$(download_verified "fftw-${fftw_version}.tar.gz" \
      "$fftw_url" "$fftw_sha256")"
    echo "Building FFTW $fftw_version"
    fftw_source="$build_dir/fftw-source"
    fftw_build="$build_dir/fftw-build"
    extract_source "$fftw_archive" "$fftw_source"
    remove_tree "$fftw_build"
    mkdir -p -- "$fftw_build"
    (
      cd "$fftw_build"
      # FFTW otherwise adds -mtune=native automatically.  Keep the local
      # runtime usable on heterogeneous nodes of the same architecture.
      CFLAGS="${APRL_FFTW_CFLAGS:--O3 -fomit-frame-pointer}" \
        "$fftw_source/configure" \
        --prefix="$prefix" \
        --libdir="$prefix/lib" \
        --enable-shared \
        --disable-static \
        --enable-threads \
        --disable-fortran
      make -j "$jobs"
      make install
    )
    find "$state_dir" -maxdepth 1 -type f -name 'fftw-*.installed' \
      ! -name "$(basename -- "$fftw_marker")" -delete
    : > "$fftw_marker"
  fi
fi

runtime_env="$prefix/share/aprl/runtime-env.sh"
cat > "$runtime_env" <<'EOF'
#!/usr/bin/env bash
# Source this file before running an APRL installation whose site MPI module
# is not already active.  APRL binaries themselves also carry relative RPATHs.
_aprl_runtime_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
export PATH="${_aprl_runtime_root}/bin:${PATH}"
export LD_LIBRARY_PATH="${_aprl_runtime_root}/lib:${_aprl_runtime_root}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export CMAKE_PREFIX_PATH="${_aprl_runtime_root}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"
export PKG_CONFIG_PATH="${_aprl_runtime_root}/lib/pkgconfig:${_aprl_runtime_root}/lib64/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
unset _aprl_runtime_root
EOF
chmod +x "$runtime_env"

lock_file="$prefix/share/aprl/dependencies.lock"
cat > "$lock_file" <<EOF
generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
host_os=$host_os
host_arch=$host_arch
profile=$profile
prefix=$prefix
cmake_version=$cmake_version
cmake_url=$cmake_url
cmake_sha256=$cmake_sha256
yaml_cpp_version=$yaml_cpp_version
yaml_cpp_url=$yaml_cpp_url
yaml_cpp_sha256=$yaml_cpp_sha256
hdf5_version=$hdf5_version
hdf5_url=$hdf5_url
hdf5_sha256=$hdf5_sha256
fftw_version=$([[ "$profile" == "full" ]] && echo "$fftw_version" || echo disabled)
fftw_url=$fftw_url
fftw_sha256=$fftw_sha256
mpi_provider=$($bundle_mpi && echo bundled-mpich || echo platform)
mpi_cc=$mpi_cc
mpi_cxx=$mpi_cxx
mpi_launcher=$mpi_launcher
mpi_identity=$mpi_identity
mpich_version=$($bundle_mpi && echo "$mpich_version" || echo disabled)
mpich_url=$mpich_url
mpich_sha256=$mpich_sha256
EOF

if [[ -x "$prefix/bin/h5pcc" ]]; then
  echo "=== local parallel-HDF5 verification ==="
  "$prefix/bin/h5pcc" -showconfig 2>/dev/null |
    grep -E 'HDF5 Version|Parallel HDF5|C Compiler|Installation point' || true
fi

if ! $keep_build; then
  remove_tree "$build_dir"
fi

echo "RESULT=PASS: local dependencies are ready."
echo "  prefix: $prefix"
echo "  lock: $lock_file"
echo "  environment: source $runtime_env"
