#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: ./check_dependencies.sh [options]

Compile and run a small dependency probe. No package is installed.

Options:
  --profile core|full       core: simulator; full: include FFTW tools (default)
  --report FILE             write the complete diagnostic report to FILE
  --mpi-cxx COMMAND         MPI C++ wrapper (default: $MPICXX or mpic++)
  --mpi-launcher COMMAND    launcher accepting '-n 2' (default: $MPIEXEC or mpiexec)
  --allow-serial-hdf5       allow a single-rank-only HDF5 installation
  --skip-mpi-run            skip the two-rank runtime/I/O check
  --sdds-root DIR           validate an existing official SDDS build
  --require-sdds            fail unless --sdds-root is supplied and valid
  --keep-probe              retain the temporary probe build directory
  -h, --help                show this help
EOF
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
profile="full"
report_file=""
mpi_cxx="${MPICXX:-mpic++}"
mpi_launcher="${MPIEXEC:-mpiexec}"
allow_serial_hdf5=false
skip_mpi_run=false
sdds_root="${SDDS_ROOT:-}"
require_sdds=false
keep_probe=false

while (($#)); do
  case "$1" in
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --report)
      [[ $# -ge 2 ]] || { echo "--report needs a value" >&2; exit 2; }
      report_file="$2"; shift 2 ;;
    --mpi-cxx)
      [[ $# -ge 2 ]] || { echo "--mpi-cxx needs a value" >&2; exit 2; }
      mpi_cxx="$2"; shift 2 ;;
    --mpi-launcher)
      [[ $# -ge 2 ]] || { echo "--mpi-launcher needs a value" >&2; exit 2; }
      mpi_launcher="$2"; shift 2 ;;
    --allow-serial-hdf5) allow_serial_hdf5=true; shift ;;
    --skip-mpi-run) skip_mpi_run=true; shift ;;
    --sdds-root)
      [[ $# -ge 2 ]] || { echo "--sdds-root needs a value" >&2; exit 2; }
      sdds_root="$2"; shift 2 ;;
    --require-sdds) require_sdds=true; shift ;;
    --keep-probe) keep_probe=true; shift ;;
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

if [[ -z "$report_file" ]]; then
  report_file="$(pwd -P)/dependency-report-$(date -u +%Y%m%dT%H%M%SZ).log"
elif [[ "$report_file" != /* ]]; then
  report_file="$(pwd -P)/$report_file"
fi
mkdir -p -- "$(dirname -- "$report_file")"
touch -- "$report_file"
exec > >(tee "$report_file") 2>&1

echo "=== dependency report ==="
echo "generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "source_root=$script_dir"
echo "profile=$profile"
echo "mpi_cxx_request=$mpi_cxx"
echo "mpi_launcher_request=$mpi_launcher"
echo "parallel_hdf5_required=$($allow_serial_hdf5 && echo no || echo yes)"
echo "sdds_requested=$([[ -n "$sdds_root" ]] && echo yes || echo no)"
echo "report=$report_file"
echo

echo "=== host ==="
uname -a || true
if [[ -r /etc/os-release ]]; then
  sed -n '1,40p' /etc/os-release
else
  echo "/etc/os-release is unavailable"
fi
echo

echo "=== administrator-facing requirements ==="
cat <<EOF
Required core interfaces:
  - CMake >= 3.16
  - C++11 compiler and standard build tools
  - one MPI implementation with matching compiler wrapper and launcher
  - HDF5 C >= 1.10; parallel HDF5 is required for multi-rank particle input
  - yaml-cpp >= 0.6
Full profile additionally requires:
  - pkg-config
  - FFTW3 double precision and fftw3_threads
Optional native Elegant conversion additionally requires:
  - a built official SDDS tree (SDDS1, rpnlib, mdbmth, mdblib)
  - zlib, liblzma and GSL development files

Do not mix MPI families. The MPI wrapper, parallel HDF5 library and runtime
launcher must come from the same compiler/MPI module stack. On a managed
cluster, ask the administrator for a consistent module collection rather
than combining system, Conda and user-built MPI/HDF5 libraries.

Common package requests (names vary by release/site):
  Debian/Ubuntu:
    build-essential cmake ninja-build pkg-config openmpi-bin libopenmpi-dev
    libhdf5-openmpi-dev libyaml-cpp-dev libfftw3-dev zlib1g-dev
    liblzma-dev libgsl-dev
  Fedora/RHEL/Rocky/Alma:
    gcc gcc-c++ make cmake ninja-build pkgconf-pkg-config openmpi-devel
    hdf5-openmpi-devel yaml-cpp-devel fftw-devel zlib-devel xz-devel gsl-devel
  openSUSE/SLES:
    gcc-c++ make cmake ninja pkg-config openmpi-devel hdf5-openmpi-devel
    yaml-cpp-devel fftw3-devel zlib-devel liblzma-devel gsl-devel
EOF
echo

missing=0
required_commands=(cmake "$mpi_cxx")
if [[ "$profile" == "full" ]]; then
  required_commands+=(pkg-config)
fi
if ! $skip_mpi_run && ! $allow_serial_hdf5; then
  required_commands+=("$mpi_launcher")
fi

echo "=== command inventory ==="
for command_name in "${required_commands[@]}"; do
  if command_path="$(command -v -- "$command_name" 2>/dev/null)"; then
    echo "$command_name=$command_path"
  else
    echo "MISSING: $command_name"
    missing=1
  fi
done
for optional_command in h5pcc h5pcc.openmpi h5cc ninja ldd; do
  if command_path="$(command -v -- "$optional_command" 2>/dev/null)"; then
    echo "$optional_command=$command_path"
  fi
done
echo

if ((missing)); then
  echo "RESULT=FAIL: required commands are missing."
  echo "Give this report to the system administrator on an offline cluster."
  exit 1
fi

echo "=== tool versions and wrappers ==="
cmake --version | sed -n '1,3p'
"$mpi_cxx" --version | sed -n '1,3p' || true
"$mpi_cxx" --showme:command 2>/dev/null ||
  "$mpi_cxx" -show 2>/dev/null || true
if command -v h5pcc >/dev/null 2>&1; then
  h5pcc -showconfig 2>/dev/null |
    grep -E 'HDF5 Version|Parallel HDF5|C Compiler|Installation point' || true
elif command -v h5pcc.openmpi >/dev/null 2>&1; then
  h5pcc.openmpi -showconfig 2>/dev/null |
    grep -E 'HDF5 Version|Parallel HDF5|C Compiler|Installation point' || true
elif command -v h5cc >/dev/null 2>&1; then
  h5cc -showconfig 2>/dev/null |
    grep -E 'HDF5 Version|Parallel HDF5|C Compiler|Installation point' || true
fi
echo

probe_parent="${TMPDIR:-/tmp}"
probe_dir="$(mktemp -d "$probe_parent/fel-dependency-check.XXXXXX")"
cleanup_probe() {
  if $keep_probe; then
    echo "probe_directory_retained=$probe_dir"
  elif [[ -n "$probe_dir" && "$probe_dir" == "$probe_parent"/fel-dependency-check.* ]]; then
    rm -rf -- "$probe_dir"
  fi
}
trap cleanup_probe EXIT

probe_full=OFF
[[ "$profile" == "full" ]] && probe_full=ON
probe_parallel=ON
$allow_serial_hdf5 && probe_parallel=OFF
probe_sdds=OFF
[[ -n "$sdds_root" ]] && probe_sdds=ON

configure_command=(
  cmake
  -S "$script_dir/maintenance/dependency_probe"
  -B "$probe_dir/build"
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_CXX_COMPILER="$mpi_cxx"
  -DFEL_PROBE_FULL="$probe_full"
  -DFEL_PROBE_REQUIRE_PARALLEL_HDF5="$probe_parallel"
  -DFEL_PROBE_SDDS="$probe_sdds"
)
if [[ -n "$sdds_root" ]]; then
  configure_command+=("-DSDDS_ROOT=$sdds_root")
fi

echo "=== CMake discovery ==="
if ! "${configure_command[@]}"; then
  echo "RESULT=FAIL: CMake could not resolve a consistent dependency set."
  echo "Inspect the selected include/library paths above; do not bypass this by mixing MPI stacks."
  exit 1
fi
echo

echo "=== compile and link probe ==="
if ! cmake --build "$probe_dir/build" --parallel 2; then
  echo "RESULT=FAIL: the discovered headers and libraries do not compile/link together."
  exit 1
fi
probe_executable="$probe_dir/build/fel_dependency_probe"
if command -v ldd >/dev/null 2>&1; then
  echo "--- resolved shared libraries ---"
  ldd "$probe_executable" || true
fi
echo

echo "=== single-rank runtime probe ==="
if ! "$probe_executable" "$probe_dir/single-rank.h5"; then
  echo "RESULT=FAIL: single-rank dependency runtime probe failed."
  exit 1
fi
echo

if $allow_serial_hdf5; then
  echo "=== two-rank runtime probe ==="
  echo "SKIPPED: serial-HDF5 mode is explicitly limited to single-rank runs."
elif $skip_mpi_run; then
  echo "=== two-rank runtime probe ==="
  echo "SKIPPED by request. Run this checker inside an allocation before production use."
else
  echo "=== two-rank collective MPI/HDF5 runtime probe ==="
  if ! "$mpi_launcher" -n 2 "$probe_executable" "$probe_dir/two-rank.h5"; then
    echo "RESULT=FAIL: MPI launch or collective parallel-HDF5 I/O failed."
    echo "If the login node forbids MPI launch, rerun with --skip-mpi-run and then run the checker inside a scheduler allocation."
    exit 1
  fi
fi
echo

if [[ -z "$sdds_root" ]]; then
  echo "=== optional official SDDS converter ==="
  echo "NOT CHECKED: supply --sdds-root /path/to/built/SDDS."
  echo "The simulator and HDF5-v3 text converter remain available without it."
  if $require_sdds; then
    echo "RESULT=FAIL: SDDS was required."
    exit 1
  fi
  echo
fi

echo "RESULT=PASS: compile, link and requested runtime dependency probes succeeded."
