#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage: sudo ./install_dependencies.sh [options]

Install packaged build dependencies. This script must run as root and never
downloads or deletes the project source. Official SDDS itself is not installed.

Options:
  --profile core|full   install core or all post-processing dependencies
  --skip-refresh        do not refresh package-manager metadata
  --report FILE         package installation log
  -h, --help            show this help
EOF
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
profile="full"
skip_refresh=false
report_file=""

while (($#)); do
  case "$1" in
    --profile)
      [[ $# -ge 2 ]] || { echo "--profile needs a value" >&2; exit 2; }
      profile="$2"; shift 2 ;;
    --skip-refresh) skip_refresh=true; shift ;;
    --report)
      [[ $# -ge 2 ]] || { echo "--report needs a value" >&2; exit 2; }
      report_file="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$profile" != "core" && "$profile" != "full" ]]; then
  echo "--profile must be core or full" >&2
  exit 2
fi
if [[ -z "$report_file" ]]; then
  report_file="$(pwd -P)/dependency-install-$(date -u +%Y%m%dT%H%M%SZ).log"
elif [[ "$report_file" != /* ]]; then
  report_file="$(pwd -P)/$report_file"
fi
mkdir -p -- "$(dirname -- "$report_file")"
touch -- "$report_file"
exec > >(tee "$report_file") 2>&1

echo "=== dependency package installation ==="
echo "generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "profile=$profile"
echo "report=$report_file"

if ((EUID != 0)); then
  echo "RESULT=FAIL: administrator privileges are required."
  echo "No package operation was attempted. Give this log and the output of"
  echo "./check_dependencies.sh to the cluster administrator."
  exit 1
fi

manager=""
for candidate in apt-get dnf yum zypper pacman; do
  if command -v "$candidate" >/dev/null 2>&1; then
    manager="$candidate"
    break
  fi
done
if [[ -z "$manager" ]]; then
  echo "RESULT=FAIL: no supported package manager was found."
  echo "Use check_dependencies.sh to generate the manual requirement report."
  exit 1
fi

packages=()
case "$manager" in
  apt-get)
    packages=(build-essential cmake ninja-build pkg-config openmpi-bin
      libopenmpi-dev libhdf5-openmpi-dev libyaml-cpp-dev zlib1g-dev
      liblzma-dev libgsl-dev)
    [[ "$profile" == "full" ]] && packages+=(libfftw3-dev)
    if ! $skip_refresh; then
      apt-get update
    fi
    DEBIAN_FRONTEND=noninteractive apt-get install -y "${packages[@]}"
    ;;
  dnf|yum)
    packages=(gcc gcc-c++ make cmake ninja-build pkgconf-pkg-config
      openmpi openmpi-devel hdf5-openmpi hdf5-openmpi-devel yaml-cpp-devel
      zlib-devel xz-devel gsl-devel)
    [[ "$profile" == "full" ]] && packages+=(fftw-devel)
    "$manager" install -y "${packages[@]}"
    ;;
  zypper)
    packages=(gcc-c++ make cmake ninja pkg-config openmpi-devel
      hdf5-openmpi-devel yaml-cpp-devel zlib-devel liblzma-devel gsl-devel)
    [[ "$profile" == "full" ]] && packages+=(fftw3-devel)
    if ! $skip_refresh; then
      zypper --non-interactive refresh
    fi
    zypper --non-interactive install "${packages[@]}"
    ;;
  pacman)
    packages=(base-devel cmake ninja pkgconf openmpi hdf5-openmpi yaml-cpp
      zlib xz gsl)
    [[ "$profile" == "full" ]] && packages+=(fftw)
    refresh_flag=""
    $skip_refresh || refresh_flag="-Sy"
    if [[ -n "$refresh_flag" ]]; then
      pacman "$refresh_flag" --needed --noconfirm "${packages[@]}"
    else
      pacman -S --needed --noconfirm "${packages[@]}"
    fi
    ;;
esac

echo
echo "Packaged dependencies installed. Official SDDS is deliberately not"
echo "downloaded: build the official tree separately and pass SDDS_ROOT."
verification_report="${report_file%.log}-verification.log"
echo "Running a root-safe compile/link check; the two-rank launch is skipped."
"$script_dir/check_dependencies.sh" --profile "$profile" \
  --skip-mpi-run --report "$verification_report"
echo "RESULT=PASS: package installation and compile/link verification succeeded."
echo "Before production, rerun check_dependencies.sh as the ordinary job user"
echo "inside an allocation so the two-rank MPI/HDF5 runtime check is exercised."

