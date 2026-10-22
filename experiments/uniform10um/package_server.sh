#!/usr/bin/env bash
# Package source and experiment recipes only; never include generated particles.
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
default_name="aprl-uniform10um-source-$(date +%Y%m%d-%H%M%S).tar.gz"
output="${1:-${repo_root}/dist/${default_name}}"
if [[ "${output}" != /* ]]; then
  output="${PWD}/${output}"
fi
if [[ -e "${output}" ]]; then
  echo "Refusing to overwrite existing archive: ${output}" >&2
  exit 2
fi
mkdir -p -- "$(dirname -- "${output}")"

tar -C "${repo_root}" -czf "${output}" \
  --transform='s|^|aprl-uniform10um/|' \
  --exclude='*.h5' --exclude='*.hdf5' --exclude='*.log' \
  --exclude='*.pyc' --exclude='__pycache__' \
  CMakeLists.txt LICENSE README.md README.zh-CN.md \
  bootstrap_dependencies.sh check_dependencies.sh \
  install.sh install_dependencies.sh install_online.sh uninstall.sh \
  src utilities verification postprocess config examples docs \
  experiments/uniform10um

echo "Source package: ${output}"
du -h -- "${output}"
sha256sum -- "${output}"
echo "The HDF5 bunch is intentionally absent; regenerate it on the server."
