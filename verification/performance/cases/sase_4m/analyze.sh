#!/usr/bin/env bash
set -Eeuo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export CAMPAIGN_DIR
set -a
# shellcheck source=/dev/null
source "${CAMPAIGN_DIR}/campaign.env"
set +a

PYTHON="${CAMPAIGN_DIR}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]] || ! "${PYTHON}" -c 'import h5py, numpy, matplotlib' >/dev/null 2>&1; then
  echo "ERROR: 请先运行 ./prepare_python_env.sh。" >&2
  exit 1
fi
export MPLCONFIGDIR="${CAMPAIGN_DIR}/.venv/.matplotlib"
mkdir -p "${MPLCONFIGDIR}"
exec "${PYTHON}" "${CAMPAIGN_DIR}/scripts/analyze_campaign.py"
