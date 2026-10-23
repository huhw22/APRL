#!/usr/bin/env bash
set -Eeuo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
VENV_DIR="${CAMPAIGN_DIR}/.venv"
export MPLCONFIGDIR="${VENV_DIR}/.matplotlib"
mkdir -p "${MPLCONFIGDIR}"

command -v python3 >/dev/null 2>&1 || {
  echo "ERROR: 找不到 python3。" >&2
  exit 1
}
if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
  python3 -m venv --system-site-packages "${VENV_DIR}"
fi
if ! "${VENV_DIR}/bin/python" -c 'import numpy, h5py, matplotlib' >/dev/null 2>&1; then
  "${VENV_DIR}/bin/python" -m pip install -r "${CAMPAIGN_DIR}/requirements.txt"
fi
"${VENV_DIR}/bin/python" - <<'PY'
import h5py
import matplotlib
import numpy
print("Python environment ready:",
      "numpy=" + numpy.__version__,
      "h5py=" + h5py.__version__,
      "matplotlib=" + matplotlib.__version__)
PY
