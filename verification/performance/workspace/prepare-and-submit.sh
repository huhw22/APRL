#!/usr/bin/env bash
# Install in <workspace>/tests/performance; no source-tree outputs are needed.
set -Eeuo pipefail
workspace="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
source "$workspace/cluster.env"
export PERFORMANCE_OUTPUT_DIR="${PERFORMANCE_OUTPUT_DIR:-$workspace/results/performance/run-$(date -u +%Y%m%dT%H%M%SZ)-$$}"
exec bash "$APRL_SOURCE_DIR/verification/performance/submit_cluster_test.sh" \
  "$@" "$workspace/source/aprl-baseline" "$APRL_DEPENDENCY_PREFIX"
