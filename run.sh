#!/usr/bin/env bash
# Portable COMPASS-MEND launcher; invoke with bash run.sh [site].
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${ROOT}/scripts/run_compass_linux.sh" "${1:-changwu}"
