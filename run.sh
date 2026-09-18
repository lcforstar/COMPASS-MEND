#!/usr/bin/env bash
# COMPASS 最简启动入口：./run.sh [站点名]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${ROOT}/scripts/run_compass_linux.sh" "${1:-changwu}"
