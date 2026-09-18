#!/usr/bin/env bash
# COMPASS 进度 GUI：./monitor.sh [站点] [日志文件] [PID]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SITE="${1:-changwu}"
LOG="${2:-}"
PID="${3:-}"
GUI="${ROOT}/scripts/compass_monitor_gui.py"
if [[ ! -f "${GUI}" ]]; then
    echo "找不到 GUI: ${GUI}" >&2
    exit 1
fi
if command -v python3 >/dev/null 2>&1; then
    PY=python3
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    echo "未找到 python3" >&2
    exit 1
fi
if [[ -n "${LOG}" && -n "${PID}" ]]; then
    exec "${PY}" "${GUI}" "${SITE}" "${LOG}" "${PID}"
elif [[ -n "${LOG}" ]]; then
    exec "${PY}" "${GUI}" "${SITE}" "${LOG}"
else
    exec "${PY}" "${GUI}" "${SITE}"
fi
