#!/usr/bin/env bash
# Based on the COMPASS launcher by Jizixing Zhu.
# Change: 2026-10-31 Wenjie Liu -- portable independent-copy startup.
# Usage: run_compass_linux.sh [site-name-or-site-directory]
# Env: MATLAB_BIN, COMPASS_OUTPUT, COMPASS_LOG_DIR, COMPASS_LOG_WINDOW=0,
#      COMPASS_LAUNCH_WAIT=1. This script never changes scientific inputs.
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
SITE="${1:-changwu}"
if [[ -f "${SITE}/site.toml" ]]; then
    SITE_DIR="$(cd "${SITE}" && pwd -P)"
else
    SITE_DIR="${PROJECT_ROOT}/data/sites/${SITE}"
fi
if [[ ! -f "${SITE_DIR}/site.toml" ]]; then
    printf 'Site configuration not found: %s/site.toml\n' "${SITE_DIR}" >&2
    exit 1
fi
SITE_DIR="$(cd "${SITE_DIR}" && pwd -P)"
if [[ -n "${MATLAB_BIN:-}" ]]; then
    matlab_executable="$(command -v "${MATLAB_BIN}" || true)"
else
    matlab_executable="$(command -v matlab || true)"
fi
if [[ -z "${matlab_executable}" || ! -x "${matlab_executable}" ]]; then
    printf 'MATLAB not found. Set MATLAB_BIN or add matlab to PATH.\n' >&2
    exit 1
fi
section=''
configured_output=''
section_re='^[[:space:]]*\[([^]]+)\][[:space:]]*(#.*)?$'
double_re='^[[:space:]]*output_dir[[:space:]]*=[[:space:]]*"([^"]*)"[[:space:]]*(#.*)?$'
single_re="^[[:space:]]*output_dir[[:space:]]*=[[:space:]]*'([^']*)'[[:space:]]*(#.*)?$"
while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%$'\r'}"
    if [[ "${line}" =~ ${section_re} ]]; then section="${BASH_REMATCH[1]}"; continue; fi
    if [[ "${section}" == 'output' ]]; then
        if [[ "${line}" =~ ${double_re} || "${line}" =~ ${single_re} ]]; then
            configured_output="${BASH_REMATCH[1]}"
            break
        fi
    fi
done < "${SITE_DIR}/site.toml"
if [[ -z "${configured_output}" ]]; then
    printf 'site.toml must provide a quoted [output].output_dir.\n' >&2
    exit 1
fi
output_source='site.toml:output.output_dir'
output_base="${SITE_DIR}"
if [[ -n "${COMPASS_OUTPUT:-}" ]]; then
    configured_output="${COMPASS_OUTPUT}"; output_base="${PROJECT_ROOT}"; output_source='env:COMPASS_OUTPUT'
elif [[ -n "${STEMMUS_OUTPUT:-}" ]]; then
    configured_output="${STEMMUS_OUTPUT}"; output_base="${PROJECT_ROOT}"; output_source='env:STEMMUS_OUTPUT'
fi
if [[ "${configured_output}" = /* ]]; then RESULT_DIR="${configured_output}"; else RESULT_DIR="${output_base}/${configured_output}"; fi
mkdir -p "${RESULT_DIR}"
RESULT_DIR="$(cd "${RESULT_DIR}" && pwd -P)"
LOG_DIR="${COMPASS_LOG_DIR:-${RESULT_DIR}/logs}"
if [[ "${LOG_DIR}" != /* ]]; then LOG_DIR="${PROJECT_ROOT}/${LOG_DIR}"; fi
mkdir -p "${LOG_DIR}"
LOG_DIR="$(cd "${LOG_DIR}" && pwd -P)"
site_label="$(basename "${SITE_DIR}")"
site_label="$(printf '%s' "${site_label}" | tr -cs '[:alnum:]_.-' '_')"
run_id="${site_label}_$(date -u +%Y%m%d_%H%M%S)_$$_${RANDOM}"
log_file="${LOG_DIR}/compass_${run_id}.log"
run_receipt="${log_file%.log}.run.json"
export COMPASS_LAUNCH_SITE="${SITE_DIR}"
export COMPASS_RUN_ID="${run_id}"
export COMPASS_RUN_RECEIPT="${run_receipt}"
export COMPASS_LAUNCH_LOG="${log_file}"
export COMPASS_LAUNCH_OUTPUT_SOURCE="${output_source}"
export COMPASS_OUTPUT="${RESULT_DIR}"
if [[ -z "${COMPASS_DIAGNOSTIC_PROFILE:-}${STEMMUS_DIAGNOSTIC_PROFILE:-}" ]]; then export COMPASS_DIAGNOSTIC_PROFILE=full; fi
# Keep the JVM: the unchanged coupled residual ledger uses Java SHA-256.
nohup "${matlab_executable}" -sd "${PROJECT_ROOT}" -batch "addpath('scripts'); compass_batch_entry;" >"${log_file}" 2>&1 < /dev/null &
launched_pid=$!
printf 'COMPASS-MEND launched in background.\nLaunched PID: %s\nLog: %s\nReceipt: %s\n' "${launched_pid}" "${log_file}" "${run_receipt}"
if [[ "${COMPASS_LAUNCH_WAIT:-0}" == '1' ]]; then
    set +e
    wait "${launched_pid}"
    exit_code=$?
    set -e
    if [[ "${exit_code}" -ne 0 ]]; then printf 'MATLAB failed; see %s\n' "${log_file}" >&2; exit "${exit_code}"; fi
    if [[ ! -f "${run_receipt}" ]] || ! grep -Eq '"status"[[:space:]]*:[[:space:]]*"completed"' "${run_receipt}"; then
        printf 'MATLAB returned without a completed run receipt: %s\n' "${run_receipt}" >&2
        exit 1
    fi
    exit 0
fi
sleep 2
if ! kill -0 "${launched_pid}" 2>/dev/null; then
    set +e
    wait "${launched_pid}"
    exit_code=$?
    set -e
    if [[ "${exit_code}" -ne 0 ]] || [[ ! -f "${run_receipt}" ]] || ! grep -Eq '"status"[[:space:]]*:[[:space:]]*"completed"' "${run_receipt}"; then
        printf 'MATLAB failed during startup; see %s\n' "${log_file}" >&2
        tail -30 "${log_file}" >&2 || true
        exit 1
    fi
    exit 0
fi
if [[ "${COMPASS_LOG_WINDOW:-1}" != '0' && -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" && "${SITE_DIR}" == "${PROJECT_ROOT}/data/sites/${site_label}" ]] && command -v python3 >/dev/null 2>&1; then
    gui_script="${PROJECT_ROOT}/scripts/compass_monitor_gui.py"
    if [[ -f "${gui_script}" ]]; then
        nohup python3 "${gui_script}" "${site_label}" "${log_file}" "${launched_pid}" >/dev/null 2>&1 < /dev/null &
    fi
fi
