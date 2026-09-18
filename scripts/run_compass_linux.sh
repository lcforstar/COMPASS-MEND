#!/usr/bin/env bash
# COMPASS Linux 后台启动器
# 用法：./scripts/run_compass_linux.sh [站点名]
# 环境变量：COMPASS_LOG_WINDOW=0 可关闭自动弹出的日志窗口
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATLAB_BIN="${MATLAB_BIN:-/home/zjzx/MATLAB/R2025a/bin/matlab}"
SITE="${1:-changwu}"
SITE_DIR="${PROJECT_ROOT}/data/sites/${SITE}"

if [[ ! -d "${SITE_DIR}" || ! -f "${SITE_DIR}/site.toml" ]]; then
    echo "找不到站点配置：${SITE_DIR}/site.toml" >&2
    echo "可用站点：$(find "${PROJECT_ROOT}/data/sites" -mindepth 1 -maxdepth 1 -type f -name site.toml -printf '%h\n' | xargs -r -n1 basename | sort | tr '\n' ' ')" >&2
    exit 1
fi

# 必须读取 site.toml 的 [output].output_dir；结果目录由用户在 data 中
# 明确配置并提前创建，启动器不猜测、不自动创建结果目录。
configured_output="$(sed -n 's/^[[:space:]]*output_dir[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "${SITE_DIR}/site.toml" | head -1)"
if [[ -z "${configured_output}" ]]; then
    echo "site.toml 未配置 [output].output_dir：${SITE_DIR}/site.toml" >&2
    echo "请在 [output] 中填写已存在的输出目录后再运行。" >&2
    exit 1
elif [[ "${configured_output}" = /* ]]; then
    RESULT_DIR="${configured_output}"
else
    RESULT_DIR="$(realpath -m "${SITE_DIR}/${configured_output}")"
fi

if [[ ! -d "${RESULT_DIR}" ]]; then
    echo "输出目录不存在：${RESULT_DIR}" >&2
    echo "请编辑 ${SITE_DIR}/site.toml 的 [output].output_dir，填写已存在的目录。" >&2
    exit 1
fi

LOG_DIR="${COMPASS_LOG_DIR:-${RESULT_DIR}/logs}"

if [[ ! -x "${MATLAB_BIN}" ]]; then
    echo "找不到可执行的 MATLAB：${MATLAB_BIN}" >&2
    exit 1
fi

mkdir -p "${LOG_DIR}"
timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="${LOG_DIR}/compass_${SITE}_${timestamp}.log"

matlab_code="cd('${PROJECT_ROOT}'); run_site('${SITE}');"
# R2025a 在本机 Ubuntu 26.04 上使用 -batch 会触发 ddux_matlab 段错误；
# -r + 显式 exit 可避开该路径，同时保留无桌面后台运行。
nohup "${MATLAB_BIN}" -nojvm -nodisplay -nosplash -r "${matlab_code} exit" >"${log_file}" 2>&1 < /dev/null &
pid=$!

# 给 MATLAB 几秒完成进程级启动；若立即崩溃，直接报错并展示日志尾部。
sleep 3
if ! kill -0 "${pid}" 2>/dev/null; then
    echo "COMPASS 启动失败（MATLAB PID ${pid} 已退出）" >&2
    echo "日志：${log_file}" >&2
    tail -30 "${log_file}" >&2 || true
    exit 1
fi

echo "COMPASS 已在 Linux 后台启动"
echo "PID: ${pid}"
echo "站点: ${SITE}"
echo "结果目录: ${RESULT_DIR}"
echo "日志: ${log_file}"
echo "查看进度: tail -f '${log_file}'"

# 在桌面环境中打开一个详细实时进度窗口；服务器/无图形环境自动跳过。
if [[ "${COMPASS_LOG_WINDOW:-1}" != "0" ]]; then
    gui_script="${PROJECT_ROOT}/scripts/compass_monitor_gui.py"
    if command -v python3 >/dev/null 2>&1 && [[ -f "${gui_script}" ]] && [[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]]; then
        setsid nohup python3 "${gui_script}" "${SITE}" "${log_file}" "${pid}" >/dev/null 2>&1 < /dev/null &
        exit 0
    fi
    monitor_script="${PROJECT_ROOT}/scripts/compass_monitor.sh"
    viewer_cmd="$(printf '%q' "${monitor_script}") $(printf '%q' "${SITE}") $(printf '%q' "${log_file}") $(printf '%q' "${pid}")"
    if command -v gnome-terminal >/dev/null 2>&1; then
        gnome-terminal --title="COMPASS Land System Model - ${SITE} 实时进度" -- bash -lc "${viewer_cmd}" >/dev/null 2>&1 &
    elif command -v konsole >/dev/null 2>&1; then
        konsole --new-tab -p tabtitle="COMPASS Land System Model - ${SITE} 实时进度" -e bash -lc "${viewer_cmd}" >/dev/null 2>&1 &
    elif command -v ptyxis >/dev/null 2>&1; then
        # Ubuntu 近期版本的 x-terminal-emulator 通常指向 Ptyxis，使用其专用参数。
        ptyxis -s -T "COMPASS Land System Model - ${SITE} 实时进度" -- bash -lc "${viewer_cmd}" >/dev/null 2>&1 &
    elif command -v x-terminal-emulator >/dev/null 2>&1; then
        x-terminal-emulator -T "COMPASS Land System Model - ${SITE} 实时进度" -e bash -lc "${viewer_cmd}" >/dev/null 2>&1 &
    else
        echo "未找到终端程序，未自动打开日志窗口；可手动执行 tail -f 查看。"
    fi
fi
