#!/usr/bin/env bash
# COMPASS 实时进度监控器（Linux 版）
set -u

SITE="${1:-changwu}"
LOG_FILE="${2:-}"
PID="${3:-}"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_DIR="${PROJECT_ROOT}/data/sites/${SITE}"
CFG="${SITE_DIR}/site.toml"
REFRESH_SEC="${COMPASS_MONITOR_REFRESH:-10}"

cfg_value() {
    local key="$1"
    sed -n "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*\"\([^\"]*\)\".*/\1/p" "${CFG}" 2>/dev/null | head -1
}

cfg_number() {
    local key="$1"
    sed -n "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*\([^ #]*\).*/\1/p" "${CFG}" 2>/dev/null | head -1
}

SITE_NAME="$(cfg_value name)"; SITE_NAME="${SITE_NAME:-${SITE}}"
LAT="$(cfg_number latitude)"; LON="$(cfg_number longitude)"; ELEV="$(cfg_number elevation)"
START="$(cfg_value start)"; END="$(cfg_value end)"; DT="$(cfg_number dt)"
OUT_CFG="$(cfg_value output_dir)"
if [[ -n "${OUT_CFG}" ]]; then
    if [[ "${OUT_CFG}" = /* ]]; then RESULT_DIR="${OUT_CFG}"; else RESULT_DIR="$(realpath -m "${SITE_DIR}/${OUT_CFG}")"; fi
else
    RESULT_DIR="${PROJECT_ROOT}/output/${SITE}"
fi

if [[ -z "${LOG_FILE}" ]]; then
    LOG_FILE="$(find "${RESULT_DIR}/logs" -maxdepth 1 -type f -name "compass_${SITE}_*.log" -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk 'NR==1{$1=""; sub(/^ /,""); print}')"
fi

TOTAL=""
if [[ -n "${START}" && -n "${END}" && -n "${DT}" ]]; then
    start_epoch="$(date -d "${START}" +%s 2>/dev/null || true)"
    end_epoch="$(date -d "${END}" +%s 2>/dev/null || true)"
    if [[ -n "${start_epoch}" && -n "${end_epoch}" && "${DT}" =~ ^[0-9]+$ && "${DT}" -gt 0 ]]; then
        TOTAL="$(( (end_epoch - start_epoch) / DT ))"
    fi
fi
TOTAL="${TOTAL:-8760}"

while true; do
    clear
    echo "=========================================================================="
    echo "  COMPASS Land System Model - Full-Year Run Progress"
    echo "=========================================================================="
    echo "  Project   : ${PROJECT_ROOT}"
    echo "  Site      : ${SITE_NAME} (${SITE})"
    echo "  Location  : ${LAT:-?} N, ${LON:-?} E, ${ELEV:-?} m a.s.l."
    echo "  Period    : ${START:-?}  ->  ${END:-?}   (dt ${DT:-?} s)"
    echo "  Output    : ${RESULT_DIR}"
    echo "  Log       : ${LOG_FILE:-not-created}"
    if [[ -n "${PID}" ]]; then
        if kill -0 "${PID}" 2>/dev/null; then
            echo "  Process   : RUNNING (PID ${PID})"
        else
            echo "  Process   : EXITED (PID ${PID})"
        fi
    fi
    echo "  Updated   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "--------------------------------------------------------------------------"

    latest_kt=0
    if [[ -n "${LOG_FILE}" && -f "${LOG_FILE}" ]]; then
        latest_line="$(tail -300 "${LOG_FILE}" 2>/dev/null | grep -Eo 'KT=[0-9]+/[0-9]+' | tail -1 || true)"
        if [[ "${latest_line}" =~ KT=([0-9]+)/([0-9]+) ]]; then
            latest_kt="${BASH_REMATCH[1]}"; TOTAL="${BASH_REMATCH[2]}"
        fi
    fi
    ckpt_kt="$(find "${RESULT_DIR}/checkpoints" -maxdepth 1 -type f -name 'checkpoint_KT*.mat' -printf '%f\n' 2>/dev/null | sed -n 's/.*KT\([0-9][0-9]*\)\.mat/\1/p' | sort -n | tail -1)"
    # checkpoint 文件名通常是 KT00168；10# 防止 Bash 把前导 0 当成八进制。
    if [[ -n "${ckpt_kt}" ]]; then ckpt_kt=$((10#${ckpt_kt})); fi
    if [[ -n "${ckpt_kt}" && "${ckpt_kt}" -gt "${latest_kt}" ]]; then latest_kt="${ckpt_kt}"; fi
    pct="$(awk -v k="${latest_kt}" -v t="${TOTAL}" 'BEGIN { if(t>0) printf "%.2f",100*k/t; else print "0.00" }')"
    filled="$(awk -v p="${pct}" 'BEGIN { n=int(p/2); if(n<0)n=0;if(n>50)n=50;print n }')"
    bar="$(printf '%*s' "${filled}" '' | tr ' ' '#')$(printf '%*s' "$((50-filled))" '' | tr ' ' '.')"
    echo "  Progress  : KT=${latest_kt}/${TOTAL} (${pct}%)"
    echo "  [${bar}]"

    if [[ -n "${LOG_FILE}" && -f "${LOG_FILE}" ]]; then
        log_bytes="$(stat -c %s "${LOG_FILE}" 2>/dev/null || echo 0)"
        log_time="$(stat -c '%y' "${LOG_FILE}" 2>/dev/null | cut -d. -f1 || echo '?')"
        echo "  Log size  : ${log_bytes} bytes; last write: ${log_time}"
        latest_status="$(tail -300 "${LOG_FILE}" | grep -E '\[COMPASS\] KT=' | tail -1 || true)"
        [[ -n "${latest_status}" ]] && echo "  Latest    : ${latest_status}"
        if [[ "${latest_kt}" -ge "${TOTAL}" && "${TOTAL}" -gt 0 ]]; then
            echo "  Speed     : completed; ETA: 0.0 min"
        elif [[ "${latest_kt}" -gt 0 ]]; then
            # 用最近两个检查点的写盘时间计算速度，避免日志创建时间、恢复运行
            # 或长时间缓冲造成虚假的超高速。
            mapfile -t recent_ckpts < <(find "${RESULT_DIR}/checkpoints" -maxdepth 1 -type f -name 'checkpoint_KT*.mat' -printf '%T@ %f\n' 2>/dev/null | sort -n | tail -2)
            rate=""
            if [[ "${#recent_ckpts[@]}" -eq 2 ]]; then
                read -r t1 f1 <<< "${recent_ckpts[0]}"
                read -r t2 f2 <<< "${recent_ckpts[1]}"
                k1="$(sed -n 's/.*KT\([0-9][0-9]*\)\.mat/\1/p' <<< "${f1}")"
                k2="$(sed -n 's/.*KT\([0-9][0-9]*\)\.mat/\1/p' <<< "${f2}")"
                if [[ -n "${k1}" && -n "${k2}" ]]; then
                    k1=$((10#${k1})); k2=$((10#${k2}))
                    delta_seconds="$(awk -v a="${t1}" -v b="${t2}" 'BEGIN {print b-a}')"
                    rate="$(awk -v dk="$((k2-k1))" -v dt="${delta_seconds}" 'BEGIN {if(dt>0) printf "%.2f",dk*60/dt; else print ""}')"
                fi
            fi
            if [[ -n "${rate}" && "${rate}" != "0.00" ]]; then
                remain=$(( TOTAL - latest_kt )); [[ "${remain}" -lt 0 ]] && remain=0
                eta_min="$(awk -v r="${remain}" -v rate="${rate}" 'BEGIN {if(rate>0) printf "%.1f",r/rate; else print "?"}')"
                echo "  Speed     : ${rate} KT/min; ETA: ${eta_min} min"
            else
                echo "  Speed     : calculating from checkpoints; ETA: pending"
            fi
        fi
        nonconv="$(find "${RESULT_DIR}/checkpoints" -maxdepth 1 -type f -name 'nonconv_diag_KT*.mat' 2>/dev/null | wc -l)"
        rollback="$(find "${RESULT_DIR}/checkpoints" -maxdepth 1 -type f -name 'resume_before_KT*.mat' 2>/dev/null | wc -l)"
        fatal="$(grep -Eic 'FATAL|MinDtFailed|OUTER-STEP-ROLLBACK|error:' "${LOG_FILE}" 2>/dev/null || true)"
        echo "  Health    : nonconv=${nonconv}; rollback=${rollback}; fatal/error=${fatal}"
    else
        echo "  Status    : waiting for log file..."
    fi

    echo "--------------------------------------------------------------------------"
    echo "  Model     : The COMPASS Land System Model"
    echo "  Authors   : Jizixing Zhu*, Xiaoya Shao, Xiaodong Gao, Wenjie Liu"
    echo "  Contact   : * Corresponding author | 3100362274@nwafu.edu.cn"
    echo "  Refresh   : every ${REFRESH_SEC}s; press Ctrl+C to exit"
    sleep "${REFRESH_SEC}"
done
