#!/usr/bin/env python3
"""COMPASS progress GUI (Linux + Windows) with project logo and live log.

Change: 2026-10-31 Wenjie Liu -- bind progress/PID to the launch-specific
sidecar and accepted log rows; ignore unrelated processes and stale checkpoints.

Usage:
  python scripts/compass_monitor_gui.py [site] [log_file] [pid]
  Double-click monitor.bat on Windows (defaults site=changwu, latest log).
"""
from __future__ import annotations

import os
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import tkinter as tk
from tkinter import font as tkfont
from tkinter import ttk

try:
    from PIL import Image, ImageTk
except Exception:  # pragma: no cover - optional
    Image = ImageTk = None

ROOT = Path(__file__).resolve().parents[1]
SITE = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] else "changwu"
# Empty string means "auto-pick newest log"; only set when caller passed a path.
LOG_ARG = Path(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].strip() else None
PID = int(sys.argv[3]) if len(sys.argv) > 3 and str(sys.argv[3]).isdigit() else 0
CFG = ROOT / "data" / "sites" / SITE / "site.toml"
LOGO = ROOT / "docs" / "figures" / "compass_logo.png"
REFRESH_MS = int(os.environ.get("COMPASS_MONITOR_REFRESH_MS", "5000"))


def ui_fonts():
    """Pick fonts that exist on Windows and Linux."""
    families = set(tkfont.families())
    if os.name == "nt":
        ui = "Segoe UI" if "Segoe UI" in families else "Microsoft YaHei UI"
        mono = "Consolas" if "Consolas" in families else "Courier New"
    else:
        ui = "DejaVu Sans" if "DejaVu Sans" in families else "Sans"
        mono = "DejaVu Sans Mono" if "DejaVu Sans Mono" in families else "Monospace"
    return ui, mono


def cfg(key, quoted=False):
    if not CFG.exists():
        return ""
    pat = rf'^\s*{re.escape(key)}\s*=\s*' + (r'"([^"]*)"' if quoted else r'([^ #]+)')
    for line in CFG.read_text(errors="replace").splitlines():
        m = re.match(pat, line)
        if m:
            return m.group(1).strip()
    return ""


def output_dir() -> Path:
    value = cfg("output_dir", True)
    if not value:
        return ROOT / "output" / SITE
    p = Path(value)
    return p if p.is_absolute() else (CFG.parent / p).resolve()


def resolve_log() -> Path:
    """Prefer explicit log path; else newest compass_<site>_*.log under output/logs."""
    if LOG_ARG is not None:
        return LOG_ARG
    log_dir = output_dir() / "logs"
    if not log_dir.is_dir():
        return log_dir / f"compass_{SITE}_(waiting).log"
    candidates = sorted(
        log_dir.glob(f"compass_{SITE}_*.log"),
        key=lambda p: p.stat().st_mtime if p.exists() else 0,
        reverse=True,
    )
    if candidates:
        return candidates[0]
    return log_dir / f"compass_{SITE}_(waiting).log"


def _pid_running(pid: int) -> bool:
    if not pid:
        return False
    if os.name == "nt":
        try:
            out = subprocess.check_output(
                ["tasklist", "/FI", f"PID eq {pid}", "/NH"],
                text=True,
                errors="replace",
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            return str(pid) in out and "INFO:" not in out.upper()
        except Exception:
            return False
    return os.path.exists(f"/proc/{pid}")


def read_run_receipt(log_path: Path):
    """Read only this log's sidecar; never infer activity from other MATLABs."""
    path = log_path.with_suffix(".run.json")
    try:
        record = json.loads(path.read_text(encoding="utf-8-sig"))
        if Path(record.get("project_root", "")).resolve() != ROOT.resolve():
            return {}
        if Path(record.get("log_file", "")).resolve() != log_path.resolve():
            return {}
        return record
    except (OSError, ValueError, TypeError):
        return {}


def read_log_tail(path: Path, n_lines: int = 400) -> str:
    if not path or not path.exists():
        return ""
    try:
        # Binary read + decode avoids Windows shared-write issues better than Path.read_text
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            chunk = min(size, 256_000)
            fh.seek(-chunk, os.SEEK_END)
            data = fh.read().decode("utf-8", errors="replace")
        lines = data.splitlines()
        return "\n".join(lines[-n_lines:])
    except Exception as exc:
        return f"(log read error: {exc})"


def latest_state(log_path: Path):
    total = 8760
    try:
        import datetime as dt
        start = dt.datetime.fromisoformat(cfg("start", True))
        end = dt.datetime.fromisoformat(cfg("end", True))
        total = int((end - start).total_seconds() / float(cfg("dt")))
    except (ValueError, TypeError, ZeroDivisionError):
        pass
    text = read_log_tail(log_path)
    matched = re.findall(r"\[COMPASS\] KT=(\d+)/(\d+)[^\r\n]*accepted=1(?:\D|$)", text)
    kt = int(matched[-1][0]) if matched else 0
    if matched:
        total = int(matched[-1][1])
    record = read_run_receipt(log_path)
    for key in ("target_steps", "configured_steps"):
        value = record.get(key)
        if isinstance(value, (int, float)) and value > 0:
            total = int(value)
            break
    if record.get("status") in ("completed", "failed"):
        kt = int(record.get("last_completed_kt", kt))
    run_pid = record.get("matlab_pid", PID)
    running = bool(run_pid and _pid_running(int(run_pid)))
    if record.get("status") in ("completed", "failed"):
        running = False
    rate = None
    # Rate uses only checkpoints written since this run began.
    try:
        import datetime as dt
        started = dt.datetime.fromisoformat(record["started_utc"].replace("Z", "+00:00")).timestamp()
        folder = Path(record["output_dir"]) / "checkpoints"
        points = sorted((p for p in folder.glob("checkpoint_KT*.mat") if p.stat().st_mtime >= started), key=lambda p: p.stat().st_mtime)
        if len(points) >= 2:
            first, last = points[-2:]
            a, b = re.search(r"KT0*(\d+)\.mat$", first.name), re.search(r"KT0*(\d+)\.mat$", last.name)
            elapsed = last.stat().st_mtime - first.stat().st_mtime
            if a and b and elapsed > 0:
                rate = (int(b.group(1)) - int(a.group(1))) * 60 / elapsed
    except (KeyError, ValueError, OSError, TypeError):
        pass
    return total, kt, rate, running, text



def run_status_label(log_path: Path, total: int, kt: int, running: bool):
    """Completion requires a completed receipt, including successful flush."""
    record = read_run_receipt(log_path)
    if record.get("status") == "failed":
        return "FAILED / CHECK LOG"
    if record.get("status") == "completed":
        if int(record.get("last_completed_kt", -1)) != int(record.get("target_steps", -2)):
            return "INCOMPLETE / CHECK LOG"
        return "COMPLETED" if record.get("full_window_completed") else "COMPLETED (SHORT RUN)"
    return "RUNNING" if running else "STOPPED / CHECK LOG"


def build():
    log_path = resolve_log()
    ui, mono = ui_fonts()

    root = tk.Tk()
    root.title(f"COMPASS Land System Model — {SITE} 实时进度")
    root.geometry("1280x760+40+40")
    root.minsize(960, 600)
    root.configure(bg="#102033")
    root.attributes("-topmost", True)
    root.after(1800, lambda: root.attributes("-topmost", False))
    root.lift()
    try:
        root.focus_force()
    except Exception:
        pass

    try:
        if Image and LOGO.exists():
            root.iconphoto(True, ImageTk.PhotoImage(Image.open(LOGO).resize((64, 64))))
    except Exception:
        pass

    left = tk.Frame(root, bg="#102033", width=280)
    left.pack(side="left", fill="y")
    left.pack_propagate(False)
    if Image and LOGO.exists():
        try:
            image = Image.open(LOGO).convert("RGBA")
            image.thumbnail((175, 175))
            logo = ImageTk.PhotoImage(image)
            tk.Label(left, image=logo, bg="#102033").pack(pady=(42, 14))
            left.logo = logo
        except Exception:
            pass
    tk.Label(
        left,
        text="COMPASS",
        fg="#8ee6ff",
        bg="#102033",
        font=(ui, 18, "bold"),
    ).pack(fill="x")
    tk.Label(
        left,
        text="Coupled Open Multiphysics Platform\nfor Atmosphere–Soil Systems",
        fg="#d8e6f0",
        bg="#102033",
        justify="center",
        wraplength=250,
        font=(ui, 10),
    ).pack(fill="x", pady=12)
    tk.Label(
        left,
        text="Jizixing Zhu*\nXiaoya Shao · Xiaodong Gao\nWenjie Liu",
        fg="#9eb4c7",
        bg="#102033",
        justify="center",
        wraplength=250,
        font=(ui, 9),
    ).pack(fill="x", side="bottom", pady=24)

    right = tk.Frame(root, bg="#f4f7fa", padx=24, pady=20)
    right.pack(side="right", expand=True, fill="both")

    tk.Label(
        right,
        text=f"{SITE} · Full-Year Run Progress",
        bg="#f4f7fa",
        fg="#102033",
        font=(ui, 18, "bold"),
        anchor="w",
    ).pack(fill="x")

    progress = tk.DoubleVar(value=0.0)
    ttk.Progressbar(right, variable=progress, maximum=100, length=520).pack(fill="x", pady=(18, 12))

    status = tk.StringVar()
    details = tk.StringVar()
    health = tk.StringVar()
    tk.Label(
        right,
        textvariable=status,
        bg="#f4f7fa",
        fg="#0b6b8a",
        font=(ui, 14, "bold"),
        anchor="w",
        wraplength=900,
        justify="left",
    ).pack(fill="x")
    tk.Label(
        right,
        textvariable=details,
        bg="#f4f7fa",
        fg="#263746",
        justify="left",
        anchor="nw",
        wraplength=900,
        font=(mono, 10),
        pady=8,
    ).pack(fill="x")
    tk.Label(
        right,
        textvariable=health,
        bg="#f4f7fa",
        fg="#4b5c68",
        justify="left",
        anchor="w",
        wraplength=900,
        font=(mono, 9),
    ).pack(fill="x", pady=(0, 8))

    log_frame = tk.LabelFrame(
        right,
        text=" Live log (tail) ",
        bg="#f4f7fa",
        fg="#263746",
        font=(ui, 10, "bold"),
        padx=6,
        pady=6,
    )
    log_frame.pack(fill="both", expand=True)
    log_text = tk.Text(
        log_frame,
        wrap="none",
        bg="#0d1b2a",
        fg="#c8e6c9",
        insertbackground="#c8e6c9",
        font=(mono, 9),
        height=18,
        borderwidth=0,
        highlightthickness=0,
    )
    yscroll = ttk.Scrollbar(log_frame, orient="vertical", command=log_text.yview)
    xscroll = ttk.Scrollbar(log_frame, orient="horizontal", command=log_text.xview)
    log_text.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
    log_text.grid(row=0, column=0, sticky="nsew")
    yscroll.grid(row=0, column=1, sticky="ns")
    xscroll.grid(row=1, column=0, sticky="ew")
    log_frame.rowconfigure(0, weight=1)
    log_frame.columnconfigure(0, weight=1)
    log_text.configure(state="disabled")

    state = {"log": log_path, "follow": True}

    def on_scroll(*_):
        # User scrolled up → stop auto-follow until they jump to end
        try:
            y0, y1 = log_text.yview()
            state["follow"] = y1 >= 0.999
        except Exception:
            pass

    log_text.bind("<ButtonRelease-1>", on_scroll)
    log_text.bind("<MouseWheel>", on_scroll)

    def refresh():
        # Auto mode: always follow newest site log (picks up a new run).
        if LOG_ARG is None:
            state["log"] = resolve_log()

        total, kt, rate, running, text = latest_state(state["log"])
        pct = 100.0 * kt / total if total else 0.0
        progress.set(min(100.0, max(0.0, pct)))
        run_state = run_status_label(state["log"], total, kt, running)
        run_record = read_run_receipt(state["log"])
        actual_output = Path(run_record.get("output_dir") or output_dir())
        status.set(f"KT={kt}/{total}   ({pct:.2f}%)   {run_state}")

        site_name = cfg("name", True) or SITE
        lat, lon, elev = cfg("latitude"), cfg("longitude"), cfg("elevation")
        speed = f"{rate:.2f} KT/min" if rate else "calculating from checkpoints"
        if run_state.startswith("FAILED"):
            eta = "check failure log"
        elif kt >= total and total > 0:
            eta = "0.0 min"
        elif rate:
            eta = f"{(total - kt) / rate:.1f} min"
        else:
            eta = "pending"
        details.set(
            f"Site      : {site_name}\n"
            f"Location  : {lat} N, {lon} E, {elev} m\n"
            f"Output    : {actual_output}\n"
            f"Log       : {state['log']}\n"
            f"Speed     : {speed}\n"
            f"ETA       : {eta}\n"
            f"Updated   : {time.strftime('%Y-%m-%d %H:%M:%S')}"
        )
        out = actual_output
        ckpt = out / "checkpoints"
        n_non = len(list(ckpt.glob("nonconv_diag_KT*.mat"))) if ckpt.is_dir() else 0
        n_rb = len(list(ckpt.glob("resume_before_KT*.mat"))) if ckpt.is_dir() else 0
        fatal = len(re.findall(r"FATAL|MinDtFailed|OUTER-STEP-ROLLBACK|error:", text, re.I))
        health.set(f"Health: nonconv={n_non}   rollback={n_rb}   fatal/error={fatal}")

        log_text.configure(state="normal")
        log_text.delete("1.0", "end")
        log_text.insert("end", text if text else "(waiting for log…)")
        if state["follow"]:
            log_text.see("end")
        log_text.configure(state="disabled")
        root.after(REFRESH_MS, refresh)

    refresh()
    root.mainloop()


if __name__ == "__main__":
    build()
