@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM COMPASS progress GUI. Usage: monitor.bat [site] [logfile] [pid]
REM ASCII-only for cmd.exe compatibility.
set "ROOT=%~dp0"
set "SITE=%~1"
if "!SITE!"=="" set "SITE=changwu"
set "LOG=%~2"
set "PID=%~3"
set "GUI=%ROOT%scripts\compass_monitor_gui.py"

if not exist "!GUI!" (
  echo GUI script not found: !GUI!
  pause
  exit /b 1
)

set "USE_PY="
where python >nul 2>&1
if not errorlevel 1 set "USE_PY=python"
if not defined USE_PY (
  where py >nul 2>&1
  if not errorlevel 1 set "USE_PY=py"
)
if not defined USE_PY (
  echo Python not found. Install Python 3 and check "Add python.exe to PATH".
  echo Then reopen the terminal and run monitor.bat again.
  echo GUI: !GUI!
  pause
  exit /b 1
)

echo Starting COMPASS progress GUI  site=!SITE!
if "!USE_PY!"=="py" (
  if not "!LOG!"=="" if not "!PID!"=="" (
    start "COMPASS-monitor-!SITE!" py -3 "!GUI!" "!SITE!" "!LOG!" "!PID!"
  ) else if not "!LOG!"=="" (
    start "COMPASS-monitor-!SITE!" py -3 "!GUI!" "!SITE!" "!LOG!"
  ) else (
    start "COMPASS-monitor-!SITE!" py -3 "!GUI!" "!SITE!"
  )
) else (
  if not "!LOG!"=="" if not "!PID!"=="" (
    start "COMPASS-monitor-!SITE!" python "!GUI!" "!SITE!" "!LOG!" "!PID!"
  ) else if not "!LOG!"=="" (
    start "COMPASS-monitor-!SITE!" python "!GUI!" "!SITE!" "!LOG!"
  ) else (
    start "COMPASS-monitor-!SITE!" python "!GUI!" "!SITE!"
  )
)
exit /b 0
