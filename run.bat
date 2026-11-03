@echo off
setlocal
REM COMPASS Windows entry. Usage: run.bat [site]
REM Change: 2026-10-31 Wenjie Liu -- do not pause in explicit wait/CI mode.
REM ASCII-only: UTF-8 Chinese breaks cmd.exe on many Windows systems.
set "ROOT=%~dp0"
if "%~1"=="" (
  call "%ROOT%scripts\run_compass_windows.bat" changwu
) else (
  call "%ROOT%scripts\run_compass_windows.bat" %*
)
set "EC=%ERRORLEVEL%"
if not "%EC%"=="0" (
  echo.
  echo Launch failed, exit code %EC%.
  if not "%COMPASS_LAUNCH_WAIT%"=="1" pause
)
exit /b %EC%
