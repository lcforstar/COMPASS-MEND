@echo off
setlocal
REM COMPASS Windows entry. Usage: run.bat [site]
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
  pause
)
exit /b %EC%
