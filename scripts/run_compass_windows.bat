@echo off
setlocal
REM Based on the COMPASS launcher by Jizixing Zhu.
REM Change: 2026-10-31 Wenjie Liu -- portable launcher; exact process receipt.
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_compass_windows.ps1" %*
exit /b %ERRORLEVEL%
