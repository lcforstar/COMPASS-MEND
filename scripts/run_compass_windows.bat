@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM COMPASS Windows launcher. Usage: scripts\run_compass_windows.bat [site]
REM Env: MATLAB_BIN, COMPASS_LOG_DIR, COMPASS_LOG_WINDOW=0
REM ASCII-only: UTF-8 non-ASCII breaks cmd.exe code pages.

set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%.." >nul 2>&1
set "PROJECT_ROOT=%CD%"
popd >nul 2>&1

set "SITE=%~1"
if "!SITE!"=="" set "SITE=changwu"
set "SITE_DIR=%PROJECT_ROOT%\data\sites\!SITE!"

if not exist "!SITE_DIR!\site.toml" (
  echo Site config not found: !SITE_DIR!\site.toml
  echo Available sites:
  if exist "%PROJECT_ROOT%\data\sites" (
    for /d %%D in ("%PROJECT_ROOT%\data\sites\*") do (
      if exist "%%D\site.toml" echo   %%~nxD
    )
  )
  exit /b 1
)

REM Parse output_dir from site.toml via PowerShell (reliable quoting)
set "CONFIGURED_OUTPUT="
for /f "usebackq delims=" %%T in (`powershell -NoProfile -Command "$p=Join-Path '!SITE_DIR!' 'site.toml'; $m=Select-String -Path $p -Pattern '^\s*output_dir\s*=\s*\"([^\"]*)\"' | Select-Object -First 1; if($m){$m.Matches[0].Groups[1].Value}"`) do set "CONFIGURED_OUTPUT=%%T"

if "!CONFIGURED_OUTPUT!"=="" (
  echo site.toml missing [output].output_dir: !SITE_DIR!\site.toml
  echo Set output_dir to an existing directory, then retry.
  exit /b 1
)

REM Absolute path if drive letter or UNC; else relative to site folder
set "RESULT_DIR="
echo.!CONFIGURED_OUTPUT!| findstr /R /C:"^[A-Za-z]:[\\/]" /C:"^\\\\" >nul
if !ERRORLEVEL! equ 0 (
  set "RESULT_DIR=!CONFIGURED_OUTPUT!"
) else (
  pushd "!SITE_DIR!" >nul 2>&1
  if not exist "!CONFIGURED_OUTPUT!\" (
    echo Output dir does not exist: !SITE_DIR!\!CONFIGURED_OUTPUT!
    echo Edit [output].output_dir in site.toml.
    popd >nul 2>&1
    exit /b 1
  )
  pushd "!CONFIGURED_OUTPUT!" >nul 2>&1
  set "RESULT_DIR=!CD!"
  popd >nul 2>&1
  popd >nul 2>&1
)

if not exist "!RESULT_DIR!\" (
  echo Output dir does not exist: !RESULT_DIR!
  echo Edit [output].output_dir in site.toml.
  exit /b 1
)

if defined COMPASS_LOG_DIR (
  set "LOG_DIR=!COMPASS_LOG_DIR!"
) else (
  set "LOG_DIR=!RESULT_DIR!\logs"
)

call :find_matlab
if not defined MATLAB_BIN (
  echo MATLAB not found.
  echo Install MATLAB or set MATLAB_BIN to matlab.exe full path.
  echo Example: set MATLAB_BIN=C:\Program Files\MATLAB\R2025b\bin\matlab.exe
  exit /b 1
)
if not exist "!MATLAB_BIN!" (
  echo MATLAB executable not found: !MATLAB_BIN!
  exit /b 1
)

if not exist "!LOG_DIR!\" mkdir "!LOG_DIR!" >nul 2>&1

set "TIMESTAMP="
for /f "usebackq delims=" %%I in (`powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"`) do set "TIMESTAMP=%%I"
if "!TIMESTAMP!"=="" set "TIMESTAMP=run"
set "LOG_FILE=!LOG_DIR!\compass_!SITE!_!TIMESTAMP!.log"

set "PROJECT_ROOT_FWD=!PROJECT_ROOT:\=/!"
set "MATLAB_CODE=cd('!PROJECT_ROOT_FWD!'); run_site('!SITE!');"

echo COMPASS Windows launcher > "!LOG_FILE!"
echo SITE=!SITE!>> "!LOG_FILE!"
echo RESULT_DIR=!RESULT_DIR!>> "!LOG_FILE!"
echo MATLAB=!MATLAB_BIN!>> "!LOG_FILE!"
echo.>> "!LOG_FILE!"

start "COMPASS-MATLAB-!SITE!" /b "!MATLAB_BIN!" -nosplash -nodesktop -minimize -batch "!MATLAB_CODE!" -logfile "!LOG_FILE!"

timeout /t 3 /nobreak >nul

set "MATLAB_PID="
for /f "tokens=2 delims=," %%P in ('tasklist /FI "IMAGENAME eq MATLAB.exe" /FO CSV /NH 2^>nul') do (
  if not defined MATLAB_PID set "MATLAB_PID=%%~P"
)

tasklist /FI "IMAGENAME eq MATLAB.exe" 2>nul | findstr /I "MATLAB.exe" >nul
if errorlevel 1 (
  echo COMPASS launch failed: MATLAB.exe not detected.
  echo Log: !LOG_FILE!
  echo ----- log tail -----
  powershell -NoProfile -Command "Get-Content -LiteralPath '!LOG_FILE!' -Tail 30 -ErrorAction SilentlyContinue"
  exit /b 1
)

echo COMPASS started in background
echo Site: !SITE!
echo Result dir: !RESULT_DIR!
echo Log: !LOG_FILE!
if defined MATLAB_PID echo MATLAB PID: !MATLAB_PID!
echo Live log: powershell -NoProfile -Command "Get-Content -LiteralPath '!LOG_FILE!' -Wait -Tail 50"

if /I "!COMPASS_LOG_WINDOW!"=="0" goto :eof

set "GUI_SCRIPT=%PROJECT_ROOT%\scripts\compass_monitor_gui.py"
set "HAVE_PYTHON=0"
where python >nul 2>&1
if not errorlevel 1 set "HAVE_PYTHON=1"
set "HAVE_PY=0"
where py >nul 2>&1
if not errorlevel 1 set "HAVE_PY=1"

if exist "!GUI_SCRIPT!" if "!HAVE_PYTHON!"=="1" (
  echo Opening progress GUI...
  if defined MATLAB_PID (
    start "COMPASS-monitor-!SITE!" python "!GUI_SCRIPT!" "!SITE!" "!LOG_FILE!" "!MATLAB_PID!"
  ) else (
    start "COMPASS-monitor-!SITE!" python "!GUI_SCRIPT!" "!SITE!" "!LOG_FILE!"
  )
  goto :eof
)
if exist "!GUI_SCRIPT!" if "!HAVE_PY!"=="1" (
  echo Opening progress GUI...
  if defined MATLAB_PID (
    start "COMPASS-monitor-!SITE!" py -3 "!GUI_SCRIPT!" "!SITE!" "!LOG_FILE!" "!MATLAB_PID!"
  ) else (
    start "COMPASS-monitor-!SITE!" py -3 "!GUI_SCRIPT!" "!SITE!" "!LOG_FILE!"
  )
  goto :eof
)

echo Python not found; GUI skipped. Install Python 3 or use monitor.bat later.
echo Log: !LOG_FILE!
start "COMPASS Log - !SITE!" powershell -NoProfile -NoExit -Command "Write-Host 'COMPASS log: !LOG_FILE!'; Get-Content -LiteralPath '!LOG_FILE!' -Wait -Tail 80"
goto :eof

:find_matlab
if defined MATLAB_BIN (
  if exist "!MATLAB_BIN!" exit /b 0
)
set "MATLAB_BIN="
if exist "F:\Matlab\bin\matlab.exe" (
  set "MATLAB_BIN=F:\Matlab\bin\matlab.exe"
  exit /b 0
)
for %%V in (R2025b R2025a R2024b R2024a R2023b R2023a) do (
  if exist "C:\Program Files\MATLAB\%%V\bin\matlab.exe" (
    set "MATLAB_BIN=C:\Program Files\MATLAB\%%V\bin\matlab.exe"
    exit /b 0
  )
  if exist "D:\Program Files\MATLAB\%%V\bin\matlab.exe" (
    set "MATLAB_BIN=D:\Program Files\MATLAB\%%V\bin\matlab.exe"
    exit /b 0
  )
)
where matlab >nul 2>&1
if not errorlevel 1 (
  for /f "delims=" %%M in ('where matlab 2^>nul') do (
    set "MATLAB_BIN=%%M"
    exit /b 0
  )
)
exit /b 1
