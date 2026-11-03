<#
Based on the COMPASS launchers by Jizixing Zhu.
Change: 2026-10-31 Wenjie Liu -- portable independent-copy startup.
Usage: run_compass_windows.ps1 [site-name-or-site-directory]
Env: MATLAB_BIN, COMPASS_OUTPUT, COMPASS_LOG_DIR, COMPASS_LOG_WINDOW=0,
     COMPASS_LAUNCH_WAIT=1 (wait and return this run's real completion status).
No scientific parameters or global MATLAB processes are changed.
#>
[CmdletBinding()]
param([Parameter(Position=0)][string]$SiteSelector = 'changwu')
$ErrorActionPreference = 'Stop'

function Quote-ProcessArgument([string]$Value) {
    if ($Value.Contains('"')) { throw 'A process argument contains a quote.' }
    return '"' + $Value + '"'
}

function Find-MatlabExecutable {
    if ($env:MATLAB_BIN) {
        $candidate = $env:MATLAB_BIN.Trim('"')
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
        $command = Get-Command $candidate -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($command) { return $command.Source }
        throw "MATLAB_BIN does not identify an executable: $candidate"
    }
    $command = Get-Command matlab -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($command) { return $command.Source }
    $roots = @($env:ProgramFiles, 'C:\Program Files', 'D:\Program Files') | Where-Object { $_ } | Select-Object -Unique
    foreach ($installRoot in $roots) {
        $matlabRoot = Join-Path $installRoot 'MATLAB'
        if (-not (Test-Path -LiteralPath $matlabRoot -PathType Container)) { continue }
        foreach ($versionDir in (Get-ChildItem -LiteralPath $matlabRoot -Directory | Sort-Object Name -Descending)) {
            $candidate = Join-Path $versionDir.FullName 'bin\matlab.exe'
            if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
        }
    }
    throw 'MATLAB was not found. Set MATLAB_BIN or put matlab.exe on PATH.'
}

function Read-ConfiguredOutput([string]$TomlPath) {
    $section = ''
    foreach ($line in [IO.File]::ReadAllLines($TomlPath)) {
        if ($line -match '^\s*\[([^\]]+)\]\s*(?:#.*)?$') { $section = $Matches[1]; continue }
        if ($section -ne 'output') { continue }
        if ($line -match '^\s*output_dir\s*=\s*("(?:[^"\\]|\\.)*"|''[^'']*'')\s*(?:#.*)?$') {
            $literal = $Matches[1]
            if ($literal.StartsWith('"')) { return [string]($literal | ConvertFrom-Json) }
            return $literal.Substring(1, $literal.Length - 2)
        }
    }
    throw 'site.toml must provide a quoted [output].output_dir.'
}

function Resolve-LocalDirectory([string]$Value, [string]$BaseDirectory) {
    if ([IO.Path]::IsPathRooted($Value)) { $candidate = $Value }
    else { $candidate = Join-Path $BaseDirectory $Value }
    return [IO.Path]::GetFullPath($candidate)
}

try {
    $projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
    if (-not $SiteSelector) { $SiteSelector = 'changwu' }
    if (Test-Path -LiteralPath (Join-Path $SiteSelector 'site.toml') -PathType Leaf) {
        $siteDir = (Resolve-Path -LiteralPath $SiteSelector).Path
    } else {
        $siteDir = Join-Path (Join-Path $projectRoot 'data\sites') $SiteSelector
    }
    $siteDir = [IO.Path]::GetFullPath($siteDir)
    $toml = Join-Path $siteDir 'site.toml'
    if (-not (Test-Path -LiteralPath $toml -PathType Leaf)) { throw "Site configuration not found: $toml" }
    $outputSource = 'site.toml:output.output_dir'
    $outputValue = Read-ConfiguredOutput $toml
    $outputBase = $siteDir
    if ($env:COMPASS_OUTPUT) { $outputValue = $env:COMPASS_OUTPUT; $outputBase = $projectRoot; $outputSource = 'env:COMPASS_OUTPUT' }
    elseif ($env:STEMMUS_OUTPUT) { $outputValue = $env:STEMMUS_OUTPUT; $outputBase = $projectRoot; $outputSource = 'env:STEMMUS_OUTPUT' }
    $resultDir = Resolve-LocalDirectory $outputValue $outputBase
    [IO.Directory]::CreateDirectory($resultDir) | Out-Null
    $logDir = Join-Path $resultDir 'logs'
    if ($env:COMPASS_LOG_DIR) { $logDir = Resolve-LocalDirectory $env:COMPASS_LOG_DIR $projectRoot }
    [IO.Directory]::CreateDirectory($logDir) | Out-Null
    $matlabExecutable = Find-MatlabExecutable
    $siteLabel = [IO.Path]::GetFileName($siteDir.TrimEnd('\','/')) -replace '[^\w.-]', '_'
    $runId = '{0}_{1}_{2}' -f $siteLabel, (Get-Date -Format 'yyyyMMdd_HHmmss'), ([Guid]::NewGuid().ToString('N').Substring(0,8))
    $logFile = Join-Path $logDir ('compass_' + $runId + '.log')
    $runReceipt = [IO.Path]::ChangeExtension($logFile, '.run.json')
    $launchReceipt = [IO.Path]::ChangeExtension($logFile, '.launch.json')
    $env:COMPASS_LAUNCH_SITE = $siteDir
    $env:COMPASS_RUN_ID = $runId
    $env:COMPASS_RUN_RECEIPT = $runReceipt
    $env:COMPASS_LAUNCH_LOG = $logFile
    $env:COMPASS_LAUNCH_OUTPUT_SOURCE = $outputSource
    $env:COMPASS_OUTPUT = $resultDir
    if (-not $env:COMPASS_DIAGNOSTIC_PROFILE -and -not $env:STEMMUS_DIAGNOSTIC_PROFILE) {
        $env:COMPASS_DIAGNOSTIC_PROFILE = 'full'
    }
    $batchCode = "addpath('scripts'); compass_batch_entry;"
    $argumentLine = '-nosplash -nodesktop -sd ' + (Quote-ProcessArgument $projectRoot) +
        ' -batch ' + (Quote-ProcessArgument $batchCode) + ' -logfile ' + (Quote-ProcessArgument $logFile)
    $process = Start-Process -FilePath $matlabExecutable -ArgumentList $argumentLine -WorkingDirectory $projectRoot -WindowStyle Hidden -PassThru
    $launch = [ordered]@{schema_version=1; run_id=$runId; launched_pid=$process.Id; matlab_executable=$matlabExecutable;
        project_root=$projectRoot; site_dir=$siteDir; output_dir=$resultDir; log_file=$logFile;
        run_receipt=$runReceipt; launched_utc=[DateTime]::UtcNow.ToString('o'); wait_requested=($env:COMPASS_LAUNCH_WAIT -eq '1')}
    [IO.File]::WriteAllText($launchReceipt, ($launch | ConvertTo-Json -Depth 5), [Text.UTF8Encoding]::new($false))
    Write-Host "COMPASS-MEND launched in background."
    Write-Host "Launched PID: $($process.Id)"
    Write-Host "Log: $logFile"
    Write-Host "Receipt: $runReceipt"
    if ($env:COMPASS_LAUNCH_WAIT -eq '1') {
        $process.WaitForExit()
        $exitCode = $process.ExitCode
        if ($null -eq $exitCode) { throw 'The launched process returned no exit status.' }
        if ($exitCode -ne 0) { Write-Error "MATLAB exited with code $exitCode. See $logFile" -ErrorAction Continue; exit $exitCode }
        if (-not (Test-Path -LiteralPath $runReceipt -PathType Leaf)) { throw 'MATLAB exited without this run completion receipt.' }
        $finished = [IO.File]::ReadAllText($runReceipt) | ConvertFrom-Json
        if ($finished.status -ne 'completed') { throw "The run did not complete: $($finished.status). See $runReceipt" }
        Write-Host "Completed: KT=$($finished.last_completed_kt)/$($finished.target_steps)"
        exit 0
    }
    Start-Sleep -Seconds 2
    $process.Refresh()
    if ($process.HasExited) {
        $process.WaitForExit()
        if ($process.ExitCode -ne 0) { throw "MATLAB exited during startup. See $logFile" }
        if (-not (Test-Path -LiteralPath $runReceipt)) { throw "MATLAB exited without a run receipt. See $logFile" }
        $finished = [IO.File]::ReadAllText($runReceipt) | ConvertFrom-Json
        if ($finished.status -ne 'completed') { throw "MATLAB stopped during startup. See $runReceipt" }
        Write-Host "Completed: KT=$($finished.last_completed_kt)/$($finished.target_steps)"
        exit 0
    }
    if ($env:COMPASS_LOG_WINDOW -ne '0') {
        $normalSite = Join-Path (Join-Path $projectRoot 'data\sites') $siteLabel
        $gui = Join-Path $PSScriptRoot 'compass_monitor_gui.py'
        if ($normalSite -eq $siteDir -and (Test-Path -LiteralPath $gui)) {
            $python = Get-Command python -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
            $pythonArgs = ''
            if (-not $python) { $python = Get-Command py -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1; $pythonArgs = '-3 ' }
            if ($python) {
                $guiArgs = $pythonArgs + (Quote-ProcessArgument $gui) + ' ' + (Quote-ProcessArgument $siteLabel) + ' ' + (Quote-ProcessArgument $logFile) + ' ' + $process.Id
                Start-Process -FilePath $python.Source -ArgumentList $guiArgs -WorkingDirectory $projectRoot -WindowStyle Hidden | Out-Null
            }
        }
    }
    exit 0
} catch {
    Write-Error $_ -ErrorAction Continue
    exit 1
}
