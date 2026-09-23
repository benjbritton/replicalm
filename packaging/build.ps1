<#
    Build the Replicalm Windows installer.

    Four steps, each skippable so a failed run does not repeat the slow parts:

        1  create the pinned conda environment from environment.yml
        2  conda-pack it into stage\env
        3  stage the package, launcher and entry points into stage\app
        4  call Inno Setup on replicalm.iss

    Why a packed environment rather than a frozen executable: GDAL and PDAL
    carry native data directories, proj.db chief among them, and a frozen build
    that loses them fails at write time with an unprojected raster -- after the
    processing is done. A packed environment keeps the layout those libraries
    expect, at the cost of size.

    Requirements
        conda           any recent Miniconda or Anaconda
        conda-pack      conda install -n base -c conda-forge conda-pack
        Inno Setup 6    https://jrsoftware.org/isdl.php  (for step 4 only)

    Usage
        .\build.ps1                     full build
        .\build.ps1 -SkipEnv            reuse an existing replicalm environment
        .\build.ps1 -SkipInstaller      stop after staging, for testing
#>
[CmdletBinding()]
param(
    [string]$EnvName = "replicalm",
    [string]$Version = "1.0.0",
    [switch]$SkipEnv,
    [switch]$SkipInstaller
)

$ErrorActionPreference = "Stop"
$here  = Split-Path -Parent $MyInvocation.MyCommand.Path
$root  = Split-Path -Parent $here
$stage = Join-Path $here "stage"

function Step($n, $text) { Write-Host "`n[$n] $text" -ForegroundColor Cyan }

# Resolve a tool without requiring an activated conda shell. conda-pack lives in
# the base environment's Scripts directory, which is only on PATH inside an
# Anaconda Prompt; a plain PowerShell session will not see it, and failing there
# tells the user to install something they already have.
function Need($exe, $hint) {
    $c = Get-Command $exe -ErrorAction SilentlyContinue
    if ($c) { return $c.Source }
    $roots = @()
    $conda = Get-Command conda -ErrorAction SilentlyContinue
    if ($conda) { $roots += Split-Path -Parent $conda.Source }
    if ($env:CONDA_PREFIX) { $roots += (Join-Path $env:CONDA_PREFIX "Scripts") }
    $roots += @("$env:USERPROFILE\anaconda3\Scripts",
                "$env:USERPROFILE\miniconda3\Scripts",
                "$env:USERPROFILE\Anaconda3\Scripts",
                "$env:LOCALAPPDATA\Continuum\anaconda3\Scripts",
                "C:\ProgramData\anaconda3\Scripts",
                "C:\ProgramData\miniconda3\Scripts")
    foreach ($r in ($roots | Where-Object { $_ } | Select-Object -Unique)) {
        foreach ($ext in @(".exe", ".bat", "")) {
            $p = Join-Path $r ($exe + $ext)
            if (Test-Path $p) { return $p }
        }
    }
    throw "$exe not found on PATH or beside conda. $hint"
}

# --- 1: the environment ----------------------------------------------------
if (-not $SkipEnv) {
    Step 1 "Creating the pinned conda environment '$EnvName'"
    $conda = Need "conda" "Install Miniconda, or run this from an Anaconda Prompt."
    $exists = (& $conda env list) -match "^\s*$EnvName\s"
    if ($exists) {
        Write-Host "    environment exists; updating to match environment.yml"
        & $conda env update -n $EnvName -f (Join-Path $here "environment.yml") --prune
    } else {
        & $conda env create -n $EnvName -f (Join-Path $here "environment.yml")
    }
    if ($LASTEXITCODE -ne 0) { throw "conda env creation failed" }
} else {
    Step 1 "Skipping environment creation (-SkipEnv)"
}

# --- 2: pack it ------------------------------------------------------------
Step 2 "Packing the environment"
$condaPack = Need "conda-pack" "conda install -n base -c conda-forge conda-pack"
Write-Host "    using $condaPack"
$envDir = Join-Path $stage "env"
if (Test-Path $envDir) { Remove-Item $envDir -Recurse -Force }
New-Item -ItemType Directory -Force -Path $envDir | Out-Null

# conda-pack always writes an archive; there is no option to emit a directory.
# Installers need the tree, so it is packed and then expanded. ZipFile is used
# rather than Expand-Archive because the latter is very slow at this size.
$zip = Join-Path $stage "env.zip"
if (Test-Path $zip) { Remove-Item $zip -Force }
& $condaPack -n $EnvName --output $zip --force
if ($LASTEXITCODE -ne 0) { throw "conda-pack failed" }
Write-Host ("    archive: {0:N0} MB, expanding" -f ((Get-Item $zip).Length / 1MB))
Add-Type -AssemblyName System.IO.Compression.FileSystem
[System.IO.Compression.ZipFile]::ExtractToDirectory($zip, $envDir)
Remove-Item $zip -Force

$unpack = Join-Path $envDir "Scripts\conda-unpack.exe"
if (-not (Test-Path $unpack)) {
    throw "conda-unpack.exe missing from the packed environment. Without it " +
          "the installed copy keeps this machine's paths."
}
$size = (Get-ChildItem $envDir -Recurse -File | Measure-Object Length -Sum).Sum / 1MB
Write-Host ("    packed environment: {0:N0} MB" -f $size)

# --- 3: stage the application ---------------------------------------------
Step 3 "Staging the application"
$appDir = Join-Path $stage "app"
if (Test-Path $appDir) { Remove-Item $appDir -Recurse -Force }
New-Item -ItemType Directory -Force -Path (Join-Path $appDir "src") | Out-Null

Copy-Item (Join-Path $root "src\replicalm") (Join-Path $appDir "src\replicalm") `
    -Recurse -Force -Exclude "__pycache__"
Copy-Item (Join-Path $here "gui.py")  $appDir -Force
Copy-Item (Join-Path $root "config.py") $appDir -Force
# the installed README is written for someone using the program, not for
# someone reading the repository
Copy-Item (Join-Path $here "Replicalm-README.md") (Join-Path $appDir "README.md") -Force
if (Test-Path (Join-Path $root "docs")) {
    Copy-Item (Join-Path $root "docs") (Join-Path $appDir "docs") -Recurse -Force
}

# Launcher: pythonw so no console window appears behind the GUI.
@"
@echo off
rem Replicalm launcher. Uses the bundled environment; nothing is required on PATH.
set "REPLICALM_HOME=%~dp0"
set "PATH=%REPLICALM_HOME%env;%REPLICALM_HOME%env\Library\bin;%REPLICALM_HOME%env\Scripts;%PATH%"
set "GDAL_DATA=%REPLICALM_HOME%env\Library\share\gdal"
set "PROJ_LIB=%REPLICALM_HOME%env\Library\share\proj"
start "" "%REPLICALM_HOME%env\pythonw.exe" "%REPLICALM_HOME%gui.py" %*
"@ | Set-Content (Join-Path $appDir "Replicalm.cmd") -Encoding ASCII

# Command line entry point, same environment, console attached.
@"
@echo off
set "REPLICALM_HOME=%~dp0"
set "PATH=%REPLICALM_HOME%env;%REPLICALM_HOME%env\Library\bin;%REPLICALM_HOME%env\Scripts;%PATH%"
set "GDAL_DATA=%REPLICALM_HOME%env\Library\share\gdal"
set "PROJ_LIB=%REPLICALM_HOME%env\Library\share\proj"
set "PYTHONPATH=%REPLICALM_HOME%src;%PYTHONPATH%"
"%REPLICALM_HOME%env\python.exe" -m replicalm.pipeline %*
"@ | Set-Content (Join-Path $appDir "replicalm-cli.cmd") -Encoding ASCII

# --- licences --------------------------------------------------------------
$licDir = Join-Path $stage "licenses"
New-Item -ItemType Directory -Force -Path $licDir | Out-Null
$apache = Join-Path $stage "env\Lib\site-packages\rvt_py-2.2.1.dist-info\LICENSE"
if (Test-Path $apache) {
    Copy-Item $apache (Join-Path $licDir "Apache-2.0.txt") -Force
} else {
    Write-Warning "rvt-py LICENSE not found in the packed environment. Apache 2.0 requires it to accompany the distribution; fix before shipping."
}
Copy-Item (Join-Path $root "LICENSE") (Join-Path $licDir "Replicalm-MIT.txt") -Force

# --- 4: the installer ------------------------------------------------------
if ($SkipInstaller) {
    Step 4 "Skipping the installer (-SkipInstaller). Staged in $stage"
    exit 0
}
Step 4 "Building the installer"
# winget installs Inno Setup per-user by default, which is not under Program
# Files, so the registry is consulted as well as the usual locations.
$isccCandidates = @(
    "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
    "$env:ProgramFiles\Inno Setup 6\ISCC.exe",
    "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe"
)
foreach ($hive in @(
        "HKLM:\SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall",
        "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall",
        "HKCU:\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall")) {
    Get-ChildItem -Path $hive -ErrorAction SilentlyContinue | ForEach-Object {
        $prop = Get-ItemProperty $_.PSPath -ErrorAction SilentlyContinue
        if ($prop.DisplayName -like "*Inno Setup*" -and $prop.InstallLocation) {
            $isccCandidates += (Join-Path $prop.InstallLocation "ISCC.exe")
        }
    }
}
$iscc = $isccCandidates | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1
if (-not $iscc) {
    throw "Inno Setup 6 not found. Install it (winget install --id " +
          "JRSoftware.InnoSetup -e), or re-run with -SkipInstaller."
}
Write-Host "    using $iscc"
& $iscc "/DAppVersion=$Version" (Join-Path $here "replicalm.iss")
if ($LASTEXITCODE -ne 0) { throw "ISCC failed" }

$out = Join-Path $here "dist\Replicalm-$Version-setup.exe"
Write-Host "`nInstaller: $out" -ForegroundColor Green
if (Test-Path $out) {
    Write-Host ("Size: {0:N0} MB" -f ((Get-Item $out).Length / 1MB))
}
