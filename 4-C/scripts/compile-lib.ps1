# ── Colour helpers ─────────────────────────────────────────────────────────
function Write-Info    { param($msg) Write-Host "[INFO]  $msg" -ForegroundColor Cyan }
function Write-Success { param($msg) Write-Host "[OK]    $msg" -ForegroundColor Green }
function Write-Warn    { param($msg) Write-Host "[WARN]  $msg" -ForegroundColor Yellow }
function Write-Fail    { param($msg) Write-Host "[ERROR] $msg" -ForegroundColor Red }

# ── Find and initialize MSVC environment ──────────────────────────────────
$vcvars = "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat"

if (-not (Test-Path $vcvars)) {
    Write-Error "! vcvarsall.bat not found at '$vcvars'. Please update the path in build script."
    exit 1
}

# Import the MSVC environment variables into the current PowerShell session
$envVars = cmd.exe /c "`"$vcvars`" x64 > nul 2>&1 && set" | ForEach-Object {
    if ($_ -match "^([^=]+)=(.*)$") {
        [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
    }
}

$ScriptDir = $PSScriptRoot
$WorkingDir = Join-Path $ScriptDir ".."
$BuildDir = Join-Path $WorkingDir "build\windows"
$OutDir = Join-Path $WorkingDir "lib\windows"

$clArgs = @(
    "/c",
    "/Fo$BuildDir\C_my_fft",
    "/O2", "/W4", "/MT",
    "$WorkingDir\C_my_fft.c"
)

$libArgs = @(
    "/OUT:$OutDir\libmyfft.lib",
    "$BuildDir\C_my_fft.obj"
)

& cl.exe $clArgs
if ($LASTEXITCODE -ne 0) {
    Write-Fail "Compile OBJ failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}
& lib.exe $libArgs
if ($LASTEXITCODE -ne 0) {
    Write-Fail "Archive to LIB failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Success "Build succeeded: $OutDir\libmyfft.lib"
