param(
    [Parameter(Mandatory)][string]$Src,
    [Parameter(Mandatory)][string]$OutDir,
    [string]$BuildType = "debug"
) 

# ── Colour helpers ─────────────────────────────────────────────────────────
function Write-Info    { param($msg) Write-Host "[INFO]  $msg" -ForegroundColor Cyan }
function Write-Success { param($msg) Write-Host "[OK]    $msg" -ForegroundColor Green }
function Write-Warn    { param($msg) Write-Host "[WARN]  $msg" -ForegroundColor Yellow }
function Write-Fail    { param($msg) Write-Host "[ERROR] $msg" -ForegroundColor Red }

# ── Find and initialize MSVC environment ──────────────────────────────────
$vcvars = "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat"

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

# ── Resolve build parameters ───────────────────────────────────────────────
$basename  = [System.IO.Path]::GetFileNameWithoutExtension($Src)
$libPath   = Join-Path (Split-Path $Src) "lib\windows"

$lib = if ($basename -match "v3") { "libfftw3-3.lib" } else { "libmyfft.lib" }

$flags = if ($BuildType -eq "release") {
    @("/O2", "/DNDEBUG")
} else {
    @("/Zi", "/Fd$OutDir\vc.pdb")
}

# ── Create output directory if needed ─────────────────────────────────────
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

# ── Show builing files ────────────────────────────────────────────────────
Write-Info "-------------------------"
Write-Info "- Building information: -"
Write-Info "-------------------------"
Write-Info "Source:     $basename.c"
Write-Info "Output:     $OutDir"
Write-Info "Build type: $BuildType"
Write-Info "Library:    $lib"
Write-Host ""

# ── Invoke cl.exe ──────────────────────────────────────────────────────────
$clArgs = @(
    "/EHsc", "/W4",
    $flags,
    "/Fo$OutDir\",
    "/Fe$OutDir\$basename.exe",
    $Src,
    "/link",
    "/LIBPATH:$libPath", $lib,
    "/OUT:$OutDir\$basename.exe",
    "/PDB:$OutDir\$basename.pdb",
    "/ILK:$OutDir\$basename.ilk",
    "/INCREMENTAL:NO"
)

Write-Info "Building '$basename' ($BuildType)..." -ForegroundColor Cyan
& cl.exe @clArgs

if ($LASTEXITCODE -ne 0) {
    Write-Fail "Build failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Success "Build succeeded: $OutDir\$basename.exe"