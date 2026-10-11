param([string]$Python = 'python', [string]$OutDir = '')
$ErrorActionPreference = 'Stop'
$projectDir = $PSScriptRoot
if (-not $OutDir) { $OutDir = Join-Path $projectDir ('demo-output-' + (Get-Date -Format 'yyyyMMdd-HHmmss')) }
Push-Location -LiteralPath $projectDir
try {
    & $Python -m crisprphibert --demo --outdir $OutDir
    if ($LASTEXITCODE -ne 0) { throw 'The demo did not complete.' }
    Write-Output ('Report: ' + (Join-Path $OutDir 'report.html'))
} finally { Pop-Location }
