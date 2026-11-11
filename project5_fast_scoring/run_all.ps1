# =============================================================================
# run_all.ps1 -- improvement 5: runs the tests and every experiment, one after the other.
#
# Usage (from this folder):
#   powershell -ExecutionPolicy Bypass -File .\run_all.ps1                 # everything
#   powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -SkipInstall    # do not reinstall scoth.fast
#   powershell -ExecutionPolicy Bypass -File .\run_all.ps1 -SkipTargeted   # skip steps 3 and 4 (the long ones)
#
# Needs project0_reproduction next to this folder: package scoth installed, and
# ..\project0_reproduction\output-data made by its experiments\main.R.
# =============================================================================
param([switch]$SkipInstall, [switch]$SkipTargeted)

if (-not (Get-Command Rscript.exe -ErrorAction SilentlyContinue)) {
  $RBin = "$env:LOCALAPPDATA\Programs\R\R-4.6.2\bin"
  if (-not (Test-Path "$RBin\Rscript.exe")) { $RBin = "$env:ProgramFiles\R\R-4.6.2\bin" }
  $env:Path = "$RBin;" + $env:Path
}
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force logs, results | Out-Null

function Invoke-Step($title, $command) {
  Write-Host ""
  Write-Host "==== $title ====" -ForegroundColor Cyan
  $t0 = Get-Date
  cmd.exe /c $command
  $code = $LASTEXITCODE
  $sec = [math]::Round(((Get-Date) - $t0).TotalSeconds)
  if ($code -ne 0) { Write-Host "FAILED (exit code $code) after $sec s" -ForegroundColor Red; exit $code }
  Write-Host "OK ($sec s)" -ForegroundColor Green
}

# ---- what must exist first ---------------------------------------------------
& Rscript.exe -e "if (!requireNamespace('scoth', quietly = TRUE)) quit(status = 1)"
if ($LASTEXITCODE -ne 0) {
  Write-Host "Package scoth not found: install project0_reproduction first." -ForegroundColor Red; exit 1
}
if (-not (Test-Path "..\project0_reproduction\output-data\summarizing_results_TCPD\summary.RDS")) {
  Write-Host "..\project0_reproduction\output-data not found: run experiments\main.R in project0_reproduction first." -ForegroundColor Red; exit 1
}

# ---- steps -------------------------------------------------------------------
if (-not $SkipInstall) {
  Invoke-Step "[1/8] install scoth.fast (next to scoth)" "R.exe CMD INSTALL --no-multiarch . > logs\install.log 2>&1"
}
Invoke-Step "[2/8] unit tests (21 original + 9 new)" "Rscript.exe experiments\run_tests.R"
if (-not $SkipTargeted) {
  Invoke-Step "[3/8] experiment A: scoring time vs length T (log: logs\targeted_runtime.log)" "Rscript.exe experiments\targeted_runtime.R > logs\targeted_runtime.log 2>&1"
  Invoke-Step "[4/8] experiment B: exactness on random series (log: logs\exactness_random.log)" "Rscript.exe experiments\exactness_random.R > logs\exactness_random.log 2>&1"
}
Invoke-Step "[5/8] exactness and ties on TCPD (log: logs\exactness_ties.log)" "Rscript.exe experiments\exactness_ties.R > logs\exactness_ties.log 2>&1"
Invoke-Step "[6/8] full TCPD pipeline with fast scoring (log: logs\main_fast_run.log)" "Rscript.exe experiments\main_fast.R > logs\main_fast_run.log 2>&1"
Invoke-Step "[7/8] three evaluation protocols" "Rscript.exe experiments\extra_protocols.R"
Invoke-Step "[8/8] side-by-side comparison with project0 (log: logs\compare_with_baseline.log)" "Rscript.exe experiments\compare_with_baseline.R > logs\compare_with_baseline.log 2>&1"

Write-Host ""
Write-Host "---- results\comparison_summary.txt ----"; Get-Content -Encoding UTF8 results\comparison_summary.txt
Write-Host "---- results\exactness_by_exact_max_len.csv ----"; Get-Content -Encoding UTF8 results\exactness_by_exact_max_len.csv
if (Test-Path results\runtime_summary.csv) { Write-Host "Figure: results\runtime_loglog.png" }
Write-Host "ALL STEPS FINISHED" -ForegroundColor Green
