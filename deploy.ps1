# deploy.ps1 - Deploy TaxoTox Shiny app to shinyapps.io
#
# Files deployed:
#   Code/app.R                     - main Shiny application
#   Code/taxotox-service.json      - Google Sheets service account key
#   Data/taxotox_reference.fst     - pre-aggregated LC50 / benchmark reference table
#   Data/Known_CAS.fst             - curated compound name -> CASRN lookup
#
# Files intentionally NOT deployed:
#   Data/final_ecotox_data.fst     - raw ECOTOX data (install-time only, 12 MB)
#   Data/Known_CAS_backup.fst      - backup copy
#   Data/dtxsid_map.csv            - API cache (install-time only)
#   Data/sample_*                  - test files
#   Code/taxotox_install*.R        - install scripts (not needed at runtime)
#   Code/taxotox_curate.R          - curation script (not needed at runtime)

$ErrorActionPreference = "Stop"
$projectRoot = $PSScriptRoot
$tempDir     = Join-Path $projectRoot "taxotox_deploy"
$tempDirForR = $tempDir.Replace("\", "/")
$rScript     = 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe'

# 1. Prepare staging directory
if (Test-Path $tempDir) {
    Remove-Item -Path $tempDir -Recurse -Force
}
New-Item -Path $tempDir        -ItemType Directory | Out-Null
New-Item -Path "$tempDir\Data" -ItemType Directory | Out-Null

# 2. Copy app files (flat — shinyapps.io serves from app root)
Copy-Item -Path (Join-Path $projectRoot "Code\app.R")               -Destination $tempDir
Copy-Item -Path (Join-Path $projectRoot "Code\taxotox-service.json") -Destination $tempDir

# 3. Copy required data files only
$dataFiles = @(
    "Data\taxotox_reference.fst",
    "Data\Known_CAS.fst"
)
foreach ($f in $dataFiles) {
    $src = Join-Path $projectRoot $f
    if (Test-Path $src) {
        Copy-Item -Path $src -Destination "$tempDir\Data\"
        Write-Host "  Copied: $f"
    } else {
        Write-Warning "  MISSING: $f - deploy may fail"
    }
}

# 4. Report staging size
$size = (Get-ChildItem $tempDir -Recurse | Measure-Object -Property Length -Sum).Sum
Write-Host ("Staging directory: {0:N0} KB" -f ($size / 1KB))

# 5. Deploy
$deployCmd = "options(renv.consent = TRUE, rsconnect.packrat.ignore.failures = TRUE); rsconnect::deployApp(appDir = '$tempDirForR', appName = 'TaxoTox', forceUpdate = TRUE, account = 'yairsuari')"
& "$rScript" -e $deployCmd

# 6. Clean up
Remove-Item -Path $tempDir -Recurse -Force
Write-Host "Deployment complete."
