# This script deploys the TaxoTox Shiny app to shinyapps.io.
# It creates a temporary directory, copies only the necessary files to it,
# and then deploys the app from that temporary directory.

# 1. Set up variables
$ErrorActionPreference = "Stop"
$projectRoot = $PSScriptRoot
$tempDir = Join-Path $projectRoot "taxotox_deploy"
$tempDirForR = $tempDir.Replace("\", "/")

# 2. Create a temporary directory
if (Test-Path $tempDir) {
    Remove-Item -Path $tempDir -Recurse -Force
}
New-Item -Path $tempDir -ItemType Directory

# 3. Copy the necessary files
Copy-Item -Path (Join-Path $projectRoot "Code/app.R") -Destination $tempDir
Copy-Item -Path (Join-Path $projectRoot "Data") -Destination $tempDir -Recurse

# 4. Deploy the app
$rScriptPath = "C:/Program Files/R/R-4.5.0/bin/Rscript.exe"
$deployCommand = "rsconnect::deployApp(appDir = '$tempDirForR', appName = 'TaxoTox', forceUpdate = TRUE, account = 'yairsuari')"

& $rScriptPath -e $deployCommand

# 5. Clean up the temporary directory
Remove-Item -Path $tempDir -Recurse -Force

Write-Host "Deployment complete."
