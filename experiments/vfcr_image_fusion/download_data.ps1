$ErrorActionPreference = 'Stop'

$root = $PSScriptRoot
$data = Join-Path $root 'data'
$raw = Join-Path $data 'raw'
$bsds = Join-Path $data 'BSDS300'
$sipi = Join-Path $data 'SIPI'
$brain = Join-Path $data 'BrainWeb'

New-Item -ItemType Directory -Force -Path $raw,$bsds,$sipi,$brain | Out-Null

$bfile = Join-Path $raw 'BSDS300-images.tgz'
$sfile = Join-Path $raw 'misc.zip'
$mfile = Join-Path $brain 't1_icbm_normal_1mm_pn0_rf0.rawb'

if (-not (Test-Path -LiteralPath $bfile) -or (Get-Item -LiteralPath $bfile).Length -lt 1000000) {
    & curl.exe -L --fail --retry 3 --connect-timeout 20 'https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/BSDS300-images.tgz' -o $bfile
    if ($LASTEXITCODE -ne 0) { throw 'BSDS300 download failed.' }
}

if (-not (Test-Path -LiteralPath $sfile) -or (Get-Item -LiteralPath $sfile).Length -lt 1000000) {
    & curl.exe -L --fail --retry 3 --connect-timeout 20 'https://sipi.usc.edu/database/misc.zip' -o $sfile
    if ($LASTEXITCODE -ne 0) { throw 'SIPI download failed.' }
}

if (-not (Test-Path -LiteralPath $mfile) -or (Get-Item -LiteralPath $mfile).Length -ne 7109137) {
    $body = @{
        do_download_alias = 'T1 ICBM normal 1mm pn0 rf0'
        format_value = 'raw_byte'
        zip_value = 'none'
        who_name = ''
        who_institution = ''
        who_email = ''
        download_for_real = '[Start download!]'
    }
    Invoke-WebRequest -Method Post -Uri 'https://brainweb.bic.mni.mcgill.ca/cgi/brainweb1' -Body $body -OutFile $mfile
}

if (-not (Get-ChildItem -LiteralPath $bsds -Recurse -Filter '*.jpg' -ErrorAction SilentlyContinue)) {
    tar -xzf $bfile -C $bsds
}

if (-not (Test-Path -LiteralPath (Join-Path $sipi 'boat.512.tiff'))) {
    Expand-Archive -LiteralPath $sfile -DestinationPath $sipi -Force
}

$boat = Get-ChildItem -LiteralPath $sipi -Recurse -Filter 'boat.512.tiff' | Select-Object -First 1
if ($boat -and $boat.DirectoryName -ne $sipi) {
    Copy-Item -LiteralPath $boat.FullName -Destination (Join-Path $sipi 'boat.512.tiff') -Force
}

$files = @(
    Get-ChildItem -LiteralPath $bsds -Recurse -Filter '385028.jpg' | Select-Object -First 1
    Get-ChildItem -LiteralPath $bsds -Recurse -Filter '113016.jpg' | Select-Object -First 1
    Get-Item -LiteralPath (Join-Path $sipi 'boat.512.tiff')
    Get-Item -LiteralPath $mfile
)

if ($files.Count -ne 4 -or $files -contains $null) {
    throw 'Dataset files are incomplete.'
}

$files | Get-FileHash -Algorithm SHA256 | Select-Object Path,Hash
