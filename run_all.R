# Run from the project root, in RStudio or with Rscript run_all.R.
if (!file.exists('code/01_download.R')) stop('Open CellNavAI.Rproj first, or set the working directory to the project root.')
options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('ranger', quietly = TRUE)) install.packages('ranger')
if (!requireNamespace('ranger', quietly = TRUE)) stop('Install the ranger package and try again.')
for (d in c('data','figures','results','report')) dir.create(d, showWarnings = FALSE)
for (f in sort(list.files('code', pattern = '^[0-9].*\\.R$', full.names = TRUE))) {
  message('\nRunning ', f)
  source(f, local = .GlobalEnv, encoding = 'UTF-8')
}
message('\nComplete. Open README.md and results/metrics.csv; review all conclusions before submission.')
