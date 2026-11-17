# Download the official DeepLoc 2.0 training/validation data, preserving partitions.
raw_path <- 'data/Swissprot_Train_Validation_dataset.csv'
url <- 'https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/Swissprot_Train_Validation_dataset.csv'
classes <- c('Cytoplasm', 'Mitochondrion', 'Nucleus')
labels <- c('Cell membrane','Cytoplasm','Endoplasmic reticulum','Extracellular',
            'Golgi apparatus','Lysosome/Vacuole','Mitochondrion','Nucleus','Peroxisome','Plastid')
if (!file.exists(raw_path)) {
  options(timeout = max(600, getOption('timeout')))
  tmp <- tempfile(fileext = '.csv')
  tryCatch({
    status <- download.file(url, tmp, mode = 'wb', method = 'libcurl')
    if (status != 0 || !file.exists(tmp) || file.info(tmp)$size < 1000) stop('Incomplete download')
    if (!file.copy(tmp, raw_path)) stop('Could not save download')
  }, error = function(e) stop('Download failed. Open the official dataset URL in a browser, save it as ', raw_path, ', then rerun. Details: ', conditionMessage(e)))
  unlink(tmp)
}
raw <- read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c('ACC','Sequence','Partition',labels)
if (!all(required %in% names(raw))) stop('Unexpected CSV schema. Missing: ', paste(setdiff(required,names(raw)), collapse=', '))
writeLines(c(paste('Source:', url), paste('MD5:', tools::md5sum(raw_path)),
             paste('Raw rows:', nrow(raw)), paste('Columns:', paste(names(raw), collapse=', '))), 'results/data_provenance.txt')
