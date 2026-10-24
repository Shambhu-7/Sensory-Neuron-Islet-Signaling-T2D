
library(GEOquery)

raw_dir <- file.path(
    getwd(),
    "data",
    "raw",
    "GSE295206"
)

cat("Downloading GSM8943583 annotation...\n")

GEOquery::getGEOSuppFiles(
    "GSM8943583",
    baseDir = raw_dir,
    makeDirectory = TRUE,
    fetch_files = TRUE,
    filter_regex = "sample4\\.csv\\.gz$"
)

cat("Downloading GSM8943584 annotation...\n")

GEOquery::getGEOSuppFiles(
    "GSM8943584",
    baseDir = raw_dir,
    makeDirectory = TRUE,
    fetch_files = TRUE,
    filter_regex = "sample8\\.csv\\.gz$"
)

cat("\nDone.\n")

