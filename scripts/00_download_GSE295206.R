
options(timeout = 7200)

# ============================================================
# PAPER 3
# Redox-associated sensory neuron-islet signaling networks
#
# Script 00:
# Download required processed GSE295206 files
#
# IMPORTANT:
# We deliberately avoid downloading the ~9 GB GSE RAW archive
# containing large CLOUPE files.
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

library(GEOquery)

project_dir <- normalizePath(".", mustWork = TRUE)

raw_dir <- file.path(
    project_dir,
    "data",
    "raw",
    "GSE295206"
)

dir.create(
    raw_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 16 Visium sections
# ------------------------------------------------------------

gsms <- paste0(
    "GSM",
    8943583:8943598
)

# Files needed for our analysis:
#
# 1. filtered_feature_bc_matrix.h5
#       processed Visium counts
#
# 2. NgnNeuSelections.csv.gz
#       author's manually selected Nageotte / neuronal barcodes
#
# 3. tissue_positions*.csv.gz
#       spatial coordinates
# ------------------------------------------------------------

wanted_pattern <- paste0(
    "(",
    "filtered_feature_bc_matrix\\.h5$",
    "|",
    "NgnNeuSelections\\.csv\\.gz$",
    "|",
    "tissue_positions.*\\.csv\\.gz$",
    ")"
)

cat("\n============================================\n")
cat("BUILDING GEO FILE INVENTORY\n")
cat("============================================\n\n")

inventory_list <- list()

for (gsm in gsms) {

    cat("Checking:", gsm, "\n")

    x <- GEOquery::getGEOSuppFiles(
        gsm,
        fetch_files = FALSE,
        filter_regex = wanted_pattern
    )

    if (!is.null(x) && nrow(x) > 0) {

        x$GSM <- gsm

        inventory_list[[gsm]] <- x

    } else {

        warning("No requested files found for ", gsm)

    }
}

inventory <- do.call(
    rbind,
    inventory_list
)

inventory <- inventory[
    ,
    c(
        "GSM",
        setdiff(
            colnames(inventory),
            "GSM"
        )
    )
]

write.csv(
    inventory,
    file.path(
        raw_dir,
        "GSE295206_download_inventory.csv"
    ),
    row.names = FALSE
)

cat("\nInventory contains ",
    nrow(inventory),
    " files.\n",
    sep = ""
)

cat("\n============================================\n")
cat("DOWNLOADING VISIUM FILES\n")
cat("============================================\n\n")

for (gsm in gsms) {

    cat("\n--------------------------------------------\n")
    cat("Downloading:", gsm, "\n")
    cat("--------------------------------------------\n")

    GEOquery::getGEOSuppFiles(
        gsm,
        baseDir = raw_dir,
        makeDirectory = TRUE,
        fetch_files = TRUE,
        filter_regex = wanted_pattern
    )
}

# ------------------------------------------------------------
# Download study-level metadata only.
# Do NOT download GSE295206_RAW.tar.
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DOWNLOADING STUDY METADATA\n")
cat("============================================\n\n")

GEOquery::getGEOSuppFiles(
    "GSE295206",
    baseDir = raw_dir,
    makeDirectory = TRUE,
    fetch_files = TRUE,
    filter_regex = "metadata_nageotteSamples"
)

# ------------------------------------------------------------
# Validate downloaded files
# ------------------------------------------------------------

all_files <- list.files(
    raw_dir,
    recursive = TRUE,
    full.names = TRUE
)

h5_files <- grep(
    "filtered_feature_bc_matrix\\.h5$",
    all_files,
    value = TRUE
)

selection_files <- grep(
    "NgnNeuSelections\\.csv\\.gz$",
    all_files,
    value = TRUE
)

position_files <- grep(
    "tissue_positions.*\\.csv\\.gz$",
    all_files,
    value = TRUE
)

cat("\n============================================\n")
cat("DOWNLOAD VALIDATION\n")
cat("============================================\n")

cat(
    "Filtered H5 matrices :",
    length(h5_files),
    "\n"
)

cat(
    "Selection CSV files  :",
    length(selection_files),
    "\n"
)

cat(
    "Position files       :",
    length(position_files),
    "\n"
)

cat("\nExpected approximately 16 of each file type.\n")

# ------------------------------------------------------------
# Save manifest
# ------------------------------------------------------------

manifest <- data.frame(
    file = basename(all_files),
    path = all_files,
    size_bytes = file.info(all_files)$size,
    stringsAsFactors = FALSE
)

write.csv(
    manifest,
    file.path(
        raw_dir,
        "GSE295206_local_manifest.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Reproducibility information
# ------------------------------------------------------------

sink(
    file.path(
        project_dir,
        "reproducibility",
        "00_download_sessionInfo.txt"
    )
)

sessionInfo()

sink()

cat("\n============================================\n")
cat("GSE295206 DOWNLOAD COMPLETE\n")
cat("============================================\n\n")

cat(
    "Data directory:\n",
    raw_dir,
    "\n\n"
)

