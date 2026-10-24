
options(stringsAsFactors = FALSE)
options(timeout = 1200)

suppressPackageStartupMessages({
    library(GEOquery)
})

project_dir <- normalizePath(".", mustWork = TRUE)

raw_dir <- file.path(
    project_dir,
    "data",
    "raw",
    "GSE81608"
)

processed_dir <- file.path(
    project_dir,
    "data",
    "processed"
)

islet_dir <- file.path(
    project_dir,
    "results",
    "islet"
)

repro_dir <- file.path(
    project_dir,
    "reproducibility"
)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(islet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 04a: GSE81608 DOWNLOAD + VALIDATION\n")
cat("============================================\n")

# ============================================================
# 1. Download processed RPKM matrix
# ============================================================

rpkm_file <- file.path(
    raw_dir,
    "GSE81608_human_islets_rpkm.txt.gz"
)

if (!file.exists(rpkm_file)) {

    cat("\nDownloading GSE81608 processed RPKM matrix...\n")

    GEOquery::getGEOSuppFiles(
        "GSE81608",
        baseDir = raw_dir,
        makeDirectory = FALSE,
        fetch_files = TRUE,
        filter_regex = "human_islets_rpkm\\.txt\\.gz$"
    )
}

if (!file.exists(rpkm_file)) {

    candidates <- list.files(
        raw_dir,
        pattern = "human_islets_rpkm.*txt.*gz$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
    )

    if (length(candidates) != 1) {
        stop(
            "Could not uniquely identify GSE81608 RPKM file."
        )
    }

    rpkm_file <- candidates[1]
}

cat(
    "RPKM file:\n",
    rpkm_file,
    "\n"
)

cat(
    "File size:",
    round(file.info(rpkm_file)$size / 1024^2, 2),
    "MB\n"
)

# ============================================================
# 2. Download GEO Series Matrix metadata
# ============================================================

cat("\nDownloading GEO Series Matrix metadata...\n")

gse <- getGEO(
    "GSE81608",
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = raw_dir
)

if (length(gse) < 1) {
    stop("GSE81608 GEO object could not be retrieved.")
}

eset <- gse[[1]]

meta_raw <- pData(eset)

cat(
    "GEO metadata samples:",
    nrow(meta_raw),
    "\n"
)

# ============================================================
# 3. Parse GEO characteristics
# ============================================================

characteristic_cols <- grep(
    "^characteristics_ch1",
    names(meta_raw),
    value = TRUE
)

if (!length(characteristic_cols)) {
    stop("No GEO characteristics columns found.")
}

extract_characteristic <- function(df, prefix) {

    out <- rep(
        NA_character_,
        nrow(df)
    )

    for (cc in characteristic_cols) {

        z <- as.character(
            df[[cc]]
        )

        hit <- grepl(
            paste0(
                "^",
                prefix,
                "\\s*:"
            ),
            z,
            ignore.case = TRUE
        )

        out[hit] <- trimws(
            sub(
                paste0(
                    "^",
                    prefix,
                    "\\s*:"
                ),
                "",
                z[hit],
                ignore.case = TRUE
            )
        )
    }

    out
}

cell_meta <- data.frame(

    GSM = rownames(meta_raw),

    title = as.character(
        meta_raw$title
    ),

    donor_id =
        extract_characteristic(
            meta_raw,
            "donor id"
        ),

    condition =
        extract_characteristic(
            meta_raw,
            "condition"
        ),

    age =
        extract_characteristic(
            meta_raw,
            "age"
        ),

    ethnicity =
        extract_characteristic(
            meta_raw,
            "ethnicity"
        ),

    gender =
        extract_characteristic(
            meta_raw,
            "gender"
        ),

    cell_subtype =
        extract_characteristic(
            meta_raw,
            "cell subtype"
        ),

    stringsAsFactors = FALSE
)

# ============================================================
# 4. Standardize metadata conservatively
# ============================================================

cell_meta$condition_original <-
    cell_meta$condition

cell_meta$cell_subtype_original <-
    cell_meta$cell_subtype

cell_meta$condition <-
    trimws(
        cell_meta$condition
    )

cell_meta$cell_subtype <-
    tolower(
        trimws(
            cell_meta$cell_subtype
        )
    )

cell_meta$cell_subtype[
    cell_meta$cell_subtype %in%
        c(
            "alpha",
            "α"
        )
] <- "alpha"

cell_meta$cell_subtype[
    cell_meta$cell_subtype %in%
        c(
            "beta",
            "β"
        )
] <- "beta"

cell_meta$cell_subtype[
    cell_meta$cell_subtype %in%
        c(
            "delta",
            "δ"
        )
] <- "delta"

cell_meta$cell_subtype[
    cell_meta$cell_subtype %in%
        c(
            "pp",
            "gamma",
            "γ"
        )
] <- "PP"

# ============================================================
# 5. Read RPKM matrix
# ============================================================

cat("\nReading RPKM matrix...\n")

rpkm_raw <- read.delim(
    gzfile(rpkm_file),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat(
    "Raw RPKM dimensions:",
    nrow(rpkm_raw),
    "rows x",
    ncol(rpkm_raw),
    "columns\n"
)

cat("\nFirst 10 column names:\n")
print(
    head(
        names(rpkm_raw),
        10
    )
)

cat("\nFirst 10 rows / first 5 columns:\n")
print(
    head(
        rpkm_raw[
            ,
            seq_len(
                min(
                    5,
                    ncol(rpkm_raw)
                )
            ),
            drop = FALSE
        ],
        10
    )
)

# ============================================================
# 6. Identify expression columns using GSM identifiers
# ============================================================

gsm_cols <- intersect(
    names(rpkm_raw),
    cell_meta$GSM
)

cat(
    "\nRPKM columns directly matching GSM IDs:",
    length(gsm_cols),
    "\n"
)

# Some GEO processed matrices may use sample numbers/titles
# rather than GSM IDs. We inspect before making assumptions.

title_match <- intersect(
    names(rpkm_raw),
    cell_meta$title
)

cat(
    "RPKM columns directly matching GEO titles:",
    length(title_match),
    "\n"
)

# ============================================================
# 7. Metadata validation
# ============================================================

cat("\n============================================\n")
cat("METADATA VALIDATION\n")
cat("============================================\n")

cat(
    "Total GEO samples:",
    nrow(cell_meta),
    "\n"
)

cat(
    "Unique donors:",
    length(
        unique(
            na.omit(
                cell_meta$donor_id
            )
        )
    ),
    "\n"
)

cat("\nCondition counts:\n")
print(
    table(
        cell_meta$condition,
        useNA = "ifany"
    )
)

cat("\nCell subtype counts:\n")
print(
    table(
        cell_meta$cell_subtype,
        useNA = "ifany"
    )
)

cat("\nDonor x condition:\n")
print(
    table(
        cell_meta$donor_id,
        cell_meta$condition,
        useNA = "ifany"
    )
)

cat("\nDonor x cell subtype:\n")
print(
    table(
        cell_meta$donor_id,
        cell_meta$cell_subtype,
        useNA = "ifany"
    )
)

# ============================================================
# 8. Endocrine-cell subset
# ============================================================

endocrine_types <- c(
    "alpha",
    "beta",
    "delta",
    "PP"
)

cell_meta$is_endocrine <-
    cell_meta$cell_subtype %in%
    endocrine_types

cat(
    "\nEndocrine cells:",
    sum(
        cell_meta$is_endocrine,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Alpha cells:",
    sum(
        cell_meta$cell_subtype ==
            "alpha",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Beta cells:",
    sum(
        cell_meta$cell_subtype ==
            "beta",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Delta cells:",
    sum(
        cell_meta$cell_subtype ==
            "delta",
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "PP cells:",
    sum(
        cell_meta$cell_subtype ==
            "PP",
        na.rm = TRUE
    ),
    "\n"
)

# ============================================================
# 9. Missing metadata check
# ============================================================

missing_summary <- data.frame(

    field = c(
        "donor_id",
        "condition",
        "cell_subtype"
    ),

    missing = c(
        sum(
            is.na(cell_meta$donor_id) |
            !nzchar(cell_meta$donor_id)
        ),
        sum(
            is.na(cell_meta$condition) |
            !nzchar(cell_meta$condition)
        ),
        sum(
            is.na(cell_meta$cell_subtype) |
            !nzchar(cell_meta$cell_subtype)
        )
    ),

    stringsAsFactors = FALSE
)

cat("\nMissing metadata:\n")

print(
    missing_summary,
    row.names = FALSE
)

# ============================================================
# 10. Save frozen metadata
# ============================================================

write.csv(
    cell_meta,
    file.path(
        processed_dir,
        "GSE81608_cell_metadata.csv"
    ),
    row.names = FALSE
)

write.csv(
    missing_summary,
    file.path(
        islet_dir,
        "GSE81608_metadata_missingness.csv"
    ),
    row.names = FALSE
)

donor_summary <- as.data.frame(
    table(
        cell_meta$donor_id,
        cell_meta$condition,
        cell_meta$cell_subtype
    ),
    stringsAsFactors = FALSE
)

names(donor_summary) <- c(
    "donor_id",
    "condition",
    "cell_subtype",
    "n_cells"
)

donor_summary <- donor_summary[
    donor_summary$n_cells > 0,
    ,
    drop = FALSE
]

write.csv(
    donor_summary,
    file.path(
        islet_dir,
        "GSE81608_donor_condition_celltype_counts.csv"
    ),
    row.names = FALSE
)

# Save raw matrix temporarily exactly as read.
# We will determine its identifier structure before converting
# it into the receptor-analysis matrix.

saveRDS(
    rpkm_raw,
    file.path(
        processed_dir,
        "GSE81608_RPKM_raw_import.rds"
    ),
    compress = TRUE
)

# ============================================================
# 11. Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "04a_GSE81608_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nRPKM file:\n")
print(rpkm_file)

cat("\nRPKM dimensions:\n")
print(dim(rpkm_raw))

cat("\nMetadata dimensions:\n")
print(dim(cell_meta))

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04a COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

