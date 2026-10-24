
options(stringsAsFactors = FALSE)
options(timeout = 1200)

# ============================================================
# Directories
# ============================================================

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
cat("SCRIPT 04c: AUDIT CURATED 1492-CELL XIN DATA\n")
cat("============================================\n")

# ============================================================
# Required Bioconductor infrastructure
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("Biobase", quietly = TRUE)) {
    BiocManager::install(
        "Biobase",
        ask = FALSE,
        update = FALSE
    )
}

suppressPackageStartupMessages({
    library(Biobase)
})

# ============================================================
# Download MuSiC Xin object
# ============================================================

outfile <- file.path(
    raw_dir,
    "XinT2Deset_MuSiC_1492cells.rds"
)

source_url <- paste0(
    "https://",
    "xuranw.github.io",
    "/MuSiC/data/XinT2Deset.rds"
)

if (!file.exists(outfile)) {

    cat("\nDownloading curated Xin 1492-cell object...\n")

    con <- gzcon(
        url(
            source_url,
            open = "rb"
        )
    )

    xin <- readRDS(con)

    close(con)

    saveRDS(
        xin,
        outfile,
        compress = FALSE
    )

} else {

    cat("\nUsing existing local curated object.\n")

    xin <- readRDS(outfile)
}

cat(
    "Local file:",
    outfile,
    "\n"
)

cat(
    "MD5:",
    unname(
        tools::md5sum(outfile)
    ),
    "\n"
)

# ============================================================
# Validate object class
# ============================================================

cat("\n============================================\n")
cat("OBJECT STRUCTURE\n")
cat("============================================\n")

cat(
    "Class:",
    paste(
        class(xin),
        collapse = ", "
    ),
    "\n"
)

if (!inherits(xin, "ExpressionSet")) {
    stop(
        "Downloaded object is not an ExpressionSet."
    )
}

expr <- Biobase::exprs(xin)

pd <- Biobase::pData(xin)

fd <- Biobase::fData(xin)

cat(
    "Expression dimensions:",
    nrow(expr),
    "genes x",
    ncol(expr),
    "cells\n"
)

cat(
    "Phenotype dimensions:",
    nrow(pd),
    "x",
    ncol(pd),
    "\n"
)

cat(
    "Feature metadata dimensions:",
    nrow(fd),
    "x",
    ncol(fd),
    "\n"
)

cat("\nPhenotype columns:\n")

print(
    names(pd)
)

cat("\nFeature metadata columns:\n")

print(
    names(fd)
)

# ============================================================
# Cell identifiers
# ============================================================

cat("\n============================================\n")
cat("CELL IDENTIFIERS\n")
cat("============================================\n")

cat("\nFirst 10 sample names:\n")
print(
    head(
        sampleNames(xin),
        10
    )
)

cat("\nLast 10 sample names:\n")
print(
    tail(
        sampleNames(xin),
        10
    )
)

# ============================================================
# Metadata values
# ============================================================

cat("\n============================================\n")
cat("PHENOTYPE METADATA PREVIEW\n")
cat("============================================\n")

print(
    head(
        pd,
        12
    )
)

for (z in names(pd)) {

    if (
        is.character(pd[[z]]) ||
        is.factor(pd[[z]])
    ) {

        vals <- unique(
            as.character(
                pd[[z]]
            )
        )

        if (length(vals) <= 30) {

            cat(
                "\nUnique values —",
                z,
                ":\n"
            )

            print(vals)
        }
    }
}

# ============================================================
# Identify likely metadata fields
# ============================================================

find_col <- function(candidates, available) {

    hit <- candidates[
        candidates %in%
        available
    ]

    if (!length(hit)) {
        return(NA_character_)
    }

    hit[1]
}

subject_col <- find_col(
    c(
        "SubjectName",
        "subject",
        "Subject",
        "donor",
        "donor_id"
    ),
    names(pd)
)

celltype_col <- find_col(
    c(
        "cellType",
        "CellType",
        "cell_type",
        "celltype"
    ),
    names(pd)
)

disease_col <- find_col(
    c(
        "Disease",
        "disease",
        "condition",
        "Condition"
    ),
    names(pd)
)

cat("\n============================================\n")
cat("DETECTED METADATA FIELDS\n")
cat("============================================\n")

cat(
    "Subject field:",
    subject_col,
    "\n"
)

cat(
    "Cell-type field:",
    celltype_col,
    "\n"
)

cat(
    "Disease field:",
    disease_col,
    "\n"
)

if (
    is.na(subject_col) ||
    is.na(celltype_col) ||
    is.na(disease_col)
) {
    stop(
        "Could not automatically identify essential metadata fields."
    )
}

# ============================================================
# Metadata summary
# ============================================================

subject <- as.character(
    pd[[subject_col]]
)

celltype <- as.character(
    pd[[celltype_col]]
)

disease <- as.character(
    pd[[disease_col]]
)

cat("\n============================================\n")
cat("CURATED DATASET SUMMARY\n")
cat("============================================\n")

cat(
    "Cells:",
    ncol(expr),
    "\n"
)

cat(
    "Genes:",
    nrow(expr),
    "\n"
)

cat(
    "Subjects:",
    length(
        unique(subject)
    ),
    "\n"
)

cat("\nCell types:\n")
print(
    table(celltype)
)

cat("\nDisease groups:\n")
print(
    table(disease)
)

cat("\nSubject x disease:\n")
print(
    table(
        subject,
        disease
    )
)

cat("\nSubject x cell type:\n")
print(
    table(
        subject,
        celltype
    )
)

# ============================================================
# Expression-value audit
# ============================================================

cat("\n============================================\n")
cat("EXPRESSION VALUE AUDIT\n")
cat("============================================\n")

finite_values <- expr[
    is.finite(expr)
]

cat(
    "Minimum:",
    min(finite_values),
    "\n"
)

cat(
    "Median:",
    median(finite_values),
    "\n"
)

cat(
    "Maximum:",
    max(finite_values),
    "\n"
)

cat(
    "Fraction zero:",
    mean(expr == 0),
    "\n"
)

sample_values <- finite_values[
    seq_len(
        min(
            length(finite_values),
            1000000L
        )
    )
]

fraction_integer_like <- mean(
    abs(
        sample_values -
        round(sample_values)
    ) < 1e-8
)

cat(
    "Fraction integer-like in first <=1M values:",
    fraction_integer_like,
    "\n"
)

# ============================================================
# Gene identifiers
# ============================================================

cat("\n============================================\n")
cat("GENE IDENTIFIER AUDIT\n")
cat("============================================\n")

cat("\nFirst 20 row names:\n")

print(
    head(
        rownames(expr),
        20
    )
)

cat("\nFirst 20 feature rows:\n")

print(
    head(
        fd,
        20
    )
)

# ============================================================
# Marker sanity check
#
# Try row names and feature metadata.
# ============================================================

markers <- c(
    "INS",
    "GCG",
    "SST",
    "PPY"
)

symbol_vector <- NULL

symbol_candidates <- c(
    "symbol",
    "Symbol",
    "gene_symbol",
    "GeneSymbol",
    "gene"
)

symbol_hit <- symbol_candidates[
    symbol_candidates %in%
    names(fd)
]

if (length(symbol_hit)) {

    symbol_vector <- as.character(
        fd[[symbol_hit[1]]]
    )

} else {

    symbol_vector <- rownames(expr)
}

marker_results <- list()

for (g in markers) {

    idx <- which(
        symbol_vector == g
    )

    if (!length(idx)) {

        marker_results[[g]] <- data.frame(
            gene = g,
            found = FALSE,
            top_celltype = NA_character_,
            stringsAsFactors = FALSE
        )

        next
    }

    idx <- idx[1]

    means <- tapply(
        expr[idx, ],
        celltype,
        mean,
        na.rm = TRUE
    )

    marker_results[[g]] <- data.frame(
        gene = g,
        found = TRUE,
        top_celltype =
            names(
                which.max(means)
            ),
        stringsAsFactors = FALSE
    )
}

marker_check <- do.call(
    rbind,
    marker_results
)

cat("\n============================================\n")
cat("MARKER SANITY CHECK\n")
cat("============================================\n")

print(
    marker_check,
    row.names = FALSE
)

# ============================================================
# Save clean audit metadata
# ============================================================

audit_meta <- data.frame(

    sample_name =
        sampleNames(xin),

    subject =
        subject,

    cell_type =
        celltype,

    disease =
        disease,

    stringsAsFactors = FALSE
)

write.csv(
    audit_meta,
    file.path(
        processed_dir,
        "GSE81608_MuSiC_1492_cell_metadata.csv"
    ),
    row.names = FALSE
)

write.csv(
    marker_check,
    file.path(
        islet_dir,
        "GSE81608_MuSiC_1492_marker_check.csv"
    ),
    row.names = FALSE
)

saveRDS(
    expr,
    file.path(
        processed_dir,
        "GSE81608_MuSiC_1492_expression.rds"
    ),
    compress = TRUE
)

# ============================================================
# Final validation
# ============================================================

dimension_pass <-
    ncol(expr) == 1492L &&
    nrow(expr) == 39849L

donor_pass <-
    length(
        unique(subject)
    ) == 18L

marker_expected <- c(
    INS = "beta",
    GCG = "alpha",
    SST = "delta",
    PPY = "gamma"
)

marker_pass <- TRUE

for (g in names(marker_expected)) {

    z <- marker_check[
        marker_check$gene == g,
        ,
        drop = FALSE
    ]

    if (
        nrow(z) != 1 ||
        !z$found ||
        tolower(z$top_celltype) !=
            tolower(marker_expected[g])
    ) {
        marker_pass <- FALSE
    }
}

cat("\n============================================\n")
cat("FINAL 1492-REFERENCE VALIDATION\n")
cat("============================================\n")

cat(
    "1492 cells + 39849 genes:",
    ifelse(
        dimension_pass,
        "PASS",
        "CHECK"
    ),
    "\n"
)

cat(
    "18 subjects:",
    ifelse(
        donor_pass,
        "PASS",
        "CHECK"
    ),
    "\n"
)

cat(
    "Endocrine marker identities:",
    ifelse(
        marker_pass,
        "PASS",
        "CHECK"
    ),
    "\n"
)

cat(
    "\nOVERALL:",
    ifelse(
        dimension_pass &&
        donor_pass &&
        marker_pass,
        "PASS",
        "CHECK"
    ),
    "\n"
)

# ============================================================
# Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "04c_GSE81608_1492_reference_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nSource:\n")
cat("MuSiC curated Xin GSE81608 ExpressionSet\n")

cat("\nLocal MD5:\n")
print(
    tools::md5sum(outfile)
)

cat("\nDimensions:\n")
print(dim(expr))

cat("\nMetadata columns:\n")
print(names(pd))

cat("\nCell types:\n")
print(table(celltype))

cat("\nDisease:\n")
print(table(disease))

cat("\nSubjects:\n")
print(table(subject))

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04c COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

