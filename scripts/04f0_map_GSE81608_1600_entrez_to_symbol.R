
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
})

infile <- "data/processed/GSE81608_RPKM_QC_1492cells.rds"

outfile <- "data/processed/GSE81608_RPKM_1600cells_SYMBOL.rds"

mapfile <- "results/islet/GSE81608_1600_Entrez_to_SYMBOL_mapping.csv"

sessionfile <- paste0(
    "reproducibility/",
    "04f0_GSE81608_1600_entrez_to_symbol_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 04f0: ENTREZ -> HGNC SYMBOL MAPPING\n")
cat("============================================\n")

expr <- readRDS(infile)

if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
}

cat(
    "Input dimensions:",
    nrow(expr), "genes x",
    ncol(expr), "cells\n"
)

entrez <- as.character(
    rownames(expr)
)

if (
    is.null(entrez) ||
    any(!nzchar(entrez))
) {
    stop("Input matrix has invalid row identifiers.")
}

if (!all(grepl("^[0-9]+$", entrez))) {
    stop(
        "Input row names are not uniformly numeric Entrez IDs."
    )
}

cat(
    "Numeric Entrez identifiers: PASS\n"
)

# ------------------------------------------------------------
# Entrez -> official HGNC symbol
# ------------------------------------------------------------

symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = entrez,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)

symbols <- as.character(symbols)

mapping <- data.frame(
    ENTREZID = entrez,
    SYMBOL = symbols,
    stringsAsFactors = FALSE
)

mapping$mapped <-
    !is.na(mapping$SYMBOL) &
    nzchar(mapping$SYMBOL)

write.csv(
    mapping,
    mapfile,
    row.names = FALSE
)

cat("\n============================================\n")
cat("MAPPING SUMMARY\n")
cat("============================================\n")

cat(
    "Input Entrez IDs:",
    nrow(mapping),
    "\n"
)

cat(
    "Mapped to HGNC symbol:",
    sum(mapping$mapped),
    "\n"
)

cat(
    "Unmapped:",
    sum(!mapping$mapped),
    "\n"
)

cat(
    "Mapping fraction:",
    round(mean(mapping$mapped), 4),
    "\n"
)

mapped_symbols <- mapping$SYMBOL[
    mapping$mapped
]

cat(
    "Unique mapped symbols:",
    length(unique(mapped_symbols)),
    "\n"
)

cat(
    "Rows involved in duplicated symbols:",
    sum(duplicated(mapped_symbols)),
    "\n"
)

# ------------------------------------------------------------
# Retain mapped genes
# ------------------------------------------------------------

keep <- mapping$mapped

expr_mapped <- expr[
    keep,
    ,
    drop = FALSE
]

symbols_mapped <- mapping$SYMBOL[
    keep
]

# ------------------------------------------------------------
# Collapse multiple Entrez IDs mapping to same symbol
#
# This sensitivity analysis uses only expression > 0.
# Summing duplicated-symbol RPKM rows therefore preserves
# whether a gene symbol is detected in each cell.
# ------------------------------------------------------------

expr_symbol <- rowsum(
    expr_mapped,
    group = symbols_mapped,
    reorder = FALSE
)

cat("\n============================================\n")
cat("SYMBOL MATRIX\n")
cat("============================================\n")

cat(
    "Dimensions:",
    nrow(expr_symbol), "symbols x",
    ncol(expr_symbol), "cells\n"
)

cat(
    "Duplicated row names after collapse:",
    sum(duplicated(rownames(expr_symbol))),
    "\n"
)

# ------------------------------------------------------------
# Endocrine marker validation
# ------------------------------------------------------------

markers <- c(
    "INS",
    "GCG",
    "SST",
    "PPY"
)

marker_present <- markers %in%
    rownames(expr_symbol)

marker_check <- data.frame(
    gene = markers,
    present = marker_present,
    stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("ENDOCRINE MARKER CHECK\n")
cat("============================================\n")

print(
    marker_check,
    row.names = FALSE
)

if (!all(marker_present)) {
    stop(
        "One or more canonical endocrine markers failed symbol mapping."
    )
}

# ------------------------------------------------------------
# High-confidence receptor overlap
# ------------------------------------------------------------

hc_file <- paste0(
    "results/islet/",
    "Table_P3_HighConfidence_LR_Islet_Compatibility.csv"
)

hc <- read.csv(
    hc_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

candidate_receptors <- sort(
    unique(hc$receptor)
)

present_receptors <- intersect(
    candidate_receptors,
    rownames(expr_symbol)
)

absent_receptors <- setdiff(
    candidate_receptors,
    rownames(expr_symbol)
)

cat("\n============================================\n")
cat("HIGH-CONFIDENCE RECEPTOR SYMBOL CHECK\n")
cat("============================================\n")

cat(
    "Candidate HC receptors:",
    length(candidate_receptors),
    "\n"
)

cat(
    "Present after Entrez->SYMBOL mapping:",
    length(present_receptors),
    "\n"
)

cat(
    "Absent:",
    length(absent_receptors),
    "\n"
)

if (length(absent_receptors)) {
    cat(
        "\nAbsent receptors:\n",
        paste(
            absent_receptors,
            collapse = "; "
        ),
        "\n"
    )
}

# ------------------------------------------------------------
# Save derivative matrix
# ------------------------------------------------------------

saveRDS(
    expr_symbol,
    outfile,
    compress = TRUE
)

cat("\nSaved symbol matrix:\n")
cat(outfile, "\n")

cat("\nSaved mapping audit:\n")
cat(mapfile, "\n")

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

sink(sessionfile)

cat("Analysis time:\n")
print(Sys.time())

cat("\nInput file:\n")
cat(infile, "\n")

cat("\nOutput file:\n")
cat(outfile, "\n")

cat("\nOriginal dimensions:\n")
print(dim(expr))

cat("\nSymbol dimensions:\n")
print(dim(expr_symbol))

cat("\nMapping summary:\n")
print(
    c(
        input_entrez = length(entrez),
        mapped = sum(mapping$mapped),
        unmapped = sum(!mapping$mapped),
        unique_symbols = nrow(expr_symbol)
    )
)

cat("\nAnnotation package versions:\n")

cat(
    "AnnotationDbi:",
    as.character(packageVersion("AnnotationDbi")),
    "\n"
)

cat(
    "org.Hs.eg.db:",
    as.character(packageVersion("org.Hs.eg.db")),
    "\n"
)

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04f0 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

