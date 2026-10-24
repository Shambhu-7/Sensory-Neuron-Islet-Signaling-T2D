
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(Biobase)
    library(edgeR)
})

xin_file <- "data/raw/GSE81608/XinT2Deset_MuSiC_1492cells.rds"

outdir <- "results/islet"
reprodir <- "reproducibility"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reprodir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 05: DONOR-LEVEL BETA-CELL T2D RESPONSE\n")
cat("============================================\n")

# ------------------------------------------------------------
# 1. Load audited 1492-cell ExpressionSet
# ------------------------------------------------------------

xin <- readRDS(xin_file)

if (!inherits(xin, "ExpressionSet")) {
    stop("Input is not an ExpressionSet.")
}

expr <- Biobase::exprs(xin)
meta <- Biobase::pData(xin)

if (!identical(colnames(expr), rownames(meta))) {
    meta <- meta[
        match(colnames(expr), rownames(meta)),
        ,
        drop = FALSE
    ]
}

if (anyNA(meta$SubjectName) ||
    anyNA(meta$cellType) ||
    anyNA(meta$Disease)) {
    stop("Missing required phenotype metadata.")
}

cat(
    "Input:",
    nrow(expr), "genes x",
    ncol(expr), "cells\n"
)

# ------------------------------------------------------------
# 2. Assay audit
# ------------------------------------------------------------

if (anyNA(expr) || any(expr < 0)) {
    stop("Expression matrix contains NA or negative values.")
}

nz <- expr[expr != 0]

integer_fraction <-
    mean(
        abs(nz - round(nz)) < 1e-8
    )

cat(
    "Fraction non-zero values integer-like:",
    integer_fraction,
    "\n"
)

cat(
    "Fractional non-zero observations:",
    sum(abs(nz - round(nz)) >= 1e-8),
    "\n"
)

if (integer_fraction < 0.99) {
    stop(
        "Assay is insufficiently count-like for this pipeline."
    )
}

# Fractional observations are preserved, not rounded.

# ------------------------------------------------------------
# 3. Select beta cells
# ------------------------------------------------------------

beta_idx <-
    tolower(
        trimws(
            as.character(meta$cellType)
        )
    ) == "beta"

beta_expr <- expr[
    ,
    beta_idx,
    drop = FALSE
]

beta_meta <- meta[
    beta_idx,
    ,
    drop = FALSE
]

cat(
    "Beta cells:",
    ncol(beta_expr),
    "\n"
)

if (ncol(beta_expr) != 472L) {
    stop("Expected 472 beta cells.")
}

# ------------------------------------------------------------
# 4. Standardize disease labels
# ------------------------------------------------------------

disease <- as.character(
    beta_meta$Disease
)

condition <- ifelse(
    tolower(disease) %in%
        c("healthy", "non-diabetic", "non diabetic"),
    "Healthy",
    ifelse(
        toupper(disease) == "T2D",
        "T2D",
        NA_character_
    )
)

if (anyNA(condition)) {
    stop("Unrecognized disease label.")
}

beta_meta$condition <- condition
beta_meta$donor_id <- as.character(
    beta_meta$SubjectName
)

# Each donor must belong to exactly one condition.

donor_condition_check <- tapply(
    beta_meta$condition,
    beta_meta$donor_id,
    function(z) length(unique(z))
)

if (any(donor_condition_check != 1L)) {
    stop(
        "At least one donor occurs in multiple disease groups."
    )
}

# ------------------------------------------------------------
# 5. Aggregate beta-cell counts by donor
# ------------------------------------------------------------

donors <- unique(
    beta_meta$donor_id
)

pb_list <- lapply(
    donors,
    function(d) {

        ii <- which(
            beta_meta$donor_id == d
        )

        rowSums(
            beta_expr[
                ,
                ii,
                drop = FALSE
            ]
        )
    }
)

pb <- do.call(
    cbind,
    pb_list
)

rownames(pb) <- rownames(beta_expr)
colnames(pb) <- donors

sample_info <- data.frame(
    donor_id = donors,
    condition = vapply(
        donors,
        function(d) {
            unique(
                beta_meta$condition[
                    beta_meta$donor_id == d
                ]
            )
        },
        character(1)
    ),
    n_beta_cells = vapply(
        donors,
        function(d) {
            sum(
                beta_meta$donor_id == d
            )
        },
        integer(1)
    ),
    stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$donor_id

# ------------------------------------------------------------
# 6. Collapse duplicate gene symbols if present
# ------------------------------------------------------------

if (
    is.null(rownames(pb)) ||
    any(!nzchar(rownames(pb)))
) {
    stop("Gene identifiers are missing.")
}

cat(
    "Duplicated gene identifiers before collapse:",
    sum(duplicated(rownames(pb))),
    "\n"
)

if (anyDuplicated(rownames(pb))) {

    pb <- rowsum(
        pb,
        group = rownames(pb),
        reorder = FALSE
    )
}

markers <- c(
    "INS",
    "GCG",
    "SST",
    "PPY"
)

cat(
    "Canonical endocrine symbols present:",
    paste(
        markers[
            markers %in% rownames(pb)
        ],
        collapse = ";"
    ),
    "\n"
)

if (!"INS" %in% rownames(pb)) {
    stop("Gene identifiers do not appear to be symbols.")
}

# ------------------------------------------------------------
# 7. Save pseudobulk input
# ------------------------------------------------------------

saveRDS(
    pb,
    file.path(
        outdir,
        "GSE81608_1492_Beta_Donor_Pseudobulk_Counts.rds"
    )
)

write.csv(
    sample_info,
    file.path(
        outdir,
        "GSE81608_1492_Beta_Donor_Pseudobulk_Samples.csv"
    ),
    row.names = FALSE
)

cat("\n============================================\n")
cat("PSEUDOBULK SAMPLES\n")
cat("============================================\n")

print(
    sample_info,
    row.names = FALSE
)

cat(
    "\nDonors:",
    nrow(sample_info),
    "\n"
)

cat(
    "Healthy:",
    sum(sample_info$condition == "Healthy"),
    "\n"
)

cat(
    "T2D:",
    sum(sample_info$condition == "T2D"),
    "\n"
)

# ------------------------------------------------------------
# 8. edgeR analysis function
# ------------------------------------------------------------

run_edger <- function(
    counts,
    samples,
    label
) {

    samples <- samples[
        colnames(counts),
        ,
        drop = FALSE
    ]

    samples$condition <- factor(
        samples$condition,
        levels = c(
            "Healthy",
            "T2D"
        )
    )

    design <- model.matrix(
        ~ condition,
        data = samples
    )

    rownames(design) <- rownames(samples)

    y <- edgeR::DGEList(
        counts = counts,
        samples = samples
    )

    keep <- edgeR::filterByExpr(
        y,
        design = design
    )

    cat(
        "\n",
        label,
        "genes before filter:",
        nrow(y),
        "\n"
    )

    cat(
        label,
        "genes retained:",
        sum(keep),
        "\n"
    )

    y <- y[
        keep,
        ,
        keep.lib.sizes = FALSE
    ]

    y <- edgeR::calcNormFactors(
        y,
        method = "TMM"
    )

    y <- edgeR::estimateDisp(
        y,
        design,
        robust = TRUE
    )

    fit <- edgeR::glmQLFit(
        y,
        design,
        robust = TRUE
    )

    coef_name <- "conditionT2D"

    if (!coef_name %in% colnames(design)) {
        stop(
            "T2D coefficient not found."
        )
    }

    qlf <- edgeR::glmQLFTest(
        fit,
        coef = coef_name
    )

    tab <- edgeR::topTags(
        qlf,
        n = Inf,
        sort.by = "none"
    )$table

    tab$gene <- rownames(tab)

    tab <- tab[
        ,
        c(
            "gene",
            setdiff(
                names(tab),
                "gene"
            )
        ),
        drop = FALSE
    ]

    tab$direction <- ifelse(
        tab$logFC > 0,
        "T2D_up",
        ifelse(
            tab$logFC < 0,
            "T2D_down",
            "zero"
        )
    )

    tab$FDR05 <-
        tab$FDR < 0.05

    tab$FDR05_absLogFC025 <-
        tab$FDR < 0.05 &
        abs(tab$logFC) >= 0.25

    logcpm <- edgeR::cpm(
        y,
        log = TRUE,
        prior.count = 2
    )

    list(
        table = tab,
        y = y,
        design = design,
        fit = fit,
        qlf = qlf,
        logCPM = logcpm
    )
}

# ------------------------------------------------------------
# 9. PRIMARY: all 18 donors
# ------------------------------------------------------------

primary <- run_edger(
    counts = pb,
    samples = sample_info,
    label = "PRIMARY_18_DONORS"
)

write.csv(
    primary$table,
    file.path(
        outdir,
        "Table_P3_Beta_T2D_edgeR_Primary_18Donors.csv"
    ),
    row.names = FALSE
)

saveRDS(
    primary$y,
    file.path(
        outdir,
        "GSE81608_Beta_edgeR_Primary_DGEList.rds"
    )
)

write.csv(
    primary$logCPM,
    file.path(
        outdir,
        "GSE81608_Beta_Primary_TMM_logCPM.csv"
    )
)

# ------------------------------------------------------------
# 10. SENSITIVITY: donors with >=10 beta cells
# ------------------------------------------------------------

sens_donors <- sample_info$donor_id[
    sample_info$n_beta_cells >= 10L
]

pb_sens <- pb[
    ,
    sens_donors,
    drop = FALSE
]

sample_sens <- sample_info[
    sens_donors,
    ,
    drop = FALSE
]

cat("\n============================================\n")
cat(">=10 BETA-CELL DONOR SENSITIVITY\n")
cat("============================================\n")

cat(
    "Donors retained:",
    nrow(sample_sens),
    "\n"
)

cat(
    "Healthy:",
    sum(sample_sens$condition == "Healthy"),
    "\n"
)

cat(
    "T2D:",
    sum(sample_sens$condition == "T2D"),
    "\n"
)

cat(
    "Excluded donors:",
    paste(
        setdiff(
            sample_info$donor_id,
            sens_donors
        ),
        collapse = ";"
    ),
    "\n"
)

sensitivity <- run_edger(
    counts = pb_sens,
    samples = sample_sens,
    label = "SENSITIVITY_GE10_CELLS"
)

write.csv(
    sensitivity$table,
    file.path(
        outdir,
        "Table_P3_Beta_T2D_edgeR_Sensitivity_GE10cells.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 11. Compare primary and sensitivity
# ------------------------------------------------------------

p <- primary$table
s <- sensitivity$table

names(s)[
    names(s) != "gene"
] <- paste0(
    "sens_",
    names(s)[
        names(s) != "gene"
    ]
)

comparison <- merge(
    p,
    s,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

comparison$same_direction_sensitivity <-
    !is.na(comparison$sens_logFC) &
    sign(comparison$logFC) ==
    sign(comparison$sens_logFC)

comparison$primary_FDR05 <-
    comparison$FDR < 0.05

comparison$sensitivity_FDR05 <-
    !is.na(comparison$sens_FDR) &
    comparison$sens_FDR < 0.05

comparison$FDR05_both <-
    comparison$primary_FDR05 &
    comparison$sensitivity_FDR05

comparison$robust_primary_DE <-
    comparison$primary_FDR05 &
    comparison$same_direction_sensitivity

comparison$robust_high_effect <-
    comparison$primary_FDR05 &
    abs(comparison$logFC) >= 0.25 &
    comparison$same_direction_sensitivity

write.csv(
    comparison,
    file.path(
        outdir,
        "Table_P3_Beta_T2D_Primary_vs_GE10Sensitivity.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PRIMARY T2D BETA-CELL RESPONSE\n")
cat("============================================\n")

cat(
    "Tested genes:",
    nrow(primary$table),
    "\n"
)

cat(
    "FDR < 0.05:",
    sum(primary$table$FDR05),
    "\n"
)

cat(
    "T2D-up FDR < 0.05:",
    sum(
        primary$table$FDR05 &
        primary$table$logFC > 0
    ),
    "\n"
)

cat(
    "T2D-down FDR < 0.05:",
    sum(
        primary$table$FDR05 &
        primary$table$logFC < 0
    ),
    "\n"
)

cat(
    "FDR < 0.05 & |logFC| >= 0.25:",
    sum(
        primary$table$FDR05_absLogFC025
    ),
    "\n"
)

cat("\n============================================\n")
cat("SENSITIVITY CONCORDANCE\n")
cat("============================================\n")

sig <- comparison$primary_FDR05

cat(
    "Primary FDR<0.05 genes with same direction:",
    sum(
        comparison$same_direction_sensitivity[
            sig
        ],
        na.rm = TRUE
    ),
    "/",
    sum(sig),
    "\n"
)

cat(
    "Primary FDR<0.05 genes also FDR<0.05 in sensitivity:",
    sum(
        comparison$FDR05_both,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Robust primary DE (FDR<0.05 + same sensitivity direction):",
    sum(
        comparison$robust_primary_DE,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Robust high-effect genes:",
    sum(
        comparison$robust_high_effect,
        na.rm = TRUE
    ),
    "\n"
)

# ------------------------------------------------------------
# 13. Top genes
# ------------------------------------------------------------

cat("\n============================================\n")
cat("TOP T2D-ASSOCIATED BETA-CELL GENES\n")
cat("Positive logFC = higher in T2D\n")
cat("============================================\n")

top <- primary$table[
    order(
        primary$table$FDR,
        -abs(primary$table$logFC)
    ),
    ,
    drop = FALSE
]

print(
    head(
        top[
            ,
            c(
                "gene",
                "logFC",
                "logCPM",
                "F",
                "PValue",
                "FDR"
            )
        ],
        30
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 14. Library/sample QC
# ------------------------------------------------------------

sample_qc <- primary$y$samples

sample_qc$donor_id <- rownames(
    primary$y$samples
)

sample_qc$n_beta_cells <-
    sample_info[
        sample_qc$donor_id,
        "n_beta_cells"
    ]

write.csv(
    sample_qc,
    file.path(
        outdir,
        "GSE81608_Beta_Pseudobulk_Library_QC.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 15. MDS
# ------------------------------------------------------------

png(
    "figures/GSE81608_Beta_Pseudobulk_MDS.png",
    width = 1800,
    height = 1500,
    res = 200
)

limma::plotMDS(
    primary$y,
    labels = sample_info[
        colnames(primary$y),
        "condition"
    ]
)

dev.off()

# ------------------------------------------------------------
# 16. Freeze method
# ------------------------------------------------------------

method <- data.frame(
    dataset = "GSE81608_Xin_MuSiC_1492",
    cell_type = "beta",
    beta_cells = ncol(beta_expr),
    biological_replicates = nrow(sample_info),
    healthy_donors =
        sum(sample_info$condition == "Healthy"),
    T2D_donors =
        sum(sample_info$condition == "T2D"),
    aggregation = "sum counts within donor",
    primary_design = "~ condition",
    contrast = "T2D - Healthy",
    normalization = "TMM",
    test = "edgeR quasi-likelihood F-test",
    robust = TRUE,
    fractional_counts = "preserved; not rounded",
    sensitivity =
        "repeat after requiring >=10 beta cells per donor",
    interpretation =
        paste(
            "T2D-associated beta-cell transcriptional",
            "differences; association does not imply",
            "causation by DRG-derived ligands."
        ),
    stringsAsFactors = FALSE
)

write.csv(
    method,
    file.path(
        outdir,
        "GSE81608_Beta_T2D_edgeR_Method.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 17. Reproducibility
# ------------------------------------------------------------

sink(
    file.path(
        reprodir,
        "05_beta_T2D_pseudobulk_edgeR_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nMethod:\n")
print(method)

cat("\nPrimary samples:\n")
print(sample_info)

cat("\nSensitivity samples:\n")
print(sample_sens)

cat("\nPrimary design:\n")
print(primary$design)

cat("\nPrimary DE summary:\n")
print(
    summary(
        limma::decideTests(
            primary$qlf,
            adjust.method = "BH",
            p.value = 0.05
        )
    )
)

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 05 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

