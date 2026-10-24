#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(edgeR)
})

cat("\n============================================\n")
cat("SCRIPT 05B: BETA-CELL DONOR ROBUSTNESS AUDIT\n")
cat("============================================\n")

dir.create("results/islet", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", recursive = TRUE, showWarnings = FALSE)

dge_file <- "results/islet/GSE81608_Beta_edgeR_Primary_DGEList.rds"

primary_file <-
    "results/islet/Table_P3_Beta_T2D_edgeR_Primary_18Donors.csv"

sensitivity_file <-
    "results/islet/Table_P3_Beta_T2D_Primary_vs_GE10Sensitivity.csv"

stopifnot(
    file.exists(dge_file),
    file.exists(primary_file),
    file.exists(sensitivity_file)
)

dge <- readRDS(dge_file)

primary <- read.csv(
    primary_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sensitivity <- read.csv(
    sensitivity_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!inherits(dge, "DGEList")) {
    stop("Saved primary object is not a DGEList.")
}

required_primary <- c(
    "gene", "logFC", "logCPM",
    "F", "PValue", "FDR"
)

if (!all(required_primary %in% colnames(primary))) {
    stop(
        "Primary table is missing required columns: ",
        paste(
            setdiff(required_primary, colnames(primary)),
            collapse = ";"
        )
    )
}

cat("Saved DGEList:", nrow(dge), "genes x",
    ncol(dge), "donors\n")

cat("Primary result rows:", nrow(primary), "\n")

# ------------------------------------------------------------
# 1. Reconstruct donor condition labels
# ------------------------------------------------------------

donor_id <- colnames(dge$counts)

normalize_condition <- function(x) {

    x2 <- tolower(trimws(as.character(x)))

    out <- rep(NA_character_, length(x2))

    out[
        grepl(
            "healthy|non[ _-]*t2d|non[ _-]*diabetic",
            x2
        )
    ] <- "Healthy"

    out[
        is.na(out) &
        grepl("(^|[ _-])t2d($|[ _-])|type[ _-]*2", x2)
    ] <- "T2D"

    out
}

condition <- rep(NA_character_, length(donor_id))

if ("condition" %in% colnames(dge$samples)) {

    tmp <- normalize_condition(dge$samples$condition)

    if (all(!is.na(tmp))) {
        condition <- tmp
    }
}

if (
    any(is.na(condition)) &&
    "group" %in% colnames(dge$samples)
) {

    tmp <- normalize_condition(dge$samples$group)

    if (all(!is.na(tmp))) {
        condition <- tmp
    }
}

# Fallback to donor identifiers.
if (any(is.na(condition))) {

    condition <- ifelse(
        grepl("^T2D", trimws(donor_id), ignore.case = TRUE),
        "T2D",
        "Healthy"
    )
}

condition <- factor(
    condition,
    levels = c("Healthy", "T2D")
)

sample_info <- data.frame(
    donor_id = donor_id,
    condition = condition,
    stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("DONOR STRUCTURE\n")
cat("============================================\n")

print(sample_info, row.names = FALSE)

cat("\nHealthy donors:",
    sum(condition == "Healthy"), "\n")

cat("T2D donors:",
    sum(condition == "T2D"), "\n")

if (
    sum(condition == "Healthy") != 12 ||
    sum(condition == "T2D") != 6
) {
    stop(
        "Expected 12 Healthy and 6 T2D donors."
    )
}

# ------------------------------------------------------------
# 2. Identify the primary 27 genes
# ------------------------------------------------------------

primary_sig <- primary[
    !is.na(primary$FDR) &
    primary$FDR < 0.05,
    ,
    drop = FALSE
]

primary_sig <- primary_sig[
    order(primary_sig$FDR),
    ,
    drop = FALSE
]

sig_genes <- primary_sig$gene

cat("\n============================================\n")
cat("PRIMARY RECEIVER SET\n")
cat("============================================\n")

cat("Primary FDR < 0.05 genes:",
    length(sig_genes), "\n")

if (length(sig_genes) != 27) {
    warning(
        "Expected 27 primary significant genes; found ",
        length(sig_genes)
    )
}

missing_from_dge <- setdiff(
    sig_genes,
    rownames(dge$counts)
)

if (length(missing_from_dge) > 0) {
    stop(
        "Primary genes missing from saved DGEList: ",
        paste(missing_from_dge, collapse = ";")
    )
}

# ------------------------------------------------------------
# 3. Donor-level normalized CPM
# ------------------------------------------------------------

cpm_matrix <- edgeR::cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = FALSE
)

cpm_sig <- cpm_matrix[
    sig_genes,
    ,
    drop = FALSE
]

donor_cpm_long <- do.call(
    rbind,
    lapply(
        sig_genes,
        function(g) {

            data.frame(
                gene = g,
                donor_id = donor_id,
                condition = as.character(condition),
                CPM = as.numeric(cpm_sig[g, ]),
                log2_CPM_plus_0.25 =
                    log2(as.numeric(cpm_sig[g, ]) + 0.25),
                stringsAsFactors = FALSE
            )
        }
    )
)

write.csv(
    donor_cpm_long,
    "results/islet/Table_P3_Beta_T2D_27genes_DonorLevel_CPM.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 4. Expression-support audit
#
# CPM >= 1 is used ONLY as a descriptive donor-support
# threshold. It does not determine differential expression.
# ------------------------------------------------------------

expression_support <- do.call(
    rbind,
    lapply(
        sig_genes,
        function(g) {

            pfc <- primary_sig$logFC[
                match(g, primary_sig$gene)
            ]

            h <- as.numeric(
                cpm_sig[g, condition == "Healthy"]
            )

            t <- as.numeric(
                cpm_sig[g, condition == "T2D"]
            )

            if (pfc > 0) {
                enriched_values <- t
                enriched_group <- "T2D"
            } else {
                enriched_values <- h
                enriched_group <- "Healthy"
            }

            total_enriched <- sum(enriched_values)

            dominance <- if (
                is.finite(total_enriched) &&
                total_enriched > 0
            ) {
                max(enriched_values) / total_enriched
            } else {
                NA_real_
            }

            data.frame(
                gene = g,

                primary_logFC = pfc,

                enriched_group = enriched_group,

                Healthy_n_donors =
                    length(h),

                T2D_n_donors =
                    length(t),

                Healthy_n_CPM_GE1 =
                    sum(h >= 1),

                T2D_n_CPM_GE1 =
                    sum(t >= 1),

                Healthy_fraction_CPM_GE1 =
                    mean(h >= 1),

                T2D_fraction_CPM_GE1 =
                    mean(t >= 1),

                Healthy_median_CPM =
                    median(h),

                T2D_median_CPM =
                    median(t),

                Healthy_max_CPM =
                    max(h),

                T2D_max_CPM =
                    max(t),

                enriched_group_n_CPM_GE1 =
                    sum(enriched_values >= 1),

                enriched_group_fraction_CPM_GE1 =
                    mean(enriched_values >= 1),

                enriched_group_max_CPM_share =
                    dominance,

                sparse_enriched_group =
                    mean(enriched_values >= 1) < 0.5,

                single_donor_dominance_GT50pct =
                    !is.na(dominance) &&
                    dominance > 0.50,

                stringsAsFactors = FALSE
            )
        }
    )
)

write.csv(
    expression_support,
    "results/islet/Table_P3_Beta_T2D_27genes_DonorExpressionSupport.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 5. edgeR fitting function
#
# This deliberately mirrors Script 05:
#   TMM
#   ~ condition
#   Healthy reference
#   robust estimateDisp
#   robust glmQLFit
#   T2D - Healthy
#
# The same primary-filtered gene universe is retained.
# ------------------------------------------------------------

run_edger <- function(y, cond) {

    cond <- factor(
        cond,
        levels = c("Healthy", "T2D")
    )

    design <- model.matrix(
        ~ condition,
        data = data.frame(
            condition = cond
        )
    )

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

    qlf <- edgeR::glmQLFTest(
        fit,
        coef = "conditionT2D"
    )

    tt <- edgeR::topTags(
        qlf,
        n = Inf,
        sort.by = "none",
        adjust.method = "BH"
    )$table

    tt <- data.frame(
        gene = rownames(tt),
        tt,
        row.names = NULL,
        check.names = FALSE
    )

    list(
        table = tt,
        design = design,
        fit = fit,
        qlf = qlf
    )
}

# ------------------------------------------------------------
# 6. Full-model reproducibility check
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FULL-MODEL REPRODUCIBILITY CHECK\n")
cat("============================================\n")

full_fit <- run_edger(
    dge,
    condition
)

full_check <- merge(
    primary[
        ,
        c(
            "gene",
            "logFC",
            "PValue",
            "FDR"
        )
    ],
    full_fit$table[
        ,
        c(
            "gene",
            "logFC",
            "PValue",
            "FDR"
        )
    ],
    by = "gene",
    suffixes = c("_script05", "_refit"),
    sort = FALSE
)

full_check$abs_logFC_difference <-
    abs(
        full_check$logFC_script05 -
        full_check$logFC_refit
    )

max_logfc_difference <-
    max(
        full_check$abs_logFC_difference,
        na.rm = TRUE
    )

cat(
    "Genes compared:",
    nrow(full_check),
    "\n"
)

cat(
    "Maximum absolute logFC difference:",
    format(max_logfc_difference, digits = 8),
    "\n"
)

write.csv(
    full_check,
    "results/islet/GSE81608_Beta_T2D_Script05_Refit_Check.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. Leave-one-donor-out analysis
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LEAVE-ONE-DONOR-OUT ANALYSIS\n")
cat("============================================\n")

lodo_list <- vector(
    "list",
    length(donor_id)
)

for (i in seq_along(donor_id)) {

    dropped_donor <- donor_id[i]

    cat(
        sprintf(
            "[%02d/%02d] Dropping: %s\n",
            i,
            length(donor_id),
            dropped_donor
        )
    )

    keep <- seq_along(donor_id) != i

    # Preserve original library sizes, but recompute
    # TMM normalization factors among retained donors.
    y_i <- dge[
        ,
        keep,
        keep.lib.sizes = TRUE
    ]

    cond_i <- condition[keep]

    if (
        sum(cond_i == "Healthy") < 2 ||
        sum(cond_i == "T2D") < 2
    ) {
        stop(
            "Insufficient biological replicates after dropping ",
            dropped_donor
        )
    }

    fit_i <- run_edger(
        y_i,
        cond_i
    )

    tab_i <- fit_i$table

    tab_i <- tab_i[
        match(sig_genes, tab_i$gene),
        ,
        drop = FALSE
    ]

    if (any(is.na(tab_i$gene))) {
        stop(
            "Failed to recover all primary genes in LODO fit: ",
            dropped_donor
        )
    }

    lodo_list[[i]] <- data.frame(
        dropped_donor = dropped_donor,
        dropped_condition =
            as.character(condition[i]),
        retained_Healthy =
            sum(cond_i == "Healthy"),
        retained_T2D =
            sum(cond_i == "T2D"),
        gene = tab_i$gene,
        logFC = tab_i$logFC,
        logCPM = tab_i$logCPM,
        F = tab_i$F,
        PValue = tab_i$PValue,
        FDR = tab_i$FDR,
        stringsAsFactors = FALSE
    )
}

lodo <- do.call(
    rbind,
    lodo_list
)

rownames(lodo) <- NULL

write.csv(
    lodo,
    "results/islet/Table_P3_Beta_T2D_27genes_LODO_AllFits.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 8. LODO robustness summary
# ------------------------------------------------------------

lodo_summary <- do.call(
    rbind,
    lapply(
        sig_genes,
        function(g) {

            z <- lodo[
                lodo$gene == g,
                ,
                drop = FALSE
            ]

            primary_fc <-
                primary_sig$logFC[
                    match(
                        g,
                        primary_sig$gene
                    )
                ]

            same_direction <-
                sign(z$logFC) ==
                sign(primary_fc)

            data.frame(
                gene = g,

                primary_logFC =
                    primary_fc,

                n_LODO =
                    nrow(z),

                n_LODO_same_direction =
                    sum(
                        same_direction,
                        na.rm = TRUE
                    ),

                fraction_LODO_same_direction =
                    mean(
                        same_direction,
                        na.rm = TRUE
                    ),

                same_direction_all_LODO =
                    all(
                        same_direction,
                        na.rm = TRUE
                    ),

                n_LODO_P_lt_0.05 =
                    sum(
                        z$PValue < 0.05,
                        na.rm = TRUE
                    ),

                fraction_LODO_P_lt_0.05 =
                    mean(
                        z$PValue < 0.05,
                        na.rm = TRUE
                    ),

                n_LODO_FDR_lt_0.05 =
                    sum(
                        z$FDR < 0.05,
                        na.rm = TRUE
                    ),

                fraction_LODO_FDR_lt_0.05 =
                    mean(
                        z$FDR < 0.05,
                        na.rm = TRUE
                    ),

                min_LODO_logFC =
                    min(
                        z$logFC,
                        na.rm = TRUE
                    ),

                max_LODO_logFC =
                    max(
                        z$logFC,
                        na.rm = TRUE
                    ),

                min_abs_LODO_logFC =
                    min(
                        abs(z$logFC),
                        na.rm = TRUE
                    ),

                median_LODO_logFC =
                    median(
                        z$logFC,
                        na.rm = TRUE
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

# ------------------------------------------------------------
# 9. Add >=10-cell sensitivity evidence
# ------------------------------------------------------------

sens_columns <- c(
    "gene",
    "sens_logFC",
    "sens_PValue",
    "sens_FDR",
    "same_direction_sensitivity",
    "primary_FDR05",
    "sensitivity_FDR05",
    "FDR05_both",
    "robust_primary_DE",
    "robust_high_effect"
)

available_sens_columns <-
    intersect(
        sens_columns,
        colnames(sensitivity)
    )

sens_sub <- sensitivity[
    ,
    available_sens_columns,
    drop = FALSE
]

# Ensure one row per gene.
sens_sub <- sens_sub[
    !duplicated(sens_sub$gene),
    ,
    drop = FALSE
]

robustness <- merge(
    primary_sig[
        ,
        c(
            "gene",
            "logFC",
            "logCPM",
            "PValue",
            "FDR"
        )
    ],
    expression_support,
    by = "gene",
    sort = FALSE
)

robustness <- merge(
    robustness,
    lodo_summary[
        ,
        setdiff(
            colnames(lodo_summary),
            "primary_logFC"
        )
    ],
    by = "gene",
    sort = FALSE
)

robustness <- merge(
    robustness,
    sens_sub,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

# Restore primary FDR ordering.
robustness <- robustness[
    match(
        sig_genes,
        robustness$gene
    ),
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# 10. Receiver robustness classification
#
# Primary receiver robustness is based on directional
# consistency, not on repeated FDR significance.
#
# Requiring FDR < 0.05 in every LODO fit would be an
# unnecessarily power-dependent criterion.
# ------------------------------------------------------------

if (
    "same_direction_sensitivity" %in%
    colnames(robustness)
) {

    robustness$receiver_primary_robust <-
        robustness$same_direction_all_LODO &
        robustness$same_direction_sensitivity

} else {

    robustness$receiver_primary_robust <-
        robustness$same_direction_all_LODO
}

if (
    "sensitivity_FDR05" %in%
    colnames(robustness)
) {

    robustness$receiver_high_confidence_core <-
        robustness$receiver_primary_robust &
        robustness$sensitivity_FDR05

} else {

    robustness$receiver_high_confidence_core <-
        FALSE
}

robustness$receiver_class <- ifelse(
    robustness$receiver_high_confidence_core,
    "Highest_confidence_core",
    ifelse(
        robustness$receiver_primary_robust,
        "Robust_primary_receiver",
        "Directionally_unstable"
    )
)

# Expression sparsity and donor dominance remain
# descriptive caution flags and DO NOT determine the
# receiver classification.

write.csv(
    robustness,
    "results/islet/Table_P3_Beta_T2D_ReceiverGene_Robustness.csv",
    row.names = FALSE
)

receiver_primary <- robustness[
    robustness$receiver_primary_robust,
    ,
    drop = FALSE
]

receiver_core <- robustness[
    robustness$receiver_high_confidence_core,
    ,
    drop = FALSE
]

write.csv(
    receiver_primary,
    "results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv",
    row.names = FALSE
)

write.csv(
    receiver_core,
    "results/islet/Table_P3_Beta_T2D_NicheNet_HighConfidenceCore.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 11. LODO logFC matrix
# ------------------------------------------------------------

lodo_logfc_matrix <- matrix(
    NA_real_,
    nrow = length(sig_genes),
    ncol = length(donor_id),
    dimnames = list(
        sig_genes,
        donor_id
    )
)

for (i in seq_along(donor_id)) {

    z <- lodo[
        lodo$dropped_donor == donor_id[i],
        ,
        drop = FALSE
    ]

    lodo_logfc_matrix[
        z$gene,
        donor_id[i]
    ] <- z$logFC
}

saveRDS(
    lodo_logfc_matrix,
    "results/islet/GSE81608_Beta_T2D_27genes_LODO_logFC_matrix.rds"
)

# ------------------------------------------------------------
# 12. Figures
# ------------------------------------------------------------

# Donor-level expression heatmap.
plot_matrix <- log2(
    cpm_sig + 0.25
)

png(
    "figures/GSE81608_Beta_T2D_27genes_DonorCPM_Heatmap.png",
    width = 2200,
    height = 2600,
    res = 220
)

stats::heatmap(
    plot_matrix,
    Rowv = NA,
    Colv = NA,
    scale = "row",
    margins = c(10, 10),
    cexRow = 0.75,
    cexCol = 0.70,
    main =
        "T2D-associated beta-cell genes: donor-level expression"
)

dev.off()

# LODO logFC heatmap.
png(
    "figures/GSE81608_Beta_T2D_27genes_LODO_LogFC_Heatmap.png",
    width = 2200,
    height = 2600,
    res = 220
)

stats::heatmap(
    lodo_logfc_matrix,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(10, 10),
    cexRow = 0.75,
    cexCol = 0.70,
    main =
        "Leave-one-donor-out T2D vs Healthy log2 fold changes"
)

dev.off()

# ------------------------------------------------------------
# 13. Method record
# ------------------------------------------------------------

method <- data.frame(
    dataset = "GSE81608_Xin_MuSiC_1492",
    cell_type = "beta",
    primary_donors = 18,
    healthy_donors = 12,
    T2D_donors = 6,
    primary_gene_universe = nrow(dge),
    primary_significant_genes =
        length(sig_genes),
    LODO_models =
        length(donor_id),
    normalization = "TMM",
    design = "~ condition",
    contrast = "T2D - Healthy",
    test =
        "edgeR quasi-likelihood F-test",
    robust = TRUE,
    donor_expression_threshold =
        "CPM >= 1; descriptive only",
    receiver_primary_definition =
        "Primary FDR<0.05 plus same effect direction in every LODO fit and >=10-cell donor sensitivity analysis",
    high_confidence_core_definition =
        "Primary receiver robustness plus FDR<0.05 in >=10-cell donor sensitivity analysis",
    interpretation =
        paste(
            "LODO assesses donor influence on the",
            "T2D-associated beta-cell transcriptional",
            "program. It does not establish causality",
            "or DRG-to-islet signaling."
        ),
    stringsAsFactors = FALSE
)

write.csv(
    method,
    "results/islet/GSE81608_Beta_T2D_DonorRobustness_Method.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 14. Console summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LODO ROBUSTNESS SUMMARY\n")
cat("============================================\n")

cat(
    "Primary FDR<0.05 genes:",
    length(sig_genes),
    "\n"
)

cat(
    "Same direction in all 18 LODO fits:",
    sum(
        robustness$same_direction_all_LODO,
        na.rm = TRUE
    ),
    "/",
    nrow(robustness),
    "\n"
)

cat(
    "FDR<0.05 in all 18 LODO fits:",
    sum(
        robustness$n_LODO_FDR_lt_0.05 == 18,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "FDR<0.05 in >=15/18 LODO fits:",
    sum(
        robustness$n_LODO_FDR_lt_0.05 >= 15,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Primary robust receiver genes:",
    sum(
        robustness$receiver_primary_robust,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Highest-confidence core:",
    sum(
        robustness$receiver_high_confidence_core,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Genes with <50% enriched-group donors at CPM>=1:",
    sum(
        robustness$sparse_enriched_group,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Genes with >50% enriched-group CPM from one donor:",
    sum(
        robustness$single_donor_dominance_GT50pct,
        na.rm = TRUE
    ),
    "\n"
)

cat("\n============================================\n")
cat("27-GENE ROBUSTNESS TABLE\n")
cat("============================================\n")

display_columns <- c(
    "gene",
    "logFC",
    "FDR",
    "enriched_group_fraction_CPM_GE1",
    "enriched_group_max_CPM_share",
    "n_LODO_same_direction",
    "n_LODO_P_lt_0.05",
    "n_LODO_FDR_lt_0.05",
    "sens_FDR",
    "receiver_class"
)

display_columns <- intersect(
    display_columns,
    colnames(robustness)
)

print(
    robustness[
        ,
        display_columns,
        drop = FALSE
    ],
    row.names = FALSE
)

cat("\n============================================\n")
cat("OUTPUT FILES\n")
cat("============================================\n")

cat(
    "results/islet/Table_P3_Beta_T2D_27genes_DonorLevel_CPM.csv\n"
)

cat(
    "results/islet/Table_P3_Beta_T2D_27genes_DonorExpressionSupport.csv\n"
)

cat(
    "results/islet/Table_P3_Beta_T2D_27genes_LODO_AllFits.csv\n"
)

cat(
    "results/islet/Table_P3_Beta_T2D_ReceiverGene_Robustness.csv\n"
)

cat(
    "results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv\n"
)

cat(
    "results/islet/Table_P3_Beta_T2D_NicheNet_HighConfidenceCore.csv\n"
)

cat(
    "figures/GSE81608_Beta_T2D_27genes_DonorCPM_Heatmap.png\n"
)

cat(
    "figures/GSE81608_Beta_T2D_27genes_LODO_LogFC_Heatmap.png\n"
)

capture.output(
    sessionInfo(),
    file =
        "reproducibility/05b_beta_T2D_donor_robustness_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 05B COMPLETED SUCCESSFULLY\n")
cat("============================================\n")
