#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(nichenetr)
})

cat(strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16a — GSE221156 INDEPENDENT NICHE-NET REPLICATION\n")
cat(strrep("=", 78), "\n", sep = "")

cat("nichenetr version:", as.character(packageVersion("nichenetr")), "\n")

# =============================================================================
# 1. Frozen inputs
# =============================================================================

de_file <- paste0(
    "results/validation/GSE221156/",
    "Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
)

lt_file <- "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"

outdir <- "results/validation/GSE221156/NicheNet_replication"
reprodir <- "reproducibility"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reprodir, recursive = TRUE, showWarnings = FALSE)

required <- c(de_file, lt_file)

missing <- required[!file.exists(required)]

if (length(missing)) {
    stop(
        "Missing required file(s):\n",
        paste(missing, collapse = "\n")
    )
}

cat("\nINPUT FILES: PASS\n")
for (f in required) {
    cat("  ", f, "\n", sep = "")
}

# =============================================================================
# 2. Load and verify frozen NicheNet prior
# =============================================================================

expected_md5 <- "b09606b04b2d4490418d9028c0e58b9f"

observed_md5 <- unname(tools::md5sum(lt_file))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FROZEN NICHE-NET PRIOR\n")
cat(strrep("=", 78), "\n", sep = "")

cat("File:", lt_file, "\n")
cat("Observed MD5:", observed_md5, "\n")
cat("Expected MD5:", expected_md5, "\n")

if (!identical(observed_md5, expected_md5)) {
    stop(
        "NicheNet prior MD5 does not match the frozen Script-06 prior."
    )
}

ligand_target_matrix <- readRDS(lt_file)

if (
    is.null(rownames(ligand_target_matrix)) ||
    is.null(colnames(ligand_target_matrix))
) {
    stop("Ligand-target matrix lacks row/column names.")
}

cat(
    "Matrix:",
    nrow(ligand_target_matrix),
    "targets x",
    ncol(ligand_target_matrix),
    "ligands\n"
)

if (
    nrow(ligand_target_matrix) != 33354 ||
    ncol(ligand_target_matrix) != 1226
) {
    stop("Unexpected dimensions for frozen NicheNet prior.")
}

cat("Frozen prior verification: PASS\n")

# =============================================================================
# 3. Read independent GSE221156 donor-level edgeR results
# =============================================================================

de <- read.csv(
    de_file,
    check.names = FALSE
)

required_cols <- c(
    "gene",
    "logFC",
    "FDR"
)

if (!all(required_cols %in% names(de))) {
    stop(
        "GSE221156 edgeR table must contain: ",
        paste(required_cols, collapse = ", ")
    )
}

de$gene <- trimws(as.character(de$gene))

de <- de[
    !is.na(de$gene) &
    de$gene != "" &
    !duplicated(de$gene),
    ,
    drop = FALSE
]

cat("\n", strrep("=", 78), "\n", sep = "")
cat("INDEPENDENT GSE221156 RECEIVER DEFINITION\n")
cat(strrep("=", 78), "\n", sep = "")

cat("Genes tested by edgeR:", nrow(de), "\n")

# Prespecified independent receiver:
# FDR < 0.05 AND |log2FC| >= 0.25
receiver_raw <- unique(
    de$gene[
        !is.na(de$FDR) &
        !is.na(de$logFC) &
        de$FDR < 0.05 &
        abs(de$logFC) >= 0.25
    ]
)

receiver_up <- unique(
    de$gene[
        !is.na(de$FDR) &
        !is.na(de$logFC) &
        de$FDR < 0.05 &
        abs(de$logFC) >= 0.25 &
        de$logFC > 0
    ]
)

receiver_down <- unique(
    de$gene[
        !is.na(de$FDR) &
        !is.na(de$logFC) &
        de$FDR < 0.05 &
        abs(de$logFC) >= 0.25 &
        de$logFC < 0
    ]
)

cat("Independent receiver raw:", length(receiver_raw), "\n")
cat("  T2D-up:", length(receiver_up), "\n")
cat("  T2D-down:", length(receiver_down), "\n")

if (length(receiver_raw) != 796) {
    warning(
        "Expected 796 receiver genes from audit, observed ",
        length(receiver_raw),
        "."
    )
}

# =============================================================================
# 4. Define independent NicheNet background
# =============================================================================

# Background = all genes actually tested by edgeR in GSE221156,
# restricted to genes represented in the frozen NicheNet prior.

background_raw <- unique(de$gene)

background_nn <- intersect(
    background_raw,
    rownames(ligand_target_matrix)
)

receiver_nn <- intersect(
    receiver_raw,
    background_nn
)

receiver_absent <- setdiff(
    receiver_raw,
    receiver_nn
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("RECEIVER / BACKGROUND OVERLAP\n")
cat(strrep("=", 78), "\n", sep = "")

cat("GSE221156 edgeR tested genes:", length(background_raw), "\n")
cat(
    "NicheNet-compatible background:",
    length(background_nn),
    "\n"
)

cat("Independent receiver raw:", length(receiver_raw), "\n")
cat(
    "Independent receiver in NicheNet:",
    length(receiver_nn),
    "\n"
)

cat(
    "Independent receiver absent from NicheNet:",
    length(receiver_absent),
    "\n"
)

if (length(receiver_nn) < 20) {
    stop("Too few independent receiver genes represented in NicheNet.")
}

if (!all(receiver_nn %in% background_nn)) {
    stop("Receiver is not fully contained in NicheNet background.")
}

# =============================================================================
# 5. FIXED sensory-source candidate list
# =============================================================================

# IMPORTANT:
# This list is frozen from the discovery analysis.
# It is NOT redefined using GSE221156.

sensory13 <- c(
    "RPH3A",
    "EDIL3",
    "SPP1",
    "SEMA7A",
    "PPIA",
    "FGF13",
    "FGF12",
    "SST",
    "CNTN2",
    "LGALS8",
    "RSPO2",
    "FLRT3",
    "PGF"
)

candidate_nn <- sensory13[
    sensory13 %in% colnames(ligand_target_matrix)
]

candidate_absent <- setdiff(
    sensory13,
    candidate_nn
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FIXED LIGAND OVERLAP\n")
cat(strrep("=", 78), "\n", sep = "")

cat("Sensory-source candidates raw:", length(sensory13), "\n")
cat(
    "Sensory-source candidates in NicheNet:",
    length(candidate_nn),
    "\n"
)

cat(
    "Sensory-source candidates absent from NicheNet:",
    length(candidate_absent),
    "\n"
)

if (length(candidate_absent)) {
    cat(
        "Absent sensory ligands:",
        paste(candidate_absent, collapse = ";"),
        "\n"
    )
}

expected_absent <- c(
    "RPH3A",
    "PPIA",
    "LGALS8"
)

if (!setequal(candidate_absent, expected_absent)) {
    warning(
        "Candidate overlap differs from frozen Script-06 analysis."
    )
}

# =============================================================================
# 6. Independent GSE221156 NicheNet ligand activity
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("RUNNING INDEPENDENT NICHE-NET REPLICATION\n")
cat(strrep("=", 78), "\n", sep = "")

ligand_activity <- nichenetr::predict_ligand_activities(
    geneset = receiver_nn,
    background_expressed_genes = background_nn,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = candidate_nn
)

ligand_activity <- as.data.frame(
    ligand_activity,
    stringsAsFactors = FALSE
)

required_activity_cols <- c(
    "test_ligand",
    "auroc",
    "aupr",
    "aupr_corrected",
    "pearson"
)

if (!all(required_activity_cols %in% names(ligand_activity))) {
    stop("Unexpected NicheNet ligand-activity output.")
}

ligand_activity <- ligand_activity[
    order(
        ligand_activity$aupr_corrected,
        decreasing = TRUE
    ),
    ,
    drop = FALSE
]

ligand_activity$rank_aupr_corrected <- seq_len(
    nrow(ligand_activity)
)

ligand_activity$positive_aupr_corrected <-
    ligand_activity$aupr_corrected > 0

ligand_activity$analysis <-
    "GSE221156_Independent_T2D_Receiver796"

rownames(ligand_activity) <- NULL

cat("\nINDEPENDENT NICHE-NET RANKING\n")
print(
    ligand_activity,
    row.names = FALSE
)

# =============================================================================
# 7. Restore all 13 candidates, including non-evaluable ligands
# =============================================================================

all13 <- data.frame(
    test_ligand = sensory13,
    stringsAsFactors = FALSE
)

all13 <- merge(
    all13,
    ligand_activity,
    by = "test_ligand",
    all.x = TRUE,
    sort = FALSE
)

# Restore original frozen order
all13 <- all13[
    match(sensory13, all13$test_ligand),
    ,
    drop = FALSE
]

all13$in_NicheNet_prior <-
    all13$test_ligand %in% candidate_nn

rownames(all13) <- NULL

# =============================================================================
# 8. Compare with frozen GSE81608 discovery ranking
# =============================================================================

# These values are frozen from Script 06.
# They are included only for rank-replication comparison.
# They are NOT used to define or alter the GSE221156 analysis.

discovery <- data.frame(
    test_ligand = c(
        "PGF",
        "FLRT3",
        "SST",
        "SEMA7A",
        "SPP1",
        "FGF12",
        "FGF13",
        "EDIL3",
        "RSPO2",
        "CNTN2"
    ),
    discovery_rank = 1:10,
    discovery_aupr_corrected = c(
        7.624951e-04,
        6.961810e-04,
        5.565271e-04,
        5.222301e-04,
        4.570261e-04,
        3.698300e-04,
        2.975972e-04,
        2.597106e-04,
        1.323636e-04,
       -1.365005e-05
    ),
    stringsAsFactors = FALSE
)

comparison <- merge(
    discovery,
    ligand_activity[
        ,
        c(
            "test_ligand",
            "rank_aupr_corrected",
            "aupr_corrected",
            "positive_aupr_corrected"
        ),
        drop = FALSE
    ],
    by = "test_ligand",
    all.x = TRUE,
    sort = FALSE
)

names(comparison)[
    names(comparison) == "rank_aupr_corrected"
] <- "validation_rank"

names(comparison)[
    names(comparison) == "aupr_corrected"
] <- "validation_aupr_corrected"

names(comparison)[
    names(comparison) == "positive_aupr_corrected"
] <- "validation_positive_aupr_corrected"

comparison <- comparison[
    order(comparison$discovery_rank),
    ,
    drop = FALSE
]

comparison$rank_change_validation_minus_discovery <-
    comparison$validation_rank -
    comparison$discovery_rank

comparison$discovery_positive_aupr_corrected <-
    comparison$discovery_aupr_corrected > 0

comparison$positive_in_both <-
    comparison$discovery_positive_aupr_corrected &
    comparison$validation_positive_aupr_corrected

# =============================================================================
# 9. Rank replication statistics
# =============================================================================

tested_comparison <- comparison[
    !is.na(comparison$validation_rank),
    ,
    drop = FALSE
]

if (nrow(tested_comparison) >= 3) {

    rank_test <- suppressWarnings(
        cor.test(
            tested_comparison$discovery_rank,
            tested_comparison$validation_rank,
            method = "spearman",
            exact = FALSE
        )
    )

    spearman_rho <- unname(rank_test$estimate)
    spearman_p <- rank_test$p.value

} else {

    spearman_rho <- NA_real_
    spearman_p <- NA_real_
}

discovery_top5 <- discovery$test_ligand[
    discovery$discovery_rank <= 5
]

validation_top5 <- ligand_activity$test_ligand[
    ligand_activity$rank_aupr_corrected <= 5
]

top5_overlap <- intersect(
    discovery_top5,
    validation_top5
)

main_candidates <- c(
    "SEMA7A",
    "EDIL3",
    "FLRT3",
    "PGF",
    "SPP1"
)

main_summary <- comparison[
    comparison$test_ligand %in% main_candidates,
    ,
    drop = FALSE
]

main_summary <- main_summary[
    match(main_candidates, main_summary$test_ligand),
    ,
    drop = FALSE
]

# =============================================================================
# 10. Summary
# =============================================================================

summary_table <- data.frame(
    metric = c(
        "GSE221156_edgeR_tested_genes",
        "Independent_receiver_raw",
        "Independent_receiver_T2D_up",
        "Independent_receiver_T2D_down",
        "NicheNet_background",
        "Independent_receiver_in_NicheNet",
        "Fixed_candidate_ligands",
        "Candidates_evaluable_in_NicheNet",
        "Candidates_positive_corrected_AUPR_validation",
        "Candidates_positive_corrected_AUPR_both_datasets",
        "Spearman_rho_discovery_vs_validation_rank",
        "Spearman_P_discovery_vs_validation_rank",
        "Top5_overlap_n"
    ),
    value = c(
        length(background_raw),
        length(receiver_raw),
        length(receiver_up),
        length(receiver_down),
        length(background_nn),
        length(receiver_nn),
        length(sensory13),
        length(candidate_nn),
        sum(
            ligand_activity$aupr_corrected > 0,
            na.rm = TRUE
        ),
        sum(
            comparison$positive_in_both,
            na.rm = TRUE
        ),
        spearman_rho,
        spearman_p,
        length(top5_overlap)
    ),
    stringsAsFactors = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("DISCOVERY -> INDEPENDENT VALIDATION COMPARISON\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    comparison,
    row.names = FALSE
)

cat("\nRANK REPLICATION\n")
cat(
    "Spearman rho:",
    format(spearman_rho, digits = 5),
    "\n"
)

cat(
    "Spearman P:",
    format(spearman_p, scientific = TRUE, digits = 5),
    "\n"
)

cat(
    "Discovery top 5:",
    paste(discovery_top5, collapse = ", "),
    "\n"
)

cat(
    "Validation top 5:",
    paste(validation_top5, collapse = ", "),
    "\n"
)

cat(
    "Top-5 overlap:",
    length(top5_overlap),
    "/ 5\n"
)

cat(
    "Shared top-5 ligands:",
    ifelse(
        length(top5_overlap),
        paste(top5_overlap, collapse = ", "),
        "none"
    ),
    "\n"
)

cat("\nMAIN MANUSCRIPT CANDIDATES\n")
print(
    main_summary,
    row.names = FALSE
)

cat("\nSUMMARY\n")
print(
    summary_table,
    row.names = FALSE
)

# =============================================================================
# 11. Save outputs
# =============================================================================

activity_file <- file.path(
    outdir,
    "Table_P3_GSE221156_Independent_NicheNet_13Candidates.csv"
)

comparison_file <- file.path(
    outdir,
    "Table_P3_GSE221156_vs_GSE81608_NicheNet_RankReplication.csv"
)

main_file <- file.path(
    outdir,
    "Table_P3_GSE221156_NicheNet_MainCandidates.csv"
)

summary_file <- file.path(
    outdir,
    "Table_P3_GSE221156_NicheNet_ReplicationSummary.csv"
)

receiver_file <- file.path(
    outdir,
    "Table_P3_GSE221156_Independent_Receiver796.csv"
)

write.csv(
    all13,
    activity_file,
    row.names = FALSE
)

write.csv(
    comparison,
    comparison_file,
    row.names = FALSE
)

write.csv(
    main_summary,
    main_file,
    row.names = FALSE
)

write.csv(
    summary_table,
    summary_file,
    row.names = FALSE
)

receiver_table <- de[
    de$gene %in% receiver_raw,
    ,
    drop = FALSE
]

receiver_table <- receiver_table[
    order(
        receiver_table$FDR,
        -abs(receiver_table$logFC)
    ),
    ,
    drop = FALSE
]

write.csv(
    receiver_table,
    receiver_file,
    row.names = FALSE
)

rds_file <- file.path(
    "data/processed",
    "GSE221156_Independent_NicheNet_16a.rds"
)

saveRDS(
    list(
        receiver_raw = receiver_raw,
        receiver_up = receiver_up,
        receiver_down = receiver_down,
        receiver_nichenet = receiver_nn,
        receiver_absent_nichenet = receiver_absent,
        background_nichenet = background_nn,
        fixed_candidates = sensory13,
        evaluable_candidates = candidate_nn,
        ligand_activity = ligand_activity,
        all13 = all13,
        discovery_validation_comparison = comparison,
        main_candidates = main_summary,
        summary = summary_table,
        prior_file = lt_file,
        prior_md5 = observed_md5
    ),
    rds_file
)

session_file <- file.path(
    reprodir,
    "16a_GSE221156_independent_NicheNet_sessionInfo.txt"
)

capture.output(
    sessionInfo(),
    file = session_file
)

# =============================================================================
# 12. Final interpretation guardrail
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("OUTPUT FILES\n")
cat(strrep("=", 78), "\n", sep = "")

for (f in c(
    activity_file,
    comparison_file,
    main_file,
    summary_file,
    receiver_file,
    rds_file,
    session_file
)) {
    cat("  ", f, "\n", sep = "")
}

cat("\nINTERPRETATION GUARDRAIL\n")
cat(
    "  GSE221156 was used to define an independent beta-cell T2D receiver program.\n"
)
cat(
    "  The 13 sensory-source ligands were fixed before this validation analysis.\n"
)
cat(
    "  NicheNet activity indicates regulatory compatibility, not direct signaling.\n"
)
cat(
    "  Corrected AUPR > 0 is not a statistical significance test.\n"
)
cat(
    "  Empirical significance will be evaluated separately in Script 16c.\n"
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16a COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 78), "\n", sep = "")
