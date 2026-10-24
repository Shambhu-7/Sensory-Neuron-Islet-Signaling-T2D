#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(nichenetr)
})

cat(strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16b — BROAD LIGAND-UNIVERSE NICHE-NET REPLICATION\n")
cat(strrep("=", 78), "\n", sep = "")

cat("nichenetr version:",
    as.character(packageVersion("nichenetr")), "\n")

# =============================================================================
# 1. Frozen inputs
# =============================================================================

broad_file <-
    "results/islet/Table_P3_HC_Ligand_1492_vs_1600_Sensitivity.csv"

discovery_broad_file <-
    "results/nichenet/Table_P3_NicheNet_BroadHC_BetaSupported_LigandActivity.csv"

validation_de_file <- paste0(
    "results/validation/GSE221156/",
    "Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
)

lt_file <-
    "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"

outdir <-
    "results/validation/GSE221156/NicheNet_broad_replication"

reprodir <- "reproducibility"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reprodir, recursive = TRUE, showWarnings = FALSE)

required <- c(
    broad_file,
    discovery_broad_file,
    validation_de_file,
    lt_file
)

missing <- required[!file.exists(required)]

if (length(missing)) {
    stop(
        "Missing required file(s):\n",
        paste(missing, collapse = "\n")
    )
}

cat("\nINPUT FILES: PASS\n")
for (f in required) cat("  ", f, "\n", sep = "")

# =============================================================================
# 2. Verify frozen NicheNet prior
# =============================================================================

expected_md5 <-
    "b09606b04b2d4490418d9028c0e58b9f"

observed_md5 <-
    unname(tools::md5sum(lt_file))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FROZEN NICHE-NET PRIOR\n")
cat(strrep("=", 78), "\n", sep = "")

cat("Observed MD5:", observed_md5, "\n")
cat("Expected MD5:", expected_md5, "\n")

if (!identical(observed_md5, expected_md5)) {
    stop("Frozen NicheNet prior MD5 mismatch.")
}

ligand_target_matrix <- readRDS(lt_file)

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
    stop("Unexpected NicheNet prior dimensions.")
}

# =============================================================================
# 3. Reconstruct the exact broad ligand universe from Script 06
# =============================================================================

broad <- read.csv(
    broad_file,
    check.names = FALSE
)

if (!"ligand" %in% names(broad)) {
    stop("Broad ligand table lacks ligand column.")
}

if (
    "beta_supported_in_both_representations" %in%
    names(broad)
) {

    broad_keep <-
        broad$beta_supported_in_both_representations %in%
        c(TRUE, "TRUE", "True", 1, "1")

    broad_ligands_raw <-
        unique(
            broad$ligand[broad_keep]
        )

} else {

    stop(
        paste(
            "Expected beta_supported_in_both_representations",
            "column is absent; refusing to redefine broad universe."
        )
    )
}

broad_ligands_raw <-
    trimws(as.character(broad_ligands_raw))

broad_ligands_raw <-
    broad_ligands_raw[
        !is.na(broad_ligands_raw) &
        broad_ligands_raw != ""
    ]

broad_ligands <-
    intersect(
        broad_ligands_raw,
        colnames(ligand_target_matrix)
    )

broad_absent <-
    setdiff(
        broad_ligands_raw,
        broad_ligands
    )

cat("\n", strrep("=", 78), "\n", sep = "")
cat("BROAD LIGAND UNIVERSE\n")
cat(strrep("=", 78), "\n", sep = "")

cat(
    "Broad HC beta-supported ligands raw:",
    length(broad_ligands_raw),
    "\n"
)

cat(
    "Broad HC beta-supported ligands in NicheNet:",
    length(broad_ligands),
    "\n"
)

cat(
    "Broad ligands absent from NicheNet:",
    length(broad_absent),
    "\n"
)

if (length(broad_absent)) {
    cat(
        "Absent:",
        paste(broad_absent, collapse = ";"),
        "\n"
    )
}

# Frozen Script-06 audit expected 45 raw / 42 evaluable.
if (length(broad_ligands_raw) != 45) {
    warning(
        "Expected 45 broad ligands from Script 06; observed ",
        length(broad_ligands_raw)
    )
}

if (length(broad_ligands) != 42) {
    warning(
        "Expected 42 NicheNet-evaluable broad ligands; observed ",
        length(broad_ligands)
    )
}

# =============================================================================
# 4. Read FROZEN GSE81608 broad-universe activity
# =============================================================================

discovery <- read.csv(
    discovery_broad_file,
    check.names = FALSE
)

required_activity_cols <- c(
    "test_ligand",
    "auroc",
    "aupr",
    "aupr_corrected",
    "pearson"
)

if (!all(required_activity_cols %in% names(discovery))) {
    stop("Frozen GSE81608 broad activity table has unexpected columns.")
}

# Re-establish ranking exactly as Script 06:
# corrected AUPR, then AUPR, then Pearson.
discovery <- discovery[
    order(
        -discovery$aupr_corrected,
        -discovery$aupr,
        -discovery$pearson
    ),
    ,
    drop = FALSE
]

discovery$discovery_broad_rank <- seq_len(nrow(discovery))

rownames(discovery) <- NULL

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FROZEN GSE81608 BROAD RANKING\n")
cat(strrep("=", 78), "\n", sep = "")

cat(
    "Ligands in frozen discovery broad ranking:",
    nrow(discovery),
    "\n"
)

if (!setequal(discovery$test_ligand, broad_ligands)) {
    stop(
        paste(
            "Frozen GSE81608 broad activity ligands do not match",
            "the reconstructed 42-ligand universe."
        )
    )
}

print(
    head(
        discovery[
            ,
            c(
                "test_ligand",
                "aupr_corrected",
                "discovery_broad_rank"
            )
        ],
        15
    ),
    row.names = FALSE
)

# =============================================================================
# 5. Define independent GSE221156 receiver exactly as Script 16a
# =============================================================================

de <- read.csv(
    validation_de_file,
    check.names = FALSE
)

required_de_cols <- c(
    "gene",
    "logFC",
    "FDR"
)

if (!all(required_de_cols %in% names(de))) {
    stop("GSE221156 edgeR table has unexpected columns.")
}

de$gene <- trimws(as.character(de$gene))

de <- de[
    !is.na(de$gene) &
    de$gene != "" &
    !duplicated(de$gene),
    ,
    drop = FALSE
]

receiver_raw <-
    unique(
        de$gene[
            !is.na(de$FDR) &
            !is.na(de$logFC) &
            de$FDR < 0.05 &
            abs(de$logFC) >= 0.25
        ]
    )

background_raw <- unique(de$gene)

background_nn <-
    intersect(
        background_raw,
        rownames(ligand_target_matrix)
    )

receiver_nn <-
    intersect(
        receiver_raw,
        background_nn
    )

cat("\n", strrep("=", 78), "\n", sep = "")
cat("INDEPENDENT GSE221156 RECEIVER\n")
cat(strrep("=", 78), "\n", sep = "")

cat(
    "edgeR-tested genes:",
    length(background_raw),
    "\n"
)

cat(
    "NicheNet-compatible background:",
    length(background_nn),
    "\n"
)

cat(
    "Independent receiver raw:",
    length(receiver_raw),
    "\n"
)

cat(
    "Independent receiver in NicheNet:",
    length(receiver_nn),
    "\n"
)

if (length(receiver_raw) != 796) {
    warning(
        "Expected 796 raw receiver genes; observed ",
        length(receiver_raw)
    )
}

if (length(receiver_nn) != 697) {
    warning(
        "Expected 697 NicheNet-compatible receiver genes; observed ",
        length(receiver_nn)
    )
}

# =============================================================================
# 6. Run broad-universe NicheNet in GSE221156
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("GSE221156 BROAD-UNIVERSE NICHE-NET ANALYSIS\n")
cat(strrep("=", 78), "\n", sep = "")

validation <-
    nichenetr::predict_ligand_activities(
        geneset = receiver_nn,
        background_expressed_genes = background_nn,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = broad_ligands
    )

validation <- as.data.frame(
    validation,
    stringsAsFactors = FALSE
)

if (!all(required_activity_cols %in% names(validation))) {
    stop("Unexpected GSE221156 NicheNet output.")
}

validation <- validation[
    order(
        -validation$aupr_corrected,
        -validation$aupr,
        -validation$pearson
    ),
    ,
    drop = FALSE
]

validation$validation_broad_rank <-
    seq_len(nrow(validation))

validation$validation_positive <-
    validation$aupr_corrected > 0

rownames(validation) <- NULL

cat("\nTOP 20 GSE221156 BROAD LIGANDS\n")

print(
    head(
        validation[
            ,
            c(
                "test_ligand",
                "aupr_corrected",
                "validation_broad_rank"
            )
        ],
        20
    ),
    row.names = FALSE
)

# =============================================================================
# 7. Compare broad rankings across cohorts
# =============================================================================

comparison <- merge(
    discovery[
        ,
        c(
            "test_ligand",
            "aupr_corrected",
            "discovery_broad_rank"
        )
    ],
    validation[
        ,
        c(
            "test_ligand",
            "aupr_corrected",
            "validation_broad_rank"
        )
    ],
    by = "test_ligand",
    suffixes = c("_GSE81608", "_GSE221156"),
    all = TRUE,
    sort = FALSE
)

comparison$discovery_positive <-
    comparison$aupr_corrected_GSE81608 > 0

comparison$validation_positive <-
    comparison$aupr_corrected_GSE221156 > 0

comparison$positive_in_both <-
    comparison$discovery_positive &
    comparison$validation_positive

comparison$rank_change <-
    comparison$validation_broad_rank -
    comparison$discovery_broad_rank

comparison$discovery_percentile <-
    1 -
    (
        comparison$discovery_broad_rank - 1
    ) /
    nrow(discovery)

comparison$validation_percentile <-
    1 -
    (
        comparison$validation_broad_rank - 1
    ) /
    nrow(validation)

comparison <- comparison[
    order(comparison$discovery_broad_rank),
    ,
    drop = FALSE
]

# =============================================================================
# 8. Fixed sensory-source candidates
# =============================================================================

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

sensory_broad <- comparison[
    comparison$test_ligand %in% sensory13,
    ,
    drop = FALSE
]

sensory_broad <- sensory_broad[
    order(sensory_broad$validation_broad_rank),
    ,
    drop = FALSE
]

main_candidates <- c(
    "SEMA7A",
    "EDIL3",
    "FLRT3",
    "PGF",
    "SPP1"
)

main_broad <- comparison[
    comparison$test_ligand %in% main_candidates,
    ,
    drop = FALSE
]

main_broad <- main_broad[
    match(
        main_candidates,
        main_broad$test_ligand
    ),
    ,
    drop = FALSE
]

# =============================================================================
# 9. Cross-cohort broad-rank concordance
# =============================================================================

rank_test <- suppressWarnings(
    cor.test(
        comparison$discovery_broad_rank,
        comparison$validation_broad_rank,
        method = "spearman",
        exact = FALSE
    )
)

rho <- unname(rank_test$estimate)
rho_p <- rank_test$p.value

top10_discovery <-
    discovery$test_ligand[
        discovery$discovery_broad_rank <= 10
    ]

top10_validation <-
    validation$test_ligand[
        validation$validation_broad_rank <= 10
    ]

top10_overlap <-
    intersect(
        top10_discovery,
        top10_validation
    )

# =============================================================================
# 10. Candidate enrichment near the top of broad universe
# =============================================================================

# Descriptive rank positions only.
# This is NOT the empirical null test; permutations are reserved for Script 16c.

main_broad$discovery_top10 <-
    main_broad$discovery_broad_rank <= 10

main_broad$validation_top10 <-
    main_broad$validation_broad_rank <= 10

main_broad$top10_both <-
    main_broad$discovery_top10 &
    main_broad$validation_top10

main_broad$discovery_top_quartile <-
    main_broad$discovery_broad_rank <=
    ceiling(nrow(discovery) / 4)

main_broad$validation_top_quartile <-
    main_broad$validation_broad_rank <=
    ceiling(nrow(validation) / 4)

main_broad$top_quartile_both <-
    main_broad$discovery_top_quartile &
    main_broad$validation_top_quartile

# =============================================================================
# 11. Print principal results
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("CROSS-COHORT BROAD-RANK CONCORDANCE\n")
cat(strrep("=", 78), "\n", sep = "")

cat(
    "Broad-universe Spearman rho:",
    format(rho, digits = 5),
    "\n"
)

cat(
    "Spearman P:",
    format(
        rho_p,
        scientific = TRUE,
        digits = 5
    ),
    "\n"
)

cat(
    "Discovery top-10:",
    paste(top10_discovery, collapse = ", "),
    "\n"
)

cat(
    "Validation top-10:",
    paste(top10_validation, collapse = ", "),
    "\n"
)

cat(
    "Top-10 overlap:",
    length(top10_overlap),
    "/ 10\n"
)

cat(
    "Shared top-10 ligands:",
    ifelse(
        length(top10_overlap),
        paste(top10_overlap, collapse = ", "),
        "none"
    ),
    "\n"
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FIXED SENSORY-SOURCE CANDIDATES IN BROAD UNIVERSE\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    sensory_broad[
        ,
        c(
            "test_ligand",
            "discovery_broad_rank",
            "validation_broad_rank",
            "aupr_corrected_GSE81608",
            "aupr_corrected_GSE221156",
            "positive_in_both"
        )
    ],
    row.names = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FIVE MAIN MANUSCRIPT CANDIDATES\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    main_broad,
    row.names = FALSE
)

# =============================================================================
# 12. Summary table
# =============================================================================

summary_table <- data.frame(
    metric = c(
        "Broad_HC_beta_supported_raw",
        "Broad_NicheNet_evaluable",
        "GSE221156_receiver_raw",
        "GSE221156_receiver_in_NicheNet",
        "Broad_rank_Spearman_rho",
        "Broad_rank_Spearman_P",
        "Top10_overlap_n",
        "Main_candidates_top10_discovery",
        "Main_candidates_top10_validation",
        "Main_candidates_top10_both",
        "Main_candidates_top_quartile_both"
    ),
    value = c(
        length(broad_ligands_raw),
        length(broad_ligands),
        length(receiver_raw),
        length(receiver_nn),
        rho,
        rho_p,
        length(top10_overlap),
        sum(main_broad$discovery_top10, na.rm = TRUE),
        sum(main_broad$validation_top10, na.rm = TRUE),
        sum(main_broad$top10_both, na.rm = TRUE),
        sum(main_broad$top_quartile_both, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 78), "\n", sep = "")

print(summary_table, row.names = FALSE)

# =============================================================================
# 13. Save outputs
# =============================================================================

validation_file <- file.path(
    outdir,
    "Table_P3_GSE221156_BroadHC_NicheNet_LigandActivity.csv"
)

comparison_file <- file.path(
    outdir,
    "Table_P3_BroadHC_NicheNet_GSE81608_vs_GSE221156.csv"
)

sensory_file <- file.path(
    outdir,
    "Table_P3_SensoryCandidates_BroadUniverse_Ranks.csv"
)

main_file <- file.path(
    outdir,
    "Table_P3_MainCandidates_BroadUniverse_Ranks.csv"
)

summary_file <- file.path(
    outdir,
    "Table_P3_BroadUniverse_ReplicationSummary.csv"
)

write.csv(
    validation,
    validation_file,
    row.names = FALSE
)

write.csv(
    comparison,
    comparison_file,
    row.names = FALSE
)

write.csv(
    sensory_broad,
    sensory_file,
    row.names = FALSE
)

write.csv(
    main_broad,
    main_file,
    row.names = FALSE
)

write.csv(
    summary_table,
    summary_file,
    row.names = FALSE
)

rds_file <- file.path(
    "data/processed",
    "BroadLigandUniverse_NicheNet_16b.rds"
)

saveRDS(
    list(
        broad_ligands_raw = broad_ligands_raw,
        broad_ligands = broad_ligands,
        discovery = discovery,
        validation = validation,
        comparison = comparison,
        sensory_candidates = sensory_broad,
        main_candidates = main_broad,
        summary = summary_table,
        prior_md5 = observed_md5
    ),
    rds_file
)

session_file <- file.path(
    reprodir,
    "16b_broad_ligand_universe_sessionInfo.txt"
)

capture.output(
    sessionInfo(),
    file = session_file
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("OUTPUT FILES\n")
cat(strrep("=", 78), "\n", sep = "")

for (f in c(
    validation_file,
    comparison_file,
    sensory_file,
    main_file,
    summary_file,
    rds_file,
    session_file
)) {
    cat("  ", f, "\n", sep = "")
}

cat("\nINTERPRETATION GUARDRAIL\n")
cat(
    "  The broad universe was defined independently of NicheNet activity.\n"
)
cat(
    "  It contains high-confidence beta-compatible ligands reproduced across\n"
)
cat(
    "  the two GSE81608 representations and represented in the frozen prior.\n"
)
cat(
    "  Broad-universe rank is descriptive and is not a significance test.\n"
)
cat(
    "  Empirical significance will be evaluated separately in Script 16c.\n"
)
cat(
    "  NicheNet activity does not establish direct neuron-to-beta-cell signaling.\n"
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16b COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 78), "\n", sep = "")
