#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(edgeR)
    library(nichenetr)
})

cat(strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16c — EXPRESSION-MATCHED NICHE-NET PERMUTATION NULL\n")
cat(strrep("=", 78), "\n", sep = "")

# =============================================================================
# 0. Runtime settings
# =============================================================================

n_perm <- as.integer(Sys.getenv("N_PERM", "10000"))
run_tag <- Sys.getenv("RUN_TAG", "final")
n_bins <- 10L

detected_cores <- parallel::detectCores(logical = FALSE)

if (is.na(detected_cores)) detected_cores <- 4L

default_cores <- max(
    1L,
    min(8L, detected_cores - 1L)
)

n_cores <- as.integer(
    Sys.getenv(
        "N_CORES",
        as.character(default_cores)
    )
)

chunk_size <- as.integer(
    Sys.getenv("CHUNK_SIZE", "250")
)

if (
    is.na(n_perm) || n_perm < 1 ||
    is.na(n_cores) || n_cores < 1 ||
    is.na(chunk_size) || chunk_size < 1
) {
    stop("Invalid runtime settings.")
}

cat("nichenetr version:",
    as.character(packageVersion("nichenetr")), "\n")
cat("Permutations per cohort:", n_perm, "\n")
cat("Expression strata:", n_bins, "\n")
cat("Parallel cores:", n_cores, "\n")
cat("Chunk size:", chunk_size, "\n")
cat("Run tag:", run_tag, "\n")

# =============================================================================
# 1. Frozen inputs
# =============================================================================

dge_file <-
    "results/islet/GSE81608_Beta_edgeR_Primary_DGEList.rds"

receiver81608_file <-
    "results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv"

broad_file <-
    "results/islet/Table_P3_HC_Ligand_1492_vs_1600_Sensitivity.csv"

obs81608_file <-
    "results/nichenet/Table_P3_NicheNet_BroadHC_BetaSupported_LigandActivity.csv"

de221156_file <- paste0(
    "results/validation/GSE221156/",
    "Table_P3_GSE221156_Beta_T2D_vs_Normal_edgeR.csv"
)

obs221156_file <- paste0(
    "results/validation/GSE221156/NicheNet_broad_replication/",
    "Table_P3_GSE221156_BroadHC_NicheNet_LigandActivity.csv"
)

lt_file <-
    "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"

outdir <- file.path(
    "results/validation/GSE221156/NicheNet_permutation",
    run_tag
)

reprodir <- "reproducibility"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reprodir, recursive = TRUE, showWarnings = FALSE)

required <- c(
    dge_file,
    receiver81608_file,
    broad_file,
    obs81608_file,
    de221156_file,
    obs221156_file,
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

# =============================================================================
# 2. Frozen NicheNet prior
# =============================================================================

expected_md5 <-
    "b09606b04b2d4490418d9028c0e58b9f"

observed_md5 <-
    unname(tools::md5sum(lt_file))

if (!identical(observed_md5, expected_md5)) {
    stop("Frozen NicheNet prior MD5 mismatch.")
}

lt <- readRDS(lt_file)

if (
    nrow(lt) != 33354 ||
    ncol(lt) != 1226
) {
    stop("Unexpected NicheNet prior dimensions.")
}

cat("\nFROZEN NICHE-NET PRIOR\n")
cat("MD5:", observed_md5, "\n")
cat(
    "Matrix:",
    nrow(lt), "targets x",
    ncol(lt), "ligands\n"
)

# =============================================================================
# 3. Reconstruct frozen 42-ligand broad universe
# =============================================================================

broad <- read.csv(
    broad_file,
    check.names = FALSE
)

if (
    !"ligand" %in% names(broad) ||
    !"beta_supported_in_both_representations" %in% names(broad)
) {
    stop("Unexpected broad-ligand table structure.")
}

keep_broad <-
    broad$beta_supported_in_both_representations %in%
    c(TRUE, "TRUE", "True", 1, "1")

broad_raw <- unique(
    trimws(
        as.character(
            broad$ligand[keep_broad]
        )
    )
)

broad_raw <- broad_raw[
    !is.na(broad_raw) &
    broad_raw != ""
]

broad_ligands <- intersect(
    broad_raw,
    colnames(lt)
)

cat("\nBROAD LIGAND UNIVERSE\n")
cat("Raw:", length(broad_raw), "\n")
cat("NicheNet-evaluable:", length(broad_ligands), "\n")

if (length(broad_raw) != 45) {
    stop("Expected 45 raw broad ligands.")
}

if (length(broad_ligands) != 42) {
    stop("Expected 42 NicheNet-evaluable broad ligands.")
}

# Fixed alphabetical order for permutation matrices
ligand_order <- sort(broad_ligands)

# =============================================================================
# 4. Helper functions
# =============================================================================

make_expression_bins <- function(
    background_genes,
    expression_values,
    n_bins = 10L
) {

    if (is.null(names(expression_values))) {
        stop("Expression vector must be named.")
    }

    missing_expr <- setdiff(
        background_genes,
        names(expression_values)
    )

    if (length(missing_expr)) {
        stop(
            "Missing expression values for ",
            length(missing_expr),
            " background genes."
        )
    }

    expr <- expression_values[background_genes]

    if (any(!is.finite(expr))) {
        stop("Non-finite expression values in background.")
    }

    # Deterministic rank-based expression strata.
    # Gene name breaks expression ties reproducibly.
    ord <- order(
        expr,
        names(expr)
    )

    n <- length(ord)

    bin_ord <- pmin(
        n_bins,
        ceiling(
            seq_len(n) * n_bins / n
        )
    )

    bins <- integer(n)
    bins[ord] <- bin_ord
    names(bins) <- names(expr)

    bins
}


sample_expression_matched <- function(
    background_genes,
    bins,
    receiver_genes
) {

    receiver_bins <- bins[receiver_genes]

    if (any(is.na(receiver_bins))) {
        stop("Receiver genes missing expression-bin assignment.")
    }

    needed <- table(receiver_bins)

    sampled <- unlist(
        lapply(
            names(needed),
            function(b) {

                pool <- background_genes[
                    bins[background_genes] ==
                    as.integer(b)
                ]

                k <- as.integer(needed[[b]])

                if (length(pool) < k) {
                    stop(
                        "Expression bin ",
                        b,
                        " has insufficient genes."
                    )
                }

                sample(
                    pool,
                    size = k,
                    replace = FALSE
                )
            }
        ),
        use.names = FALSE
    )

    sampled <- unique(sampled)

    if (length(sampled) != length(receiver_genes)) {
        stop(
            "Matched null receiver has wrong size: ",
            length(sampled),
            " vs ",
            length(receiver_genes)
        )
    }

    sampled
}


run_activity <- function(
    geneset,
    background,
    lt_sub,
    ligands
) {

    x <- nichenetr::predict_ligand_activities(
        geneset = geneset,
        background_expressed_genes = background,
        ligand_target_matrix = lt_sub,
        potential_ligands = ligands
    )

    x <- as.data.frame(
        x,
        stringsAsFactors = FALSE
    )

    v <- x$aupr_corrected
    names(v) <- x$test_ligand

    v[ligands]
}


verify_observed <- function(
    saved,
    recomputed,
    cohort,
    tolerance = 1e-10
) {

    saved <- saved[
        match(
            names(recomputed),
            saved$test_ligand
        ),
        ,
        drop = FALSE
    ]

    if (any(is.na(saved$test_ligand))) {
        stop(cohort, ": saved activity table missing ligands.")
    }

    delta <- max(
        abs(
            saved$aupr_corrected -
            as.numeric(recomputed)
        ),
        na.rm = TRUE
    )

    cat(
        cohort,
        "max |saved - recomputed|:",
        format(delta, scientific = TRUE),
        "\n"
    )

    if (!is.finite(delta) || delta > tolerance) {
        stop(
            cohort,
            ": observed activity does not reproduce saved analysis."
        )
    }

    TRUE
}


run_null <- function(
    cohort,
    receiver,
    background,
    bins,
    lt_sub,
    ligands,
    n_perm,
    n_cores,
    chunk_size,
    seed_base,
    checkpoint_file
) {

    cat("\n", strrep("-", 78), "\n", sep = "")
    cat("PERMUTATION NULL:", cohort, "\n")
    cat(strrep("-", 78), "\n", sep = "")

    cat("Receiver size:", length(receiver), "\n")
    cat("Background size:", length(background), "\n")
    cat("Permutations:", n_perm, "\n")

    receiver_bin_counts <- table(
        factor(
            bins[receiver],
            levels = seq_len(n_bins)
        )
    )

    cat("Observed receiver genes per expression bin:\n")
    print(receiver_bin_counts)

    # Fixed independent seeds make the run reproducible regardless
    # of the number of parallel workers.
    set.seed(seed_base)

    perm_seeds <- sample.int(
        .Machine$integer.max,
        size = n_perm,
        replace = FALSE
    )

    completed <- 0L
    scores <- matrix(
        numeric(0),
        nrow = 0,
        ncol = length(ligands),
        dimnames = list(NULL, ligands)
    )

    if (file.exists(checkpoint_file)) {

        old <- readRDS(checkpoint_file)

        if (
            !identical(old$n_perm, n_perm) ||
            !identical(old$ligands, ligands) ||
            !identical(old$seed_base, seed_base) ||
            !identical(old$receiver_size, length(receiver))
        ) {
            stop(
                cohort,
                ": checkpoint settings do not match current run."
            )
        }

        scores <- old$scores
        completed <- nrow(scores)

        cat(
            "Resuming from checkpoint:",
            completed,
            "permutations completed.\n"
        )
    }

    if (completed >= n_perm) {
        cat("Checkpoint already complete.\n")
        return(scores)
    }

    starts <- seq(
        completed + 1L,
        n_perm,
        by = chunk_size
    )

    for (start_i in starts) {

        end_i <- min(
            start_i + chunk_size - 1L,
            n_perm
        )

        idx <- start_i:end_i

        t0 <- proc.time()[["elapsed"]]

        ans <- parallel::mclapply(
            idx,
            function(i) {

                set.seed(perm_seeds[i])

                gs <- sample_expression_matched(
                    background_genes = background,
                    bins = bins,
                    receiver_genes = receiver
                )

                run_activity(
                    geneset = gs,
                    background = background,
                    lt_sub = lt_sub,
                    ligands = ligands
                )
            },
            mc.cores = n_cores,
            mc.preschedule = TRUE,
            mc.set.seed = FALSE
        )

        bad <- vapply(
            ans,
            inherits,
            logical(1),
            what = "try-error"
        )

        if (any(bad)) {
            stop(
                cohort,
                ": permutation worker failure in chunk ",
                start_i, "-", end_i
            )
        }

        chunk <- do.call(
            rbind,
            ans
        )

        chunk <- chunk[
            ,
            ligands,
            drop = FALSE
        ]

        scores <- rbind(
            scores,
            chunk
        )

        saveRDS(
            list(
                cohort = cohort,
                n_perm = n_perm,
                seed_base = seed_base,
                receiver_size = length(receiver),
                ligands = ligands,
                scores = scores
            ),
            checkpoint_file
        )

        elapsed <- proc.time()[["elapsed"]] - t0

        cat(
            sprintf(
                "  %s: %d/%d completed | chunk %.1f sec\n",
                cohort,
                nrow(scores),
                n_perm,
                elapsed
            )
        )
    }

    scores
}


summarize_null <- function(
    cohort,
    observed,
    null_scores
) {

    observed <- observed[
        match(
            colnames(null_scores),
            observed$test_ligand
        ),
        ,
        drop = FALSE
    ]

    obs <- observed$aupr_corrected
    names(obs) <- observed$test_ligand

    null_mean <- colMeans(
        null_scores,
        na.rm = TRUE
    )

    null_sd <- apply(
        null_scores,
        2,
        sd,
        na.rm = TRUE
    )

    null_q95 <- apply(
        null_scores,
        2,
        quantile,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
    )

    null_q99 <- apply(
        null_scores,
        2,
        quantile,
        probs = 0.99,
        na.rm = TRUE,
        names = FALSE
    )

    # Conservative one-sided empirical P:
    # null corrected AUPR >= observed corrected AUPR.
    exceed <- vapply(
        names(obs),
        function(lig) {
            sum(
                null_scores[, lig] >= obs[[lig]],
                na.rm = TRUE
            )
        },
        integer(1)
    )

    empirical_p <-
        (1 + exceed) /
        (1 + nrow(null_scores))

    empirical_fdr <- p.adjust(
        empirical_p,
        method = "BH"
    )

    z <- (
        obs - null_mean[names(obs)]
    ) / null_sd[names(obs)]

    z[
        !is.finite(z)
    ] <- NA_real_

    out <- data.frame(
        cohort = cohort,
        test_ligand = names(obs),
        observed_aupr_corrected =
            as.numeric(obs),
        null_mean =
            as.numeric(null_mean[names(obs)]),
        null_sd =
            as.numeric(null_sd[names(obs)]),
        null_q95 =
            as.numeric(null_q95[names(obs)]),
        null_q99 =
            as.numeric(null_q99[names(obs)]),
        observed_minus_null_mean =
            as.numeric(
                obs -
                null_mean[names(obs)]
            ),
        null_z =
            as.numeric(z),
        n_null_ge_observed =
            as.integer(exceed),
        n_permutations =
            nrow(null_scores),
        empirical_P =
            as.numeric(empirical_p),
        empirical_FDR =
            as.numeric(empirical_fdr),
        stringsAsFactors = FALSE
    )

    out <- out[
        order(
            out$empirical_P,
            -out$observed_aupr_corrected
        ),
        ,
        drop = FALSE
    ]

    rownames(out) <- NULL

    out
}

# =============================================================================
# 5. GSE81608 inputs
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("GSE81608 SETUP\n")
cat(strrep("=", 78), "\n", sep = "")

dge <- readRDS(dge_file)

genes81608 <- rownames(dge$counts)

if (
    is.null(genes81608) ||
    length(genes81608) != length(dge$AveLogCPM)
) {
    stop("Could not align GSE81608 AveLogCPM with genes.")
}

expr81608 <- as.numeric(dge$AveLogCPM)
names(expr81608) <- genes81608

background81608 <- intersect(
    genes81608,
    rownames(lt)
)

receiver81608_tbl <- read.csv(
    receiver81608_file,
    check.names = FALSE
)

if (!"gene" %in% names(receiver81608_tbl)) {
    stop("GSE81608 receiver table lacks gene column.")
}

receiver81608_raw <- unique(
    receiver81608_tbl$gene
)

receiver81608 <- intersect(
    receiver81608_raw,
    background81608
)

cat("edgeR-tested genes:", length(genes81608), "\n")
cat("NicheNet background:", length(background81608), "\n")
cat("Receiver raw:", length(receiver81608_raw), "\n")
cat("Receiver in NicheNet:", length(receiver81608), "\n")

if (length(background81608) != 10742) {
    warning(
        "Expected GSE81608 NicheNet background = 10742; observed ",
        length(background81608)
    )
}

if (length(receiver81608) != 22) {
    stop(
        "Expected GSE81608 NicheNet receiver = 22; observed ",
        length(receiver81608)
    )
}

bins81608 <- make_expression_bins(
    background_genes = background81608,
    expression_values = expr81608,
    n_bins = n_bins
)

lt81608 <- lt[
    background81608,
    ligand_order,
    drop = FALSE
]

obs81608 <- read.csv(
    obs81608_file,
    check.names = FALSE
)

obs81608_recomputed <- run_activity(
    geneset = receiver81608,
    background = background81608,
    lt_sub = lt81608,
    ligands = ligand_order
)

verify_observed(
    saved = obs81608,
    recomputed = obs81608_recomputed,
    cohort = "GSE81608"
)

# =============================================================================
# 6. GSE221156 inputs
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("GSE221156 SETUP\n")
cat(strrep("=", 78), "\n", sep = "")

de221156 <- read.csv(
    de221156_file,
    check.names = FALSE
)

required_de <- c(
    "gene",
    "logFC",
    "logCPM",
    "FDR"
)

if (!all(required_de %in% names(de221156))) {
    stop("Unexpected GSE221156 edgeR columns.")
}

de221156$gene <- trimws(
    as.character(de221156$gene)
)

de221156 <- de221156[
    !is.na(de221156$gene) &
    de221156$gene != "" &
    !duplicated(de221156$gene),
    ,
    drop = FALSE
]

genes221156 <- de221156$gene

expr221156 <- de221156$logCPM
names(expr221156) <- genes221156

background221156 <- intersect(
    genes221156,
    rownames(lt)
)

receiver221156_raw <- unique(
    de221156$gene[
        !is.na(de221156$FDR) &
        !is.na(de221156$logFC) &
        de221156$FDR < 0.05 &
        abs(de221156$logFC) >= 0.25
    ]
)

receiver221156 <- intersect(
    receiver221156_raw,
    background221156
)

cat("edgeR-tested genes:", length(genes221156), "\n")
cat("NicheNet background:", length(background221156), "\n")
cat("Receiver raw:", length(receiver221156_raw), "\n")
cat("Receiver in NicheNet:", length(receiver221156), "\n")

if (length(background221156) != 13712) {
    warning(
        "Expected GSE221156 NicheNet background = 13712; observed ",
        length(background221156)
    )
}

if (length(receiver221156) != 697) {
    stop(
        "Expected GSE221156 NicheNet receiver = 697; observed ",
        length(receiver221156)
    )
}

bins221156 <- make_expression_bins(
    background_genes = background221156,
    expression_values = expr221156,
    n_bins = n_bins
)

lt221156 <- lt[
    background221156,
    ligand_order,
    drop = FALSE
]

obs221156 <- read.csv(
    obs221156_file,
    check.names = FALSE
)

obs221156_recomputed <- run_activity(
    geneset = receiver221156,
    background = background221156,
    lt_sub = lt221156,
    ligands = ligand_order
)

verify_observed(
    saved = obs221156,
    recomputed = obs221156_recomputed,
    cohort = "GSE221156"
)

# =============================================================================
# 7. Expression-distribution audit
# =============================================================================

audit_expression <- function(
    cohort,
    background,
    receiver,
    expr,
    bins
) {

    tab <- data.frame(
        cohort = cohort,
        bin = seq_len(n_bins),
        background_n = as.integer(
            table(
                factor(
                    bins[background],
                    levels = seq_len(n_bins)
                )
            )
        ),
        receiver_n = as.integer(
            table(
                factor(
                    bins[receiver],
                    levels = seq_len(n_bins)
                )
            )
        ),
        stringsAsFactors = FALSE
    )

    tab
}

expr_audit <- rbind(
    audit_expression(
        "GSE81608",
        background81608,
        receiver81608,
        expr81608,
        bins81608
    ),
    audit_expression(
        "GSE221156",
        background221156,
        receiver221156,
        expr221156,
        bins221156
    )
)

cat("\nEXPRESSION-MATCHING AUDIT\n")
print(expr_audit, row.names = FALSE)

write.csv(
    expr_audit,
    file.path(
        outdir,
        "Table_P3_NicheNet_ExpressionMatching_Audit.csv"
    ),
    row.names = FALSE
)

# =============================================================================
# 8. Run / resume permutation nulls
# =============================================================================

checkpoint81608 <- file.path(
    outdir,
    paste0(
        "checkpoint_GSE81608_",
        n_perm,
        "perm.rds"
    )
)

checkpoint221156 <- file.path(
    outdir,
    paste0(
        "checkpoint_GSE221156_",
        n_perm,
        "perm.rds"
    )
)

null81608 <- run_null(
    cohort = "GSE81608",
    receiver = receiver81608,
    background = background81608,
    bins = bins81608,
    lt_sub = lt81608,
    ligands = ligand_order,
    n_perm = n_perm,
    n_cores = n_cores,
    chunk_size = chunk_size,
    seed_base = 8160801L,
    checkpoint_file = checkpoint81608
)

null221156 <- run_null(
    cohort = "GSE221156",
    receiver = receiver221156,
    background = background221156,
    bins = bins221156,
    lt_sub = lt221156,
    ligands = ligand_order,
    n_perm = n_perm,
    n_cores = n_cores,
    chunk_size = chunk_size,
    seed_base = 22115601L,
    checkpoint_file = checkpoint221156
)

# =============================================================================
# 9. Empirical P values
# =============================================================================

res81608 <- summarize_null(
    cohort = "GSE81608",
    observed = obs81608,
    null_scores = null81608
)

res221156 <- summarize_null(
    cohort = "GSE221156",
    observed = obs221156,
    null_scores = null221156
)

write.csv(
    res81608,
    file.path(
        outdir,
        "Table_P3_GSE81608_NicheNet_ExpressionMatched_EmpiricalNull.csv"
    ),
    row.names = FALSE
)

write.csv(
    res221156,
    file.path(
        outdir,
        "Table_P3_GSE221156_NicheNet_ExpressionMatched_EmpiricalNull.csv"
    ),
    row.names = FALSE
)

# =============================================================================
# 10. Cross-cohort replication summary
# =============================================================================

cross <- merge(
    res81608[
        ,
        c(
            "test_ligand",
            "observed_aupr_corrected",
            "empirical_P",
            "empirical_FDR",
            "null_z"
        )
    ],
    res221156[
        ,
        c(
            "test_ligand",
            "observed_aupr_corrected",
            "empirical_P",
            "empirical_FDR",
            "null_z"
        )
    ],
    by = "test_ligand",
    suffixes = c(
        "_GSE81608",
        "_GSE221156"
    ),
    sort = FALSE
)

# Intersection-union/conjunction summary:
# evidence must survive in BOTH cohorts.
cross$conjunction_P <- pmax(
    cross$empirical_P_GSE81608,
    cross$empirical_P_GSE221156
)

cross$conjunction_FDR <- p.adjust(
    cross$conjunction_P,
    method = "BH"
)

cross$empirical_FDR05_both <-
    cross$empirical_FDR_GSE81608 < 0.05 &
    cross$empirical_FDR_GSE221156 < 0.05

cross <- cross[
    order(
        cross$conjunction_P,
        -cross$observed_aupr_corrected_GSE221156
    ),
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

main_res <- cross[
    cross$test_ligand %in% main_candidates,
    ,
    drop = FALSE
]

main_res <- main_res[
    match(
        main_candidates,
        main_res$test_ligand
    ),
    ,
    drop = FALSE
]

write.csv(
    cross,
    file.path(
        outdir,
        "Table_P3_NicheNet_EmpiricalNull_CrossCohort.csv"
    ),
    row.names = FALSE
)

write.csv(
    main_res,
    file.path(
        outdir,
        "Table_P3_NicheNet_EmpiricalNull_MainCandidates.csv"
    ),
    row.names = FALSE
)

# =============================================================================
# 11. Print final results
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("GSE81608 EMPIRICAL NULL — TOP RESULTS\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    head(res81608, 15),
    row.names = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("GSE221156 EMPIRICAL NULL — TOP RESULTS\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    head(res221156, 15),
    row.names = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FIVE MAIN MANUSCRIPT CANDIDATES\n")
cat(strrep("=", 78), "\n", sep = "")

print(
    main_res,
    row.names = FALSE
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("EMPIRICAL SIGNIFICANCE SUMMARY\n")
cat(strrep("=", 78), "\n", sep = "")

cat(
    "GSE81608 ligands empirical FDR < 0.05:",
    sum(res81608$empirical_FDR < 0.05),
    "/ 42\n"
)

cat(
    "GSE221156 ligands empirical FDR < 0.05:",
    sum(res221156$empirical_FDR < 0.05),
    "/ 42\n"
)

cat(
    "Ligands empirical FDR < 0.05 in BOTH cohorts:",
    sum(cross$empirical_FDR05_both),
    "/ 42\n"
)

if (any(cross$empirical_FDR05_both)) {
    cat(
        "Significant in both:",
        paste(
            cross$test_ligand[
                cross$empirical_FDR05_both
            ],
            collapse = ", "
        ),
        "\n"
    )
}

# =============================================================================
# 12. Reproducibility record
# =============================================================================

analysis_record <- data.frame(
    nichenetr_version =
        as.character(packageVersion("nichenetr")),
    prior_md5 = observed_md5,
    broad_ligands_raw = length(broad_raw),
    broad_ligands_evaluable = length(broad_ligands),
    expression_bins = n_bins,
    permutations_per_cohort = n_perm,
    parallel_cores = n_cores,
    GSE81608_background = length(background81608),
    GSE81608_receiver = length(receiver81608),
    GSE221156_background = length(background221156),
    GSE221156_receiver = length(receiver221156),
    empirical_P_formula =
        "(1 + number(null >= observed)) / (1 + N)",
    multiple_testing =
        "Benjamini-Hochberg across 42 ligands within each cohort",
    null_model =
        paste(
            "same receiver size and same expression-decile",
            "composition as observed receiver"
        ),
    stringsAsFactors = FALSE
)

write.csv(
    analysis_record,
    file.path(
        outdir,
        "Table_P3_NicheNet_EmpiricalNull_AnalysisRecord.csv"
    ),
    row.names = FALSE
)

saveRDS(
    list(
        analysis_record = analysis_record,
        expression_audit = expr_audit,
        GSE81608_results = res81608,
        GSE221156_results = res221156,
        cross_cohort = cross,
        main_candidates = main_res
    ),
    file.path(
        outdir,
        "NicheNet_ExpressionMatched_EmpiricalNull_16c.rds"
    )
)

session_file <- file.path(
    reprodir,
    paste0(
        "16c_NicheNet_expression_matched_permutation_",
        run_tag,
        "_sessionInfo.txt"
    )
)

capture.output(
    sessionInfo(),
    file = session_file
)

cat("\nINTERPRETATION GUARDRAIL\n")
cat(
    "  Empirical P values test whether observed corrected AUPR exceeds\n"
)
cat(
    "  expression-matched random receiver sets of identical size.\n"
)
cat(
    "  These tests assess receiver-program specificity within the frozen\n"
)
cat(
    "  NicheNet prior; they do not demonstrate direct inter-organ signaling.\n"
)
cat(
    "  Candidate ligands and the 42-ligand universe were fixed before\n"
)
cat(
    "  permutation testing.\n"
)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SCRIPT 16c COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 78), "\n", sep = "")
