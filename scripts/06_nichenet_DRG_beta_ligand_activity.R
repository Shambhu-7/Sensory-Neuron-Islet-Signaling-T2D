#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(timeout = 3600)

suppressPackageStartupMessages({
    library(nichenetr)
})

cat("\n============================================\n")
cat("SCRIPT 06: NICHE-NET DRG -> BETA-CELL ACTIVITY\n")
cat("============================================\n")

cat(
    "nichenetr version:",
    as.character(packageVersion("nichenetr")),
    "\n"
)

dir.create(
    "data/raw/NicheNet",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "results/nichenet",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "figures",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "reproducibility",
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

receiver_file <-
    "results/islet/Table_P3_Beta_T2D_NicheNet_PrimaryReceiverGenes.csv"

core_file <-
    "results/islet/Table_P3_Beta_T2D_NicheNet_HighConfidenceCore.csv"

dge_file <-
    "results/islet/GSE81608_Beta_edgeR_Primary_DGEList.rds"

sensory_file <-
    "results/islet/Table_P3_SensorySource_BetaCell_HC_Shortlist.csv"

broad_file <-
    "results/islet/Table_P3_HC_Ligand_1492_vs_1600_Sensitivity.csv"

required_files <- c(
    receiver_file,
    core_file,
    dge_file,
    sensory_file,
    broad_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {
    stop(
        "Missing required files: ",
        paste(missing_files, collapse = "; ")
    )
}

receiver <- read.csv(
    receiver_file,
    check.names = FALSE
)

core <- read.csv(
    core_file,
    check.names = FALSE
)

sensory <- read.csv(
    sensory_file,
    check.names = FALSE
)

broad <- read.csv(
    broad_file,
    check.names = FALSE
)

dge <- readRDS(dge_file)

stopifnot(
    "gene" %in% colnames(receiver),
    "gene" %in% colnames(core),
    "ligand" %in% colnames(sensory),
    "ligand" %in% colnames(broad)
)

cat("\nPrimary receiver genes:",
    nrow(receiver), "\n")

cat("High-confidence core genes:",
    nrow(core), "\n")

cat("Sensory-source ligand candidates:",
    nrow(sensory), "\n")

# ------------------------------------------------------------
# 2. Freeze official NicheNet v2 human ligand-target matrix
# ------------------------------------------------------------

lt_file <-
    "data/raw/NicheNet/ligand_target_matrix_nsga2r_final.rds"

lt_url <-
    paste0(
        "https://zenodo.org/record/7074291/files/",
        "ligand_target_matrix_nsga2r_final.rds"
    )

if (!file.exists(lt_file)) {

    cat("\nDownloading official NicheNet ligand-target matrix...\n")

    download.file(
        lt_url,
        destfile = lt_file,
        mode = "wb",
        quiet = FALSE
    )
}

if (!file.exists(lt_file)) {
    stop("NicheNet ligand-target matrix download failed.")
}

lt_md5 <- unname(
    tools::md5sum(lt_file)
)

ligand_target_matrix <-
    readRDS(lt_file)

if (
    is.null(rownames(ligand_target_matrix)) ||
    is.null(colnames(ligand_target_matrix))
) {
    stop(
        "Ligand-target matrix lacks row/column names."
    )
}

cat("\n============================================\n")
cat("NICHE-NET PRIOR MODEL\n")
cat("============================================\n")

cat(
    "Matrix:",
    nrow(ligand_target_matrix),
    "targets x",
    ncol(ligand_target_matrix),
    "ligands\n"
)

cat(
    "MD5:",
    lt_md5,
    "\n"
)

# ------------------------------------------------------------
# 3. Define beta-cell background
#
# The background is the actual 12,200-gene universe tested
# in the primary beta-cell edgeR analysis, intersected with
# genes represented in the NicheNet target model.
# ------------------------------------------------------------

tested_beta_genes <-
    rownames(dge$counts)

background_genes <-
    intersect(
        tested_beta_genes,
        rownames(ligand_target_matrix)
    )

primary_receiver_raw <-
    unique(receiver$gene)

core_receiver_raw <-
    unique(core$gene)

primary_receiver <-
    intersect(
        primary_receiver_raw,
        background_genes
    )

core_receiver <-
    intersect(
        core_receiver_raw,
        background_genes
    )

primary_receiver_missing <-
    setdiff(
        primary_receiver_raw,
        primary_receiver
    )

core_receiver_missing <-
    setdiff(
        core_receiver_raw,
        core_receiver
    )

cat("\n============================================\n")
cat("RECEIVER / BACKGROUND OVERLAP\n")
cat("============================================\n")

cat(
    "edgeR tested genes:",
    length(tested_beta_genes),
    "\n"
)

cat(
    "NicheNet-compatible background:",
    length(background_genes),
    "\n"
)

cat(
    "Primary receiver raw:",
    length(primary_receiver_raw),
    "\n"
)

cat(
    "Primary receiver in NicheNet:",
    length(primary_receiver),
    "\n"
)

cat(
    "Primary receiver absent from NicheNet:",
    length(primary_receiver_missing),
    "\n"
)

if (length(primary_receiver_missing) > 0) {
    cat(
        "Absent primary genes:",
        paste(
            primary_receiver_missing,
            collapse = ";"
        ),
        "\n"
    )
}

cat(
    "Core receiver raw:",
    length(core_receiver_raw),
    "\n"
)

cat(
    "Core receiver in NicheNet:",
    length(core_receiver),
    "\n"
)

if (length(primary_receiver) < 5) {
    stop(
        "Too few primary receiver genes overlap ",
        "the NicheNet target matrix."
    )
}

# ------------------------------------------------------------
# 4. Candidate ligand sets
# ------------------------------------------------------------

sensory_ligands_raw <-
    unique(sensory$ligand)

sensory_ligands <-
    intersect(
        sensory_ligands_raw,
        colnames(ligand_target_matrix)
    )

sensory_missing <-
    setdiff(
        sensory_ligands_raw,
        sensory_ligands
    )

# Broader high-confidence beta-supported ligand pool.
#
# Prefer ligands reproduced as beta-supported in both
# GSE81608 representations.
# ------------------------------------------------------------

if (
    "beta_supported_in_both_representations" %in%
    colnames(broad)
) {

    broad_keep <-
        broad$beta_supported_in_both_representations %in%
        c(TRUE, "TRUE", "True", 1, "1")

    broad_ligands_raw <-
        unique(
            broad$ligand[
                broad_keep
            ]
        )

} else {

    broad_ligands_raw <-
        unique(broad$ligand)
}

broad_ligands <-
    intersect(
        broad_ligands_raw,
        colnames(ligand_target_matrix)
    )

cat("\n============================================\n")
cat("LIGAND OVERLAP\n")
cat("============================================\n")

cat(
    "Sensory-source candidates raw:",
    length(sensory_ligands_raw),
    "\n"
)

cat(
    "Sensory-source candidates in NicheNet:",
    length(sensory_ligands),
    "\n"
)

cat(
    "Sensory-source candidates absent from NicheNet:",
    length(sensory_missing),
    "\n"
)

if (length(sensory_missing) > 0) {
    cat(
        "Absent sensory ligands:",
        paste(
            sensory_missing,
            collapse = ";"
        ),
        "\n"
    )
}

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

if (length(sensory_ligands) < 2) {
    stop(
        "Too few sensory-source ligands represented ",
        "in the NicheNet ligand-target matrix."
    )
}

# ------------------------------------------------------------
# 5. Helper: ligand activity
# ------------------------------------------------------------

run_activity <- function(
    geneset,
    potential_ligands,
    analysis_name
) {

    if (length(geneset) < 2) {
        stop(
            analysis_name,
            ": fewer than 2 receiver genes."
        )
    }

    if (length(potential_ligands) < 1) {
        stop(
            analysis_name,
            ": no candidate ligands."
        )
    }

    activity <-
        nichenetr::predict_ligand_activities(
            geneset = geneset,
            background_expressed_genes =
                background_genes,
            ligand_target_matrix =
                ligand_target_matrix,
            potential_ligands =
                potential_ligands
        )

    activity <- as.data.frame(
        activity,
        stringsAsFactors = FALSE
    )

    activity <- activity[
        order(
            -activity$aupr_corrected,
            -activity$aupr,
            -activity$pearson
        ),
        ,
        drop = FALSE
    ]

    activity$rank_aupr_corrected <-
        seq_len(nrow(activity))

    activity$analysis <-
        analysis_name

    rownames(activity) <- NULL

    activity
}

# ------------------------------------------------------------
# 6. Primary NicheNet analysis
#
# Sensory-source DRG candidates
# x
# 27-gene donor-robust T2D beta-cell program
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PRIMARY NICHE-NET ANALYSIS\n")
cat("============================================\n")

primary_activity <-
    run_activity(
        geneset =
            primary_receiver,
        potential_ligands =
            sensory_ligands,
        analysis_name =
            "Sensory13_PrimaryReceiver27"
    )

write.csv(
    primary_activity,
    paste0(
        "results/nichenet/",
        "Table_P3_NicheNet_SensorySource_",
        "PrimaryReceiver_LigandActivity.csv"
    ),
    row.names = FALSE
)

print(
    primary_activity,
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. High-confidence receiver sensitivity
# ------------------------------------------------------------

cat("\n============================================\n")
cat("HIGH-CONFIDENCE RECEIVER SENSITIVITY\n")
cat("============================================\n")

if (length(core_receiver) >= 2) {

    core_activity <-
        run_activity(
            geneset =
                core_receiver,
            potential_ligands =
                sensory_ligands,
            analysis_name =
                "Sensory13_HighConfidenceCore10"
        )

    write.csv(
        core_activity,
        paste0(
            "results/nichenet/",
            "Table_P3_NicheNet_SensorySource_",
            "HighConfidenceCore_LigandActivity.csv"
        ),
        row.names = FALSE
    )

    print(
        core_activity,
        row.names = FALSE
    )

} else {

    core_activity <- data.frame()

    warning(
        "Too few high-confidence core genes overlap ",
        "NicheNet; core sensitivity skipped."
    )
}

# ------------------------------------------------------------
# 8. Broader HC beta-compatible ligand context
# ------------------------------------------------------------

cat("\n============================================\n")
cat("BROADER HC BETA-COMPATIBLE CONTEXT\n")
cat("============================================\n")

broad_activity <-
    run_activity(
        geneset =
            primary_receiver,
        potential_ligands =
            broad_ligands,
        analysis_name =
            "BroadHC_BetaSupported_PrimaryReceiver27"
    )

write.csv(
    broad_activity,
    paste0(
        "results/nichenet/",
        "Table_P3_NicheNet_BroadHC_",
        "BetaSupported_LigandActivity.csv"
    ),
    row.names = FALSE
)

print(
    head(
        broad_activity,
        20
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 9. Merge primary ranking with DRG/islet evidence
# ------------------------------------------------------------

primary_integrated <- merge(
    sensory,
    primary_activity,
    by.x = "ligand",
    by.y = "test_ligand",
    all.x = TRUE,
    sort = FALSE
)

if (nrow(core_activity) > 0) {

    core_small <- core_activity[
        ,
        c(
            "test_ligand",
            "auroc",
            "aupr",
            "aupr_corrected",
            "pearson",
            "rank_aupr_corrected"
        ),
        drop = FALSE
    ]

    colnames(core_small) <- c(
        "ligand",
        "core_auroc",
        "core_aupr",
        "core_aupr_corrected",
        "core_pearson",
        "core_rank_aupr_corrected"
    )

    primary_integrated <- merge(
        primary_integrated,
        core_small,
        by = "ligand",
        all.x = TRUE,
        sort = FALSE
    )
}

# Restore sensory shortlist ordering where possible.
primary_integrated <-
    primary_integrated[
        match(
            sensory$ligand,
            primary_integrated$ligand
        ),
        ,
        drop = FALSE
    ]

# Then produce a ranked version.
ranked_integrated <-
    primary_integrated[
        order(
            is.na(
                primary_integrated$aupr_corrected
            ),
            -primary_integrated$aupr_corrected
        ),
        ,
        drop = FALSE
    ]

write.csv(
    ranked_integrated,
    paste0(
        "results/nichenet/",
        "Table_P3_Integrated_DRG_Islet_",
        "NicheNet_Prioritization.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 10. Top predicted target links for sensory candidates
# ------------------------------------------------------------

target_link_list <- lapply(
    sensory_ligands,
    function(lig) {

        z <-
            nichenetr::get_weighted_ligand_target_links(
                ligand = lig,
                geneset = primary_receiver,
                ligand_target_matrix =
                    ligand_target_matrix,
                n = length(primary_receiver)
            )

        z <- as.data.frame(
            z,
            stringsAsFactors = FALSE
        )

        if (nrow(z) == 0) {
            return(NULL)
        }

        z
    }
)

target_links <- do.call(
    rbind,
    target_link_list
)

if (is.null(target_links)) {

    target_links <- data.frame()

} else {

    rownames(target_links) <- NULL

    write.csv(
        target_links,
        paste0(
            "results/nichenet/",
            "Table_P3_NicheNet_SensorySource_",
            "PredictedTargetLinks.csv"
        ),
        row.names = FALSE
    )
}

# ------------------------------------------------------------
# 11. Input/model audit tables
# ------------------------------------------------------------

receiver_audit <- data.frame(
    gene = unique(
        c(
            primary_receiver_raw,
            core_receiver_raw
        )
    ),
    primary_receiver =
        unique(
            c(
                primary_receiver_raw,
                core_receiver_raw
            )
        ) %in% primary_receiver_raw,
    high_confidence_core =
        unique(
            c(
                primary_receiver_raw,
                core_receiver_raw
            )
        ) %in% core_receiver_raw,
    in_edgeR_background =
        unique(
            c(
                primary_receiver_raw,
                core_receiver_raw
            )
        ) %in% tested_beta_genes,
    in_NicheNet_target_matrix =
        unique(
            c(
                primary_receiver_raw,
                core_receiver_raw
            )
        ) %in% rownames(
            ligand_target_matrix
        ),
    stringsAsFactors = FALSE
)

write.csv(
    receiver_audit,
    "results/nichenet/GSE81608_NicheNet_ReceiverGene_Audit.csv",
    row.names = FALSE
)

ligand_audit <- data.frame(
    ligand = unique(
        c(
            sensory_ligands_raw,
            broad_ligands_raw
        )
    ),
    sensory_source_candidate =
        unique(
            c(
                sensory_ligands_raw,
                broad_ligands_raw
            )
        ) %in% sensory_ligands_raw,
    broad_HC_beta_supported =
        unique(
            c(
                sensory_ligands_raw,
                broad_ligands_raw
            )
        ) %in% broad_ligands_raw,
    in_NicheNet_ligand_target_matrix =
        unique(
            c(
                sensory_ligands_raw,
                broad_ligands_raw
            )
        ) %in% colnames(
            ligand_target_matrix
        ),
    stringsAsFactors = FALSE
)

write.csv(
    ligand_audit,
    "results/nichenet/DRG_NicheNet_Ligand_Audit.csv",
    row.names = FALSE
)

model_record <- data.frame(
    nichenetr_version =
        as.character(
            packageVersion("nichenetr")
        ),
    organism = "human",
    model =
        "ligand_target_matrix_nsga2r_final.rds",
    source =
        "Zenodo record 7074291",
    local_file =
        lt_file,
    md5 =
        lt_md5,
    model_targets =
        nrow(ligand_target_matrix),
    model_ligands =
        ncol(ligand_target_matrix),
    edgeR_tested_background =
        length(tested_beta_genes),
    NicheNet_background =
        length(background_genes),
    primary_receiver_raw =
        length(primary_receiver_raw),
    primary_receiver_in_model =
        length(primary_receiver),
    core_receiver_raw =
        length(core_receiver_raw),
    core_receiver_in_model =
        length(core_receiver),
    sensory_ligands_raw =
        length(sensory_ligands_raw),
    sensory_ligands_in_model =
        length(sensory_ligands),
    broad_ligands_raw =
        length(broad_ligands_raw),
    broad_ligands_in_model =
        length(broad_ligands),
    interpretation =
        paste(
            "Ligand activity quantifies how well",
            "the NicheNet prior predicts the observed",
            "T2D-associated beta-cell receiver program.",
            "It does not demonstrate physical",
            "DRG-to-islet signaling or causality."
        ),
    stringsAsFactors = FALSE
)

write.csv(
    model_record,
    "results/nichenet/NicheNet_Model_and_Analysis_Record.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 12. Simple activity figure
# ------------------------------------------------------------

plot_df <- primary_activity

plot_df <- plot_df[
    order(
        plot_df$aupr_corrected,
        decreasing = FALSE
    ),
    ,
    drop = FALSE
]

png(
    "figures/P3_NicheNet_SensorySource_LigandActivity.png",
    width = 1800,
    height = 1600,
    res = 220
)

par(
    mar = c(5, 9, 4, 2)
)

barplot(
    plot_df$aupr_corrected,
    names.arg =
        plot_df$test_ligand,
    horiz = TRUE,
    las = 1,
    xlab = "Corrected AUPR",
    main =
        paste0(
            "NicheNet activity of DRG sensory-source ",
            "ligand candidates"
        )
)

dev.off()

# ------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PRIMARY NICHE-NET RANKING\n")
cat("============================================\n")

print(
    primary_activity[
        ,
        intersect(
            c(
                "test_ligand",
                "auroc",
                "aupr",
                "aupr_corrected",
                "pearson",
                "rank_aupr_corrected"
            ),
            colnames(primary_activity)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)

if (nrow(core_activity) > 0) {

    cat("\n============================================\n")
    cat("CORE-SENSITIVITY RANKING\n")
    cat("============================================\n")

    print(
        core_activity[
            ,
            intersect(
                c(
                    "test_ligand",
                    "aupr_corrected",
                    "rank_aupr_corrected"
                ),
                colnames(core_activity)
            ),
            drop = FALSE
        ],
        row.names = FALSE
    )
}

cat("\n============================================\n")
cat("TOP BROAD HC BETA-COMPATIBLE LIGANDS\n")
cat("============================================\n")

print(
    head(
        broad_activity[
            ,
            intersect(
                c(
                    "test_ligand",
                    "aupr_corrected",
                    "rank_aupr_corrected"
                ),
                colnames(broad_activity)
            ),
            drop = FALSE
        ],
        15
    ),
    row.names = FALSE
)

capture.output(
    sessionInfo(),
    file =
        "reproducibility/06_nichenet_DRG_beta_ligand_activity_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 06 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")
