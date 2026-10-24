
cat("====================================================\n")
cat("SCRIPT 07g: PAPER 3 FINAL CANDIDATE EVIDENCE INTEGRATION\n")
cat("====================================================\n\n")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required.")
}

library(data.table)

dir.create(
    "results/integration",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "reproducibility",
    recursive = TRUE,
    showWarnings = FALSE
)

## ============================================================
## INPUT FILES
## ============================================================

f_lr <- paste0(
    "results/islet/",
    "Table_P3_SensorySource_BetaCell_HC_Shortlist.csv"
)

f_primary <- paste0(
    "results/nichenet/",
    "Table_P3_NicheNet_SensorySource_",
    "PrimaryReceiver_LigandActivity.csv"
)

f_core <- paste0(
    "results/nichenet/",
    "Table_P3_NicheNet_SensorySource_",
    "HighConfidenceCore_LigandActivity.csv"
)

f_integrated_nichenet <- paste0(
    "results/nichenet/",
    "Table_P3_Integrated_DRG_Islet_",
    "NicheNet_Prioritization.csv"
)

f_anatomy <- paste0(
    "results/validation/",
    "Table_P3_TraceNSeq_13Ligand_",
    "AnatomicalValidation.csv"
)

f_subtype <- paste0(
    "results/validation/",
    "Table_P3_TraceNSeq_",
    "CandidateExpression_BySubtype.csv"
)

files_required <- c(
    f_lr,
    f_primary,
    f_core,
    f_integrated_nichenet,
    f_anatomy,
    f_subtype
)

missing_files <- files_required[
    !file.exists(files_required)
]

if (length(missing_files) > 0) {
    stop(
        "Missing required files:\n",
        paste(missing_files, collapse = "\n")
    )
}

cat("All required input files: PASS\n\n")

## ============================================================
## 13 SENSORY-SOURCE CANDIDATES
## ============================================================

candidate_map <- data.frame(
    HumanGene = c(
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
    ),
    MouseGene = c(
        "Rph3a",
        "Edil3",
        "Spp1",
        "Sema7a",
        "Ppia",
        "Fgf13",
        "Fgf12",
        "Sst",
        "Cntn2",
        "Lgals8",
        "Rspo2",
        "Flrt3",
        "Pgf"
    ),
    stringsAsFactors = FALSE
)

## ============================================================
## HELPERS
## ============================================================

find_gene_column <- function(x) {

    nms <- names(x)
    low <- tolower(nms)

    preferred <- c(
        "test_ligand",
        "potential_ligand",
        "ligand",
        "human_gene",
        "humangene",
        "gene",
        "genesymbol",
        "source_genesymbol",
        "source"
    )

    for (p in preferred) {

        hit <- which(low == p)

        if (length(hit) > 0) {
            return(nms[hit[1]])
        }
    }

    hit <- grep(
        "ligand|gene",
        low
    )

    if (length(hit) > 0) {
        return(nms[hit[1]])
    }

    stop(
        "Could not identify ligand/gene column in table with columns: ",
        paste(nms, collapse = ", ")
    )
}

get_column <- function(
    x,
    candidates,
    default = NA_real_
) {

    nms <- names(x)
    low <- tolower(nms)

    for (candidate in candidates) {

        hit <- which(
            low == tolower(candidate)
        )

        if (length(hit) > 0) {
            return(x[[hit[1]]])
        }
    }

    rep(default, nrow(x))
}

collapse_table_by_gene <- function(
    x,
    gene_col,
    prefix
) {

    x$HumanGene <- toupper(
        trimws(
            as.character(x[[gene_col]])
        )
    )

    x <- x[
        x$HumanGene %in%
            candidate_map$HumanGene,
        ,
        drop = FALSE
    ]

    if (nrow(x) == 0) {
        stop(
            "No candidate genes found after reading ",
            prefix
        )
    }

    other <- setdiff(
        names(x),
        c(gene_col, "HumanGene")
    )

    ## If already one row per ligand, preserve values directly.
    if (!anyDuplicated(x$HumanGene)) {

        out <- x[
            ,
            c("HumanGene", other),
            drop = FALSE
        ]

        names(out)[-1] <- paste0(
            prefix,
            names(out)[-1]
        )

        return(out)
    }

    ## If multiple receptor rows exist per ligand,
    ## collapse values without discarding information.
    dt <- as.data.table(x)

    collapsed <- dt[
        ,
        lapply(
            .SD,
            function(v) {

                v <- v[
                    !is.na(v)
                ]

                if (length(v) == 0) {
                    return(NA_character_)
                }

                vals <- unique(
                    as.character(v)
                )

                paste(
                    vals,
                    collapse = "; "
                )
            }
        ),
        by = HumanGene,
        .SDcols = other
    ]

    collapsed <- as.data.frame(
        collapsed
    )

    names(collapsed)[-1] <- paste0(
        prefix,
        names(collapsed)[-1]
    )

    collapsed
}

standardize_activity <- function(
    x,
    prefix
) {

    gene_col <- find_gene_column(x)

    out <- data.frame(
        HumanGene = toupper(
            trimws(
                as.character(x[[gene_col]])
            )
        ),
        stringsAsFactors = FALSE
    )

    out[[paste0(prefix, "AUPR")]] <-
        as.numeric(
            get_column(
                x,
                c("aupr")
            )
        )

    out[[paste0(prefix, "AUPR_corrected")]] <-
        as.numeric(
            get_column(
                x,
                c(
                    "aupr_corrected",
                    "corrected_aupr"
                )
            )
        )

    out[[paste0(prefix, "Pearson")]] <-
        as.numeric(
            get_column(
                x,
                c("pearson")
            )
        )

    out[[paste0(prefix, "Rank")]] <-
        as.numeric(
            get_column(
                x,
                c(
                    "rank_aupr_corrected",
                    "rank_corrected_aupr",
                    "rank"
                )
            )
        )

    out <- out[
        out$HumanGene %in%
            candidate_map$HumanGene,
        ,
        drop = FALSE
    ]

    out <- out[
        !duplicated(out$HumanGene),
        ,
        drop = FALSE
    ]

    out
}

## ============================================================
## READ HUMAN DRG / BETA-CELL LR SHORTLIST
## ============================================================

lr_raw <- read.csv(
    f_lr,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat(
    "Sensory-source HC shortlist:",
    nrow(lr_raw),
    "rows\n"
)

lr_gene_col <- find_gene_column(
    lr_raw
)

cat(
    "Ligand column identified as:",
    lr_gene_col,
    "\n"
)

lr <- collapse_table_by_gene(
    lr_raw,
    lr_gene_col,
    "LR_"
)

cat(
    "Unique candidate ligands in LR table:",
    nrow(lr),
    "\n\n"
)

## ============================================================
## READ NICHENET ACTIVITIES
## ============================================================

primary_raw <- read.csv(
    f_primary,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

core_raw <- read.csv(
    f_core,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

primary <- standardize_activity(
    primary_raw,
    "NicheNet_Primary_"
)

core <- standardize_activity(
    core_raw,
    "NicheNet_Core_"
)

cat(
    "NicheNet primary evaluable sensory ligands:",
    nrow(primary),
    "\n"
)

cat(
    "NicheNet core evaluable sensory ligands:",
    nrow(core),
    "\n\n"
)

## ============================================================
## TRACE-N-SEQ ANATOMICAL VALIDATION
## ============================================================

anatomy <- read.csv(
    f_anatomy,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!"Gene" %in% names(anatomy)) {
    stop(
        "Trace-n-Seq anatomical table lacks Gene column."
    )
}

anatomy$MouseGene <- as.character(
    anatomy$Gene
)

anatomy$HumanGene <- toupper(
    anatomy$MouseGene
)

keep_anatomy <- c(
    "HumanGene",
    "MouseGene",
    "DRG_Batch1_PctDetected",
    "DRG_Batch2_PctDetected",
    "DRG_Combined_PctDetected",
    "NG_PctDetected",
    "DRG_Batch1_MeanLogCPM",
    "DRG_Batch2_MeanLogCPM",
    "NG_MeanLogCPM",
    "Detected_Both_DRG_Batches",
    "Detected_DRG_and_NG"
)

keep_anatomy <- intersect(
    keep_anatomy,
    names(anatomy)
)

anatomy <- anatomy[
    ,
    keep_anatomy,
    drop = FALSE
]

cat(
    "Trace-n-Seq candidates:",
    nrow(anatomy),
    "\n\n"
)

## ============================================================
## TRACE-N-SEQ SENSORY SUBTYPE INTEGRATION
## ============================================================

subtype <- read.csv(
    f_subtype,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

needed_subtype <- c(
    "Batch",
    "Ganglion",
    "CellType",
    "Gene",
    "N_cells",
    "Pct_detected",
    "Mean_log1p_CPM"
)

if (!all(
    needed_subtype %in% names(subtype)
)) {
    stop(
        "Subtype table does not contain expected columns."
    )
}

subtype$HumanGene <- toupper(
    as.character(
        subtype$Gene
    )
)

subtype <- subtype[
    subtype$HumanGene %in%
        candidate_map$HumanGene,
    ,
    drop = FALSE
]

sub_drg <- as.data.table(
    subtype[
        subtype$Ganglion == "DRG",
        ,
        drop = FALSE
    ]
)

subtype_combined <- sub_drg[
    ,
    .(
        N_cells = sum(
            N_cells,
            na.rm = TRUE
        ),

        Pct_detected = sum(
            N_cells *
            Pct_detected,
            na.rm = TRUE
        ) /
        sum(
            N_cells,
            na.rm = TRUE
        ),

        Mean_log1p_CPM = sum(
            N_cells *
            Mean_log1p_CPM,
            na.rm = TRUE
        ) /
        sum(
            N_cells,
            na.rm = TRUE
        )
    ),
    by = .(
        HumanGene,
        CellType
    )
]

setorder(
    subtype_combined,
    HumanGene,
    -Pct_detected,
    -Mean_log1p_CPM
)

dominant_subtype <- subtype_combined[
    ,
    .SD[1],
    by = HumanGene
]

dominant_subtype <- as.data.frame(
    dominant_subtype
)

names(dominant_subtype)[
    names(dominant_subtype) ==
        "CellType"
] <- "DRG_Dominant_SensorySubtype"

names(dominant_subtype)[
    names(dominant_subtype) ==
        "Pct_detected"
] <- "DRG_DominantSubtype_PctDetected"

names(dominant_subtype)[
    names(dominant_subtype) ==
        "Mean_log1p_CPM"
] <- "DRG_DominantSubtype_MeanLogCPM"

names(dominant_subtype)[
    names(dominant_subtype) ==
        "N_cells"
] <- "DRG_DominantSubtype_Ncells"

write.csv(
    subtype_combined,
    paste0(
        "results/integration/",
        "Table_P3_Final_SensorySubtypeEvidence.csv"
    ),
    row.names = FALSE
)

## ============================================================
## BIOLOGICAL CURATION
##
## IMPORTANT:
## This is not a mathematical score.
## It records biological interpretation rules decided
## independently of NicheNet ranking.
## ============================================================

curation <- data.frame(

    HumanGene = candidate_map$HumanGene,

    Signaling_Curation = c(
        "Intracellular/synaptic-vesicle-associated; not treated as a conventional extracellular ligand",
        "Extracellular matrix/integrin ligand",
        "Secreted extracellular signaling protein",
        "GPI-anchored membrane signaling protein",
        "Non-classically secreted stress-associated protein; interpret extracellular signaling cautiously",
        "Intracellular FHF; not treated as a canonical secreted FGFR ligand",
        "Intracellular FHF; not treated as a canonical secreted FGFR ligand",
        "Secreted neuropeptide",
        "Cell-surface adhesion molecule; contact-dependent interpretation",
        "Non-classically secreted extracellular lectin",
        "Secreted extracellular signaling protein",
        "Transmembrane adhesion/signaling protein; contact-dependent interpretation",
        "Extracellular growth factor"
    ),

    Headline_Eligible = c(
        FALSE,  # RPH3A
        TRUE,   # EDIL3
        TRUE,   # SPP1
        TRUE,   # SEMA7A
        FALSE,  # PPIA
        FALSE,  # FGF13
        FALSE,  # FGF12
        TRUE,   # SST
        TRUE,   # CNTN2
        TRUE,   # LGALS8
        TRUE,   # RSPO2
        TRUE,   # FLRT3
        TRUE    # PGF
    ),

    Original_GSE295206_Highlight = c(
        FALSE,
        FALSE,
        TRUE,   # SPP1
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE
    ),

    Evidence_Class = c(
        "Mechanistic-caveat candidate",
        "Lead integrated candidate",
        "Reference/established axis",
        "Lead integrated candidate",
        "Secondary/non-classical candidate",
        "Mechanistic-caveat candidate",
        "Mechanistic-caveat candidate",
        "Secondary subtype-restricted candidate",
        "Secondary candidate",
        "Secondary candidate",
        "Secondary candidate",
        "NicheNet-prioritized candidate",
        "NicheNet-prioritized candidate"
    ),

    stringsAsFactors = FALSE
)

## ============================================================
## MERGE ALL EVIDENCE
## ============================================================

final <- candidate_map

final <- merge(
    final,
    lr,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

final <- merge(
    final,
    primary,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

final <- merge(
    final,
    core,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

final <- merge(
    final,
    anatomy,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

final <- merge(
    final,
    dominant_subtype,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

final <- merge(
    final,
    curation,
    by = "HumanGene",
    all.x = TRUE,
    sort = FALSE
)

## Restore original candidate order.
final$Candidate_Order <- match(
    final$HumanGene,
    candidate_map$HumanGene
)

final <- final[
    order(final$Candidate_Order),
    ,
    drop = FALSE
]

## NicheNet evaluability.
final$NicheNet_Evaluable <-
    !is.na(
        final$NicheNet_Primary_AUPR_corrected
    )

## Positive corrected AUPR in both receiver analyses.
final$NicheNet_Positive_Both <- (
    !is.na(
        final$NicheNet_Primary_AUPR_corrected
    ) &
    !is.na(
        final$NicheNet_Core_AUPR_corrected
    ) &
    final$NicheNet_Primary_AUPR_corrected > 0 &
    final$NicheNet_Core_AUPR_corrected > 0
)

## Anatomical replication.
final$TraceNSeq_DRG_Replicated <- (
    !is.na(
        final$DRG_Batch1_PctDetected
    ) &
    !is.na(
        final$DRG_Batch2_PctDetected
    ) &
    final$DRG_Batch1_PctDetected > 0 &
    final$DRG_Batch2_PctDetected > 0
)

## Useful descriptive categories.
final$TraceNSeq_DRG_ExpressionClass <-
    ifelse(
        is.na(
            final$DRG_Combined_PctDetected
        ),
        NA,
        ifelse(
            final$DRG_Combined_PctDetected >= 75,
            "Broad (>=75%)",
            ifelse(
                final$DRG_Combined_PctDetected >= 25,
                "Intermediate (25-<75%)",
                "Restricted (<25%)"
            )
        )
    )

## ============================================================
## EVIDENCE-BASED DISPLAY ORDER
##
## This is deliberately NOT an additive point score.
## Candidates are grouped by biological interpretation.
## ============================================================

display_order <- c(
    "SEMA7A",
    "EDIL3",
    "FLRT3",
    "PGF",
    "SPP1",
    "RSPO2",
    "SST",
    "LGALS8",
    "CNTN2",
    "PPIA",
    "RPH3A",
    "FGF13",
    "FGF12"
)

final$Display_Order <- match(
    final$HumanGene,
    display_order
)

final <- final[
    order(final$Display_Order),
    ,
    drop = FALSE
]

## ============================================================
## SAVE COMPLETE TABLE
## ============================================================

write.csv(
    final,
    paste0(
        "results/integration/",
        "Table_P3_Final_13Candidate_",
        "EvidenceIntegration.csv"
    ),
    row.names = FALSE
)

## ============================================================
## MAIN MANUSCRIPT CANDIDATES
## ============================================================

main_genes <- c(
    "SEMA7A",
    "EDIL3",
    "FLRT3",
    "PGF",
    "SPP1"
)

main <- final[
    final$HumanGene %in% main_genes,
    ,
    drop = FALSE
]

main$Main_Interpretation <- c(
    "Lead integrated cross-tissue candidate",
    "Lead integrated cross-tissue candidate",
    "NicheNet-prioritized cross-tissue candidate",
    "NicheNet-prioritized cross-tissue candidate",
    "Reference axis highlighted in original diabetic DRG study"
)[
    match(
        main$HumanGene,
        main_genes
    )
]

main$Main_Order <- match(
    main$HumanGene,
    main_genes
)

main <- main[
    order(main$Main_Order),
    ,
    drop = FALSE
]

write.csv(
    main,
    paste0(
        "results/integration/",
        "Table_P3_Final_MainCandidate_Hierarchy.csv"
    ),
    row.names = FALSE
)

## ============================================================
## COMPACT CONSOLE TABLE
## ============================================================

compact_cols <- c(
    "HumanGene",
    "Evidence_Class",
    "NicheNet_Primary_Rank",
    "NicheNet_Primary_AUPR_corrected",
    "NicheNet_Core_Rank",
    "NicheNet_Core_AUPR_corrected",
    "DRG_Batch1_PctDetected",
    "DRG_Batch2_PctDetected",
    "DRG_Combined_PctDetected",
    "NG_PctDetected",
    "DRG_Dominant_SensorySubtype",
    "DRG_DominantSubtype_PctDetected",
    "NicheNet_Evaluable",
    "NicheNet_Positive_Both",
    "TraceNSeq_DRG_Replicated",
    "Original_GSE295206_Highlight"
)

compact_cols <- intersect(
    compact_cols,
    names(final)
)

cat("\n====================================================\n")
cat("FINAL 13-CANDIDATE EVIDENCE INTEGRATION\n")
cat("====================================================\n")

print(
    final[
        ,
        compact_cols,
        drop = FALSE
    ],
    row.names = FALSE,
    digits = 5
)

cat("\n====================================================\n")
cat("MAIN MANUSCRIPT CANDIDATES\n")
cat("====================================================\n")

main_compact <- c(
    "HumanGene",
    "Main_Interpretation",
    "NicheNet_Primary_Rank",
    "NicheNet_Core_Rank",
    "DRG_Combined_PctDetected",
    "NG_PctDetected",
    "DRG_Dominant_SensorySubtype",
    "DRG_DominantSubtype_PctDetected",
    "Signaling_Curation"
)

main_compact <- intersect(
    main_compact,
    names(main)
)

print(
    main[
        ,
        main_compact,
        drop = FALSE
    ],
    row.names = FALSE,
    digits = 5
)

## ============================================================
## AUDIT
## ============================================================

cat("\n====================================================\n")
cat("AUDIT\n")
cat("====================================================\n")

cat(
    "Candidates expected:",
    nrow(candidate_map),
    "\n"
)

cat(
    "Candidates integrated:",
    nrow(final),
    "\n"
)

cat(
    "NicheNet evaluable:",
    sum(
        final$NicheNet_Evaluable,
        na.rm = TRUE
    ),
    "/",
    nrow(final),
    "\n"
)

cat(
    "Positive corrected AUPR in BOTH receiver analyses:",
    sum(
        final$NicheNet_Positive_Both,
        na.rm = TRUE
    ),
    "/",
    nrow(final),
    "\n"
)

cat(
    "Detected in BOTH pancreas-traced DRG batches:",
    sum(
        final$TraceNSeq_DRG_Replicated,
        na.rm = TRUE
    ),
    "/",
    nrow(final),
    "\n"
)

cat(
    "Main manuscript candidates:",
    paste(
        main$HumanGene,
        collapse = ", "
    ),
    "\n"
)

## ============================================================
## SAVE RDS + SESSION INFO
## ============================================================

saveRDS(
    list(
        candidate_map = candidate_map,
        final_integration = final,
        main_candidates = main,
        subtype_combined = subtype_combined,
        biological_curation = curation,
        nichenet_primary = primary,
        nichenet_core = core
    ),
    paste0(
        "data/processed/",
        "Paper3_final_candidate_",
        "evidence_integration.rds"
    )
)

writeLines(
    capture.output(
        sessionInfo()
    ),
    paste0(
        "reproducibility/",
        "07g_Paper3_candidate_",
        "integration_sessionInfo.txt"
    )
)

cat("\nOutput files:\n")

cat(
    "  results/integration/",
    "Table_P3_Final_13Candidate_EvidenceIntegration.csv\n",
    sep = ""
)

cat(
    "  results/integration/",
    "Table_P3_Final_MainCandidate_Hierarchy.csv\n",
    sep = ""
)

cat(
    "  results/integration/",
    "Table_P3_Final_SensorySubtypeEvidence.csv\n",
    sep = ""
)

cat(
    "  data/processed/",
    "Paper3_final_candidate_evidence_integration.rds\n",
    sep = ""
)

cat("\n====================================================\n")
cat("SCRIPT 07g COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

