
options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", mustWork = TRUE)
proc_dir <- file.path(project_dir, "data", "processed")
islet_dir <- file.path(project_dir, "results", "islet")
repro_dir <- file.path(project_dir, "reproducibility")

dir.create(islet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 04f: GSE81608 1600-CELL RECEPTOR SENSITIVITY\n")
cat("============================================\n")

expr_file <- file.path(
    proc_dir,
    "GSE81608_RPKM_1600cells_SYMBOL.rds"
)

meta_file <- file.path(
    proc_dir,
    "GSE81608_QC_1492_cell_metadata.csv"
)

hc_file <- file.path(
    islet_dir,
    "Table_P3_HighConfidence_LR_Islet_Compatibility.csv"
)

nearby_file <- file.path(
    islet_dir,
    "Table_P3_SensorySource_BetaCell_HC_Shortlist.csv"
)

required_files <- c(
    expr_file,
    meta_file,
    hc_file,
    nearby_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files)) {
    stop(
        "Missing required files:\n",
        paste(missing_files, collapse = "\n")
    )
}

# ============================================================
# EXACT SAME THRESHOLDS AS SCRIPT 04d
# ============================================================

MIN_CELLS_PER_DONOR <- 5L

PRIMARY_CELL_FRACTION <- 0.10
STRICT_CELL_FRACTION <- 0.20

MIN_EVALUABLE_DONORS <- 3L

PRIMARY_DONOR_FRACTION <- 0.50
STRONG_DONOR_FRACTION <- 0.75


to_bool <- function(z) {

    if (is.logical(z)) {
        return(z)
    }

    tolower(
        trimws(as.character(z))
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
}


collapse_unique <- function(z) {

    z <- as.character(z)

    z <- sort(
        unique(
            z[
                !is.na(z) &
                nzchar(z)
            ]
        )
    )

    if (!length(z)) {
        return(NA_character_)
    }

    paste(
        z,
        collapse = ";"
    )
}

# ============================================================
# 1. Load current full 1600-cell GEO representation
# ============================================================

expr <- readRDS(expr_file)

meta <- read.csv(
    meta_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
}

cat("\nExpression matrix:",
    nrow(expr), "genes x",
    ncol(expr), "cells\n"
)

cat(
    "Metadata rows:",
    nrow(meta),
    "\n"
)

if (
    ncol(expr) != 1600L ||
    nrow(meta) != 1600L
) {
    stop(
        "Expected exactly 1600 cells."
    )
}

needed_meta <- c(
    "donor_id",
    "condition",
    "cell_subtype"
)

if (!all(needed_meta %in% names(meta))) {

    stop(
        "Missing metadata columns: ",
        paste(
            setdiff(
                needed_meta,
                names(meta)
            ),
            collapse = ", "
        )
    )
}

# ============================================================
# 2. Align metadata to Sample_N columns
# ============================================================

if ("matrix_column" %in% names(meta)) {

    idx <- match(
        colnames(expr),
        meta$matrix_column
    )

} else if ("sample_number" %in% names(meta)) {

    idx <- match(
        colnames(expr),
        paste0(
            "Sample_",
            meta$sample_number
        )
    )

} else {

    sample_num <- suppressWarnings(
        as.integer(
            sub(
                "^Sample_",
                "",
                colnames(expr)
            )
        )
    )

    if (
        anyNA(sample_num) ||
        !identical(
            sort(sample_num),
            seq_len(ncol(expr))
        )
    ) {
        stop(
            "Cannot safely infer Sample_N metadata mapping."
        )
    }

    idx <- sample_num
}

if (anyNA(idx)) {
    stop(
        "Could not align all expression columns to metadata."
    )
}

meta <- meta[
    idx,
    ,
    drop = FALSE
]

rownames(meta) <- colnames(expr)

cat(
    "Expression/metadata alignment: PASS\n"
)

# ============================================================
# 3. Gene-symbol validation
# ============================================================

if (
    is.null(rownames(expr)) ||
    any(!nzchar(rownames(expr)))
) {
    stop(
        "Expression matrix lacks gene-symbol row names."
    )
}

markers <- c(
    "INS",
    "GCG",
    "SST",
    "PPY"
)

marker_present <- markers[
    markers %in% rownames(expr)
]

cat(
    "Endocrine marker symbols present:",
    paste(
        marker_present,
        collapse = ";"
    ),
    "\n"
)

if (length(marker_present) < 3L) {
    stop(
        "Gene row names do not appear to be gene symbols."
    )
}

cat(
    "Duplicated gene symbols:",
    sum(
        duplicated(
            rownames(expr)
        )
    ),
    "\n"
)

# Only zero/non-zero detection is used.
# Summing duplicate symbols therefore preserves detection status.

if (anyDuplicated(rownames(expr))) {

    expr <- rowsum(
        expr,
        group = rownames(expr),
        reorder = FALSE
    )

    cat(
        "Duplicated symbols collapsed for detection analysis.\n"
    )
}

# ============================================================
# 4. Load frozen high-confidence LR network
# ============================================================

hc <- read.csv(
    hc_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_hc <- c(
    "ligand",
    "receptor",
    "candidate_tier",
    "logFC",
    "FDR",
    "enriched_region",
    "high_confidence_sensitivity",
    "beta_ND_primary_compatible",
    "beta_T2D_primary_compatible",
    "beta_ND_strong_compatible",
    "beta_T2D_strong_compatible",
    "alpha_ND_primary_compatible",
    "alpha_T2D_primary_compatible"
)

missing_hc <- setdiff(
    required_hc,
    names(hc)
)

if (length(missing_hc)) {

    stop(
        "04e table missing columns:\n",
        paste(
            missing_hc,
            collapse = "\n"
        )
    )
}

bool_cols <- c(
    "high_confidence_sensitivity",
    "beta_ND_primary_compatible",
    "beta_T2D_primary_compatible",
    "beta_ND_strong_compatible",
    "beta_T2D_strong_compatible",
    "alpha_ND_primary_compatible",
    "alpha_T2D_primary_compatible"
)

for (cc in bool_cols) {
    hc[[cc]] <- to_bool(
        hc[[cc]]
    )
}

hc <- hc[
    hc$high_confidence_sensitivity,
    ,
    drop = FALSE
]

hc <- hc[
    !duplicated(
        paste(
            hc$ligand,
            hc$receptor,
            sep = "|||"
        )
    ),
    ,
    drop = FALSE
]

cat(
    "\nHigh-confidence LR pairs:",
    nrow(hc),
    "\n"
)

# ============================================================
# 5. Receptor identifier validation
# ============================================================

candidate_receptors <- sort(
    unique(
        hc$receptor
    )
)

receptors_present <- intersect(
    candidate_receptors,
    rownames(expr)
)

receptors_absent <- setdiff(
    candidate_receptors,
    rownames(expr)
)

cat("\n============================================\n")
cat("RECEPTOR IDENTIFIER VALIDATION\n")
cat("============================================\n")

cat(
    "HC candidate receptors:",
    length(candidate_receptors),
    "\n"
)

cat(
    "Present in 1600-cell matrix:",
    length(receptors_present),
    "\n"
)

cat(
    "Absent from 1600-cell matrix:",
    length(receptors_absent),
    "\n"
)

if (length(receptors_present) < 100L) {

    stop(
        paste(
            "Unexpectedly low receptor overlap.",
            "Inspect gene identifiers."
        )
    )
}

write.csv(
    data.frame(
        receptor = receptors_absent
    ),
    file.path(
        islet_dir,
        "GSE81608_1600_HC_receptors_absent_from_matrix.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 6. Alpha/beta cells
# ============================================================

meta$cell_type <- tolower(
    trimws(
        as.character(
            meta$cell_subtype
        )
    )
)

meta$condition <- trimws(
    as.character(
        meta$condition
    )
)

meta$condition[
    tolower(meta$condition) %in%
        c(
            "t2d",
            "type 2 diabetes",
            "type2 diabetes"
        )
] <- "T2D"

meta$condition[
    tolower(meta$condition) %in%
        c(
            "non-diabetic",
            "non diabetic",
            "nont2d",
            "non t2d"
        )
] <- "non-diabetic"

keep_cells <- meta$cell_type %in%
    c(
        "alpha",
        "beta"
    )

expr_ab <- expr[
    receptors_present,
    keep_cells,
    drop = FALSE
]

meta_ab <- meta[
    keep_cells,
    ,
    drop = FALSE
]

stopifnot(
    identical(
        colnames(expr_ab),
        rownames(meta_ab)
    )
)

cat("\n============================================\n")
cat("1600-CELL ALPHA/BETA REPRESENTATION\n")
cat("============================================\n")

cat(
    "Alpha/beta cells:",
    ncol(expr_ab),
    "\n"
)

print(
    table(
        meta_ab$condition,
        meta_ab$cell_type
    )
)

# ============================================================
# 7. Donor-level receptor detection
# ============================================================

group_key <- paste(
    meta_ab$donor_id,
    meta_ab$condition,
    meta_ab$cell_type,
    sep = "|||"
)

groups <- split(
    seq_len(nrow(meta_ab)),
    group_key
)

donor_results <- vector(
    "list",
    length(receptors_present) *
        length(groups)
)

counter <- 1L

cat(
    "\nCalculating donor-level receptor detection...\n"
)

for (receptor in receptors_present) {

    values <- expr_ab[
        receptor,
        ,
        drop = TRUE
    ]

    for (g in names(groups)) {

        ii <- groups[[g]]

        x <- values[ii]

        n_cells <- length(x)

        n_detected <- sum(
            x > 0,
            na.rm = TRUE
        )

        frac <- n_detected /
            n_cells

        evaluable <-
            n_cells >=
            MIN_CELLS_PER_DONOR

        donor_results[[counter]] <-
            data.frame(

                receptor = receptor,

                donor_id =
                    meta_ab$donor_id[
                        ii[1]
                    ],

                condition =
                    meta_ab$condition[
                        ii[1]
                    ],

                cell_type =
                    meta_ab$cell_type[
                        ii[1]
                    ],

                n_cells =
                    n_cells,

                n_detected =
                    n_detected,

                fraction_detected =
                    frac,

                evaluable =
                    evaluable,

                support_any =
                    evaluable &&
                    n_detected >= 1L,

                support_10pct =
                    evaluable &&
                    frac >=
                    PRIMARY_CELL_FRACTION,

                support_20pct =
                    evaluable &&
                    frac >=
                    STRICT_CELL_FRACTION,

                stringsAsFactors = FALSE
            )

        counter <- counter + 1L
    }
}

donor_support <- do.call(
    rbind,
    donor_results
)

rownames(donor_support) <- NULL

write.csv(
    donor_support,
    file.path(
        islet_dir,
        "GSE81608_1600_Receptor_DonorLevel_Detection.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 8. Summarize across donors
# ============================================================

summary_key <- paste(
    donor_support$receptor,
    donor_support$condition,
    donor_support$cell_type,
    sep = "|||"
)

summary_groups <- split(
    seq_len(nrow(donor_support)),
    summary_key
)

support_summary <- do.call(
    rbind,
    lapply(
        summary_groups,
        function(ii) {

            z <- donor_support[
                ii,
                ,
                drop = FALSE
            ]

            ev <- z$evaluable

            n_eval <- sum(ev)

            n_any <- sum(
                z$support_any & ev
            )

            n_10 <- sum(
                z$support_10pct & ev
            )

            n_20 <- sum(
                z$support_20pct & ev
            )

            frac_any <- if (n_eval > 0) {
                n_any / n_eval
            } else {
                NA_real_
            }

            frac_10 <- if (n_eval > 0) {
                n_10 / n_eval
            } else {
                NA_real_
            }

            frac_20 <- if (n_eval > 0) {
                n_20 / n_eval
            } else {
                NA_real_
            }

            data.frame(

                receptor =
                    z$receptor[1],

                condition =
                    z$condition[1],

                cell_type =
                    z$cell_type[1],

                n_evaluable_donors =
                    n_eval,

                n_donors_any =
                    n_any,

                n_donors_10pct =
                    n_10,

                n_donors_20pct =
                    n_20,

                donor_fraction_any =
                    frac_any,

                donor_fraction_10pct =
                    frac_10,

                donor_fraction_20pct =
                    frac_20,

                median_donor_detection_fraction =
                    if (n_eval > 0) {
                        median(
                            z$fraction_detected[ev],
                            na.rm = TRUE
                        )
                    } else {
                        NA_real_
                    },

                primary_compatible =
                    n_eval >=
                    MIN_EVALUABLE_DONORS &&
                    !is.na(frac_10) &&
                    frac_10 >=
                    PRIMARY_DONOR_FRACTION,

                strong_compatible =
                    n_eval >=
                    MIN_EVALUABLE_DONORS &&
                    !is.na(frac_20) &&
                    frac_20 >=
                    STRONG_DONOR_FRACTION,

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(support_summary) <- NULL

write.csv(
    support_summary,
    file.path(
        islet_dir,
        "GSE81608_1600_Receptor_DonorSupport_Summary.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 9. Evaluable donors
# ============================================================

eval_summary <- unique(
    donor_support[
        ,
        c(
            "donor_id",
            "condition",
            "cell_type",
            "n_cells",
            "evaluable"
        )
    ]
)

eval_table <- aggregate(
    evaluable ~
        condition +
        cell_type,
    data = eval_summary,
    FUN = sum
)

cat("\n============================================\n")
cat("EVALUABLE DONORS: 1600-CELL SENSITIVITY\n")
cat("============================================\n")

print(
    eval_table,
    row.names = FALSE
)

# ============================================================
# 10. Helper for receptor support
# ============================================================

get_support <- function(
    condition_value,
    celltype_value,
    prefix
) {

    z <- support_summary[
        support_summary$condition ==
            condition_value &
        support_summary$cell_type ==
            celltype_value,
        ,
        drop = FALSE
    ]

    keep <- c(
        "receptor",
        "n_evaluable_donors",
        "n_donors_10pct",
        "n_donors_20pct",
        "donor_fraction_10pct",
        "donor_fraction_20pct",
        "median_donor_detection_fraction",
        "primary_compatible",
        "strong_compatible"
    )

    z <- z[
        ,
        keep,
        drop = FALSE
    ]

    names(z)[-1] <- paste0(
        prefix,
        "_",
        names(z)[-1]
    )

    z
}

beta_ND_1600 <- get_support(
    "non-diabetic",
    "beta",
    "sens1600_beta_ND"
)

beta_T2D_1600 <- get_support(
    "T2D",
    "beta",
    "sens1600_beta_T2D"
)

alpha_ND_1600 <- get_support(
    "non-diabetic",
    "alpha",
    "sens1600_alpha_ND"
)

alpha_T2D_1600 <- get_support(
    "T2D",
    "alpha",
    "sens1600_alpha_T2D"
)

# ============================================================
# 11. Attach 1600-cell support to 04e HC pairs
# ============================================================

pair_compare <- hc

support_tables <- list(
    beta_ND_1600,
    beta_T2D_1600,
    alpha_ND_1600,
    alpha_T2D_1600
)

for (z in support_tables) {

    pair_compare <- merge(
        pair_compare,
        z,
        by = "receptor",
        all.x = TRUE,
        sort = FALSE
    )
}

sens_bool_cols <- grep(
    "^sens1600_.*_(primary|strong)_compatible$",
    names(pair_compare),
    value = TRUE
)

for (cc in sens_bool_cols) {

    pair_compare[[cc]][
        is.na(
            pair_compare[[cc]]
        )
    ] <- FALSE
}

# ============================================================
# 12. Primary 1492 vs sensitivity 1600 classifications
# ============================================================

pair_compare$primary1492_beta_both <-
    pair_compare$beta_ND_primary_compatible &
    pair_compare$beta_T2D_primary_compatible

pair_compare$primary1492_beta_strong_both <-
    pair_compare$beta_ND_strong_compatible &
    pair_compare$beta_T2D_strong_compatible

pair_compare$sens1600_beta_both <-
    pair_compare$sens1600_beta_ND_primary_compatible &
    pair_compare$sens1600_beta_T2D_primary_compatible

pair_compare$sens1600_beta_strong_both <-
    pair_compare$sens1600_beta_ND_strong_compatible &
    pair_compare$sens1600_beta_T2D_strong_compatible

pair_compare$primary1492_alpha_both <-
    pair_compare$alpha_ND_primary_compatible &
    pair_compare$alpha_T2D_primary_compatible

pair_compare$sens1600_alpha_both <-
    pair_compare$sens1600_alpha_ND_primary_compatible &
    pair_compare$sens1600_alpha_T2D_primary_compatible

pair_compare$beta_primary_class_agreement <-
    pair_compare$primary1492_beta_both ==
    pair_compare$sens1600_beta_both

pair_compare$beta_strong_class_agreement <-
    pair_compare$primary1492_beta_strong_both ==
    pair_compare$sens1600_beta_strong_both

pair_compare$beta_primary_reproduced <-
    pair_compare$primary1492_beta_both &
    pair_compare$sens1600_beta_both

pair_compare$beta_strong_reproduced <-
    pair_compare$primary1492_beta_strong_both &
    pair_compare$sens1600_beta_strong_both

write.csv(
    pair_compare,
    file.path(
        islet_dir,
        "Table_P3_HC_LR_1492_vs_1600_Sensitivity.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 13. Ligand-level sensitivity comparison
# ============================================================

ligand_groups <- split(
    seq_len(nrow(pair_compare)),
    pair_compare$ligand
)

ligand_compare <- do.call(
    rbind,
    lapply(
        names(ligand_groups),
        function(g) {

            z <- pair_compare[
                ligand_groups[[g]],
                ,
                drop = FALSE
            ]

            p1492 <- z$receptor[
                z$primary1492_beta_both
            ]

            s1600 <- z$receptor[
                z$sens1600_beta_both
            ]

            reproduced <- intersect(
                unique(p1492),
                unique(s1600)
            )

            p1492_strong <- z$receptor[
                z$primary1492_beta_strong_both
            ]

            s1600_strong <- z$receptor[
                z$sens1600_beta_strong_both
            ]

            reproduced_strong <- intersect(
                unique(p1492_strong),
                unique(s1600_strong)
            )

            data.frame(

                ligand =
                    g,

                candidate_tier =
                    z$candidate_tier[1],

                logFC =
                    z$logFC[1],

                FDR =
                    z$FDR[1],

                enriched_region =
                    z$enriched_region[1],

                n_HC_edges =
                    nrow(z),

                n_primary1492_beta_both_receptors =
                    length(
                        unique(p1492)
                    ),

                primary1492_beta_both_receptors =
                    collapse_unique(
                        p1492
                    ),

                n_sens1600_beta_both_receptors =
                    length(
                        unique(s1600)
                    ),

                sens1600_beta_both_receptors =
                    collapse_unique(
                        s1600
                    ),

                n_reproduced_beta_both_receptors =
                    length(
                        reproduced
                    ),

                reproduced_beta_both_receptors =
                    collapse_unique(
                        reproduced
                    ),

                n_reproduced_strong_beta_receptors =
                    length(
                        reproduced_strong
                    ),

                reproduced_strong_beta_receptors =
                    collapse_unique(
                        reproduced_strong
                    ),

                primary1492_beta_supported =
                    length(p1492) > 0,

                sens1600_beta_supported =
                    length(s1600) > 0,

                beta_supported_in_both_representations =
                    length(reproduced) > 0,

                strong_beta_supported_in_both_representations =
                    length(
                        reproduced_strong
                    ) > 0,

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(ligand_compare) <- NULL

tier_rank <- match(
    ligand_compare$candidate_tier,
    c(
        "Tier_A_strict",
        "Tier_B_LODO_robust",
        "Tier_C_directionally_robust"
    )
)

tier_rank[
    is.na(tier_rank)
] <- 99L

ligand_compare <- ligand_compare[
    order(
        !ligand_compare$
            strong_beta_supported_in_both_representations,
        !ligand_compare$
            beta_supported_in_both_representations,
        tier_rank,
        ligand_compare$FDR
    ),
    ,
    drop = FALSE
]

write.csv(
    ligand_compare,
    file.path(
        islet_dir,
        "Table_P3_HC_Ligand_1492_vs_1600_Sensitivity.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 14. Exact 13 nearby-neuron candidates
# ============================================================

nearby <- read.csv(
    nearby_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

nearby_primary_cols <- intersect(
    c(
        "ligand",
        "candidate_tier",
        "logFC",
        "FDR",
        "HC_beta_both_receptors",
        "HC_beta_strong_both_receptors"
    ),
    names(nearby)
)

nearby_compare <- merge(
    nearby[
        ,
        nearby_primary_cols,
        drop = FALSE
    ],
    ligand_compare[
        ,
        c(
            "ligand",
            "sens1600_beta_both_receptors",
            "reproduced_beta_both_receptors",
            "reproduced_strong_beta_receptors",
            "sens1600_beta_supported",
            "beta_supported_in_both_representations",
            "strong_beta_supported_in_both_representations"
        ),
        drop = FALSE
    ],
    by = "ligand",
    all.x = TRUE,
    sort = FALSE
)

nearby_compare <- nearby_compare[
    match(
        nearby$ligand,
        nearby_compare$ligand
    ),
    ,
    drop = FALSE
]

write.csv(
    nearby_compare,
    file.path(
        islet_dir,
        "Table_P3_NearbyNeuron_HC_1492_vs_1600_Sensitivity.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 15. Global sensitivity summary
# ============================================================

n_primary_pairs <- sum(
    pair_compare$primary1492_beta_both
)

n_sens_pairs <- sum(
    pair_compare$sens1600_beta_both
)

n_reproduced_pairs <- sum(
    pair_compare$beta_primary_reproduced
)

primary_ligands <- unique(
    pair_compare$ligand[
        pair_compare$primary1492_beta_both
    ]
)

sens_ligands <- unique(
    pair_compare$ligand[
        pair_compare$sens1600_beta_both
    ]
)

reproduced_ligands <- intersect(
    primary_ligands,
    sens_ligands
)

cat("\n============================================\n")
cat("1492 vs 1600 BETA-COMPATIBILITY SENSITIVITY\n")
cat("============================================\n")

cat(
    "HC LR pairs tested:",
    nrow(pair_compare),
    "\n"
)

cat(
    "Primary 1492 beta-both pairs:",
    n_primary_pairs,
    "\n"
)

cat(
    "Sensitivity 1600 beta-both pairs:",
    n_sens_pairs,
    "\n"
)

cat(
    "Primary beta-both pairs reproduced in 1600:",
    n_reproduced_pairs,
    "\n"
)

cat(
    "Fraction of primary beta-both pairs reproduced:",
    if (n_primary_pairs > 0) {
        round(
            n_reproduced_pairs /
                n_primary_pairs,
            4
        )
    } else {
        NA
    },
    "\n"
)

cat(
    "Primary 1492 beta-both ligands:",
    length(primary_ligands),
    "\n"
)

cat(
    "Sensitivity 1600 beta-both ligands:",
    length(sens_ligands),
    "\n"
)

cat(
    "Primary beta-both ligands retained in 1600:",
    length(reproduced_ligands),
    "\n"
)

cat(
    "Fraction of primary beta-both ligands retained:",
    if (length(primary_ligands) > 0) {
        round(
            length(reproduced_ligands) /
                length(primary_ligands),
            4
        )
    } else {
        NA
    },
    "\n"
)

cat(
    "Pair-level beta primary classification agreement:",
    round(
        mean(
            pair_compare$
                beta_primary_class_agreement
        ),
        4
    ),
    "\n"
)

cat(
    "Pair-level beta strong classification agreement:",
    round(
        mean(
            pair_compare$
                beta_strong_class_agreement
        ),
        4
    ),
    "\n"
)

# ============================================================
# 16. Print 13 nearby-neuron sensitivity results
# ============================================================

cat("\n============================================\n")
cat("13 NEARBY-NEURON CANDIDATES: SENSITIVITY\n")
cat("============================================\n")

display_cols <- c(
    "ligand",
    "candidate_tier",
    "HC_beta_both_receptors",
    "sens1600_beta_both_receptors",
    "reproduced_beta_both_receptors",
    "reproduced_strong_beta_receptors",
    "beta_supported_in_both_representations",
    "strong_beta_supported_in_both_representations"
)

print(
    nearby_compare[
        ,
        intersect(
            display_cols,
            names(nearby_compare)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)

# ============================================================
# 17. Freeze analysis method
# ============================================================

method <- data.frame(

    dataset =
        "GSE81608_current_GEO_1600",

    role =
        "dataset_representation_sensitivity_analysis",

    assay =
        "RPKM",

    cells =
        ncol(expr),

    genes =
        nrow(expr),

    min_cells_per_donor_celltype =
        MIN_CELLS_PER_DONOR,

    primary_within_donor_detection_fraction =
        PRIMARY_CELL_FRACTION,

    strict_within_donor_detection_fraction =
        STRICT_CELL_FRACTION,

    min_evaluable_donors =
        MIN_EVALUABLE_DONORS,

    primary_required_donor_fraction =
        PRIMARY_DONOR_FRACTION,

    strong_required_donor_fraction =
        STRONG_DONOR_FRACTION,

    interpretation =
        paste(
            "Sensitivity analysis of receptor detection",
            "using the current 1600-cell GEO",
            "representation of GSE81608;",
            "this is not an independent cohort",
            "or external validation."
        ),

    stringsAsFactors = FALSE
)

write.csv(
    method,
    file.path(
        islet_dir,
        "GSE81608_1600_ReceptorCompatibility_Method.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 18. Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "04f_GSE81608_1600_receptor_sensitivity_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nMethod:\n")
print(method)

cat("\nEvaluable donors:\n")
print(eval_table)

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nHC receptor overlap:\n")

print(
    c(
        candidate =
            length(candidate_receptors),

        present =
            length(receptors_present),

        absent =
            length(receptors_absent)
    )
)

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04f COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

