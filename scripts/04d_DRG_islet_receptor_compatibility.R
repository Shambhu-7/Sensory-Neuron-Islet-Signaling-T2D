
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(Biobase)
})

project_dir <- normalizePath(".", mustWork = TRUE)

raw_dir <- file.path(project_dir, "data", "raw", "GSE81608")
ligand_dir <- file.path(project_dir, "results", "ligands")
islet_dir <- file.path(project_dir, "results", "islet")
repro_dir <- file.path(project_dir, "reproducibility")

dir.create(islet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 04d: DRG -> ISLET RECEPTOR COMPATIBILITY\n")
cat("============================================\n")

# ============================================================
# 1. Required inputs
# ============================================================

xin_file <- file.path(
    raw_dir,
    "XinT2Deset_MuSiC_1492cells.rds"
)

ligand_file <- file.path(
    ligand_dir,
    "Table_P3_Prioritized_Robust_DRG_Ligands.csv"
)

pair_file <- file.path(
    ligand_dir,
    "Table_P3_Prioritized_Ligand_Receptor_Pairs.csv"
)

required <- c(
    xin_file,
    ligand_file,
    pair_file
)

missing <- required[!file.exists(required)]

if (length(missing)) {
    stop(
        "Missing required files:\n",
        paste(missing, collapse = "\n")
    )
}

# ============================================================
# 2. Load primary 1492-cell Xin dataset
# ============================================================

xin <- readRDS(xin_file)

if (!inherits(xin, "ExpressionSet")) {
    stop("Xin object is not an ExpressionSet.")
}

expr <- Biobase::exprs(xin)
meta <- Biobase::pData(xin)

cat(
    "Xin expression matrix:",
    nrow(expr),
    "genes x",
    ncol(expr),
    "cells\n"
)

if (
    nrow(expr) != 39849L ||
    ncol(expr) != 1492L
) {
    stop(
        "Xin object dimensions are not 39849 x 1492."
    )
}

needed_meta <- c(
    "SubjectName",
    "cellType",
    "Disease"
)

if (!all(needed_meta %in% names(meta))) {
    stop(
        "Missing metadata fields: ",
        paste(
            setdiff(needed_meta, names(meta)),
            collapse = ", "
        )
    )
}

# ============================================================
# 3. Standardize metadata
# ============================================================

meta$donor_id <- as.character(
    meta$SubjectName
)

meta$cell_type <- tolower(
    trimws(
        as.character(meta$cellType)
    )
)

meta$condition_raw <- as.character(
    meta$Disease
)

disease_lower <- tolower(
    trimws(meta$condition_raw)
)

meta$condition <- ifelse(
    grepl(
        "t2d|type.?2|diab",
        disease_lower
    ),
    "T2D",
    "non-diabetic"
)

if (!identical(
    rownames(meta),
    colnames(expr)
)) {

    m <- match(
        colnames(expr),
        rownames(meta)
    )

    if (anyNA(m)) {
        stop(
            "Could not align phenotype metadata to expression matrix."
        )
    }

    meta <- meta[m, , drop = FALSE]
}

stopifnot(
    identical(
        rownames(meta),
        colnames(expr)
    )
)

cat("\nSubjects:", length(unique(meta$donor_id)), "\n")

cat("\nDisease groups:\n")
print(table(meta$condition))

cat("\nCell types:\n")
print(table(meta$cell_type))

cat("\nDonor x disease:\n")
print(
    table(
        meta$donor_id,
        meta$condition
    )
)

# ============================================================
# 4. Load robust DRG ligands and frozen LR pairs
# ============================================================

ligands <- read.csv(
    ligand_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

pairs <- read.csv(
    pair_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

needed_pair <- c(
    "ligand",
    "receptor",
    "candidate_tier",
    "logFC",
    "FDR",
    "enriched_region"
)

if (!all(needed_pair %in% names(pairs))) {
    stop(
        "LR pair table lacks required columns."
    )
}

pairs <- pairs[
    !duplicated(
        paste(
            pairs$ligand,
            pairs$receptor,
            sep = "|||"
        )
    ),
    ,
    drop = FALSE
]

cat(
    "\nPrioritized DRG ligands:",
    length(unique(pairs$ligand)),
    "\n"
)

cat(
    "Prioritized ligand-receptor pairs:",
    nrow(pairs),
    "\n"
)

cat(
    "Unique candidate receptors:",
    length(unique(pairs$receptor)),
    "\n"
)

# ============================================================
# 5. Validate gene identifiers
# ============================================================

candidate_receptors <- sort(
    unique(
        pairs$receptor
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
    "Candidate receptors:",
    length(candidate_receptors),
    "\n"
)

cat(
    "Present in Xin expression matrix:",
    length(receptors_present),
    "\n"
)

cat(
    "Absent from Xin expression matrix:",
    length(receptors_absent),
    "\n"
)

if (length(receptors_present) < 50) {
    stop(
        paste(
            "Unexpectedly low receptor-symbol overlap.",
            "Check gene identifiers before proceeding."
        )
    )
}

write.csv(
    data.frame(
        receptor = receptors_absent
    ),
    file.path(
        islet_dir,
        "GSE81608_1492_LR_receptors_absent_from_matrix.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 6. Restrict to alpha and beta cells
# ============================================================

target_celltypes <- c(
    "alpha",
    "beta"
)

keep_cells <- meta$cell_type %in%
    target_celltypes

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
cat("ALPHA/BETA DATASET\n")
cat("============================================\n")

cat(
    "Cells:",
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
#
# Detection = count > 0.
#
# Donor evaluability:
# >=5 cells for that donor/cell-type combination.
#
# Primary donor support:
# receptor detected in >=10% cells.
#
# Strict donor support:
# receptor detected in >=20% cells.
# ============================================================

MIN_CELLS_PER_DONOR <- 5L
PRIMARY_CELL_FRACTION <- 0.10
STRICT_CELL_FRACTION <- 0.20

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

cat("\nCalculating donor-level receptor support...\n")

for (receptor in receptors_present) {

    values <- expr_ab[
        receptor,
        ,
        drop = TRUE
    ]

    for (g in names(groups)) {

        idx <- groups[[g]]

        x <- values[idx]

        n_cells <- length(x)

        n_detected <- sum(
            x > 0,
            na.rm = TRUE
        )

        fraction_detected <-
            n_detected / n_cells

        donor_results[[counter]] <-
            data.frame(

                receptor = receptor,

                donor_id =
                    meta_ab$donor_id[idx[1]],

                condition =
                    meta_ab$condition[idx[1]],

                cell_type =
                    meta_ab$cell_type[idx[1]],

                n_cells =
                    n_cells,

                n_detected =
                    n_detected,

                fraction_detected =
                    fraction_detected,

                mean_count =
                    mean(
                        x,
                        na.rm = TRUE
                    ),

                median_count =
                    median(
                        x,
                        na.rm = TRUE
                    ),

                evaluable =
                    n_cells >=
                    MIN_CELLS_PER_DONOR,

                support_any =
                    n_detected >= 1,

                support_10pct =
                    n_cells >=
                    MIN_CELLS_PER_DONOR &&
                    fraction_detected >=
                    PRIMARY_CELL_FRACTION,

                support_20pct =
                    n_cells >=
                    MIN_CELLS_PER_DONOR &&
                    fraction_detected >=
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
        paste0(
            "GSE81608_1492_Receptor_",
            "DonorLevel_Detection.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 8. Summarize across biological donors
#
# Primary compatibility:
# >=3 evaluable donors AND
# >=50% evaluable donors support receptor at >=10%.
#
# Strong compatibility:
# >=3 evaluable donors AND
# >=75% evaluable donors support receptor at >=20%.
# ============================================================

MIN_EVALUABLE_DONORS <- 3L
PRIMARY_DONOR_FRACTION <- 0.50
STRONG_DONOR_FRACTION <- 0.75

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

            eval <- z$evaluable

            n_eval <- sum(eval)

            n_any <- sum(
                z$support_any & eval
            )

            n_10 <- sum(
                z$support_10pct & eval
            )

            n_20 <- sum(
                z$support_20pct & eval
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

            median_detection <- if (
                n_eval > 0
            ) {
                median(
                    z$fraction_detected[eval],
                    na.rm = TRUE
                )
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

                n_total_donors =
                    nrow(z),

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
                    median_detection,

                primary_compatible =
                    n_eval >=
                    MIN_EVALUABLE_DONORS &&
                    frac_10 >=
                    PRIMARY_DONOR_FRACTION,

                strong_compatible =
                    n_eval >=
                    MIN_EVALUABLE_DONORS &&
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
        paste0(
            "GSE81608_1492_Receptor_",
            "DonorSupport_Summary.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 9. Evaluable donor summary
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

cat("\n============================================\n")
cat("EVALUABLE DONORS\n")
cat("============================================\n")

eval_table <- aggregate(
    evaluable ~ condition + cell_type,
    data = eval_summary,
    FUN = sum
)

print(
    eval_table,
    row.names = FALSE
)

# ============================================================
# 10. Helper to attach receptor support to LR pairs
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

alpha_ND <- get_support(
    "non-diabetic",
    "alpha",
    "alpha_ND"
)

alpha_T2D <- get_support(
    "T2D",
    "alpha",
    "alpha_T2D"
)

beta_ND <- get_support(
    "non-diabetic",
    "beta",
    "beta_ND"
)

beta_T2D <- get_support(
    "T2D",
    "beta",
    "beta_T2D"
)

pair_support <- pairs

pair_support <- merge(
    pair_support,
    alpha_ND,
    by = "receptor",
    all.x = TRUE,
    sort = FALSE
)

pair_support <- merge(
    pair_support,
    alpha_T2D,
    by = "receptor",
    all.x = TRUE,
    sort = FALSE
)

pair_support <- merge(
    pair_support,
    beta_ND,
    by = "receptor",
    all.x = TRUE,
    sort = FALSE
)

pair_support <- merge(
    pair_support,
    beta_T2D,
    by = "receptor",
    all.x = TRUE,
    sort = FALSE
)

# Receptors absent from the expression matrix remain FALSE,
# not NA, for compatibility classification.

compat_cols <- grep(
    "primary_compatible$|strong_compatible$",
    names(pair_support),
    value = TRUE
)

for (cc in compat_cols) {

    pair_support[[cc]][
        is.na(pair_support[[cc]])
    ] <- FALSE
}

# ============================================================
# 11. Cross-tissue compatibility classes
#
# These describe receptor presence only.
# They do NOT demonstrate signaling.
# They do NOT constitute T2D differential-expression tests.
# ============================================================

pair_support$beta_primary_any <-
    pair_support$beta_ND_primary_compatible |
    pair_support$beta_T2D_primary_compatible

pair_support$beta_primary_both <-
    pair_support$beta_ND_primary_compatible &
    pair_support$beta_T2D_primary_compatible

pair_support$beta_strong_any <-
    pair_support$beta_ND_strong_compatible |
    pair_support$beta_T2D_strong_compatible

pair_support$beta_strong_both <-
    pair_support$beta_ND_strong_compatible &
    pair_support$beta_T2D_strong_compatible

pair_support$alpha_primary_any <-
    pair_support$alpha_ND_primary_compatible |
    pair_support$alpha_T2D_primary_compatible

pair_support$alpha_primary_both <-
    pair_support$alpha_ND_primary_compatible &
    pair_support$alpha_T2D_primary_compatible

pair_support$alpha_strong_any <-
    pair_support$alpha_ND_strong_compatible |
    pair_support$alpha_T2D_strong_compatible

pair_support$alpha_strong_both <-
    pair_support$alpha_ND_strong_compatible &
    pair_support$alpha_T2D_strong_compatible

pair_support$islet_primary_any <-
    pair_support$beta_primary_any |
    pair_support$alpha_primary_any

pair_support$islet_strong_any <-
    pair_support$beta_strong_any |
    pair_support$alpha_strong_any

pair_support$beta_support_pattern <-
    ifelse(
        pair_support$beta_primary_both,
        "supported_both_conditions",
        ifelse(
            pair_support$beta_T2D_primary_compatible &
            !pair_support$beta_ND_primary_compatible,
            "T2D_only_descriptive",
            ifelse(
                pair_support$beta_ND_primary_compatible &
                !pair_support$beta_T2D_primary_compatible,
                "non_diabetic_only_descriptive",
                "not_primary_supported"
            )
        )
    )

pair_support$alpha_support_pattern <-
    ifelse(
        pair_support$alpha_primary_both,
        "supported_both_conditions",
        ifelse(
            pair_support$alpha_T2D_primary_compatible &
            !pair_support$alpha_ND_primary_compatible,
            "T2D_only_descriptive",
            ifelse(
                pair_support$alpha_ND_primary_compatible &
                !pair_support$alpha_T2D_primary_compatible,
                "non_diabetic_only_descriptive",
                "not_primary_supported"
            )
        )
    )

# ============================================================
# 12. Order pair table
# ============================================================

tier_rank <- match(
    pair_support$candidate_tier,
    c(
        "Tier_A_strict",
        "Tier_B_LODO_robust",
        "Tier_C_directionally_robust"
    )
)

tier_rank[
    is.na(tier_rank)
] <- 99

pair_support <- pair_support[
    order(
        !pair_support$beta_primary_both,
        !pair_support$beta_strong_any,
        !pair_support$beta_primary_any,
        !pair_support$alpha_primary_both,
        tier_rank,
        pair_support$FDR
    ),
    ,
    drop = FALSE
]

write.csv(
    pair_support,
    file.path(
        islet_dir,
        paste0(
            "Table_P3_DRG_Islet_LR_",
            "Compatibility_1492cells.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 13. Ligand-level cross-tissue summary
# ============================================================

collapse_unique <- function(x) {

    x <- as.character(x)

    x <- x[
        !is.na(x) &
        nzchar(x)
    ]

    if (!length(x)) {
        return(NA_character_)
    }

    paste(
        sort(unique(x)),
        collapse = ";"
    )
}

ligand_groups <- split(
    seq_len(nrow(pair_support)),
    pair_support$ligand
)

ligand_summary <- do.call(
    rbind,
    lapply(
        names(ligand_groups),
        function(g) {

            ii <- ligand_groups[[g]]

            z <- pair_support[
                ii,
                ,
                drop = FALSE
            ]

            beta_both_receptors <-
                z$receptor[
                    z$beta_primary_both
                ]

            beta_any_receptors <-
                z$receptor[
                    z$beta_primary_any
                ]

            beta_strong_receptors <-
                z$receptor[
                    z$beta_strong_any
                ]

            alpha_both_receptors <-
                z$receptor[
                    z$alpha_primary_both
                ]

            alpha_any_receptors <-
                z$receptor[
                    z$alpha_primary_any
                ]

            data.frame(

                ligand = g,

                candidate_tier =
                    z$candidate_tier[1],

                logFC =
                    z$logFC[1],

                FDR =
                    z$FDR[1],

                enriched_region =
                    z$enriched_region[1],

                highlighted_original_GSE295206_study =
                    if (
                        "highlighted_original_GSE295206_study" %in%
                        names(z)
                    ) {
                        z$highlighted_original_GSE295206_study[1]
                    } else {
                        NA
                    },

                n_curated_LR_pairs =
                    nrow(z),

                n_beta_primary_receptors =
                    length(
                        unique(beta_any_receptors)
                    ),

                n_beta_both_conditions_receptors =
                    length(
                        unique(beta_both_receptors)
                    ),

                n_beta_strong_receptors =
                    length(
                        unique(beta_strong_receptors)
                    ),

                beta_primary_receptors =
                    collapse_unique(
                        beta_any_receptors
                    ),

                beta_both_conditions_receptors =
                    collapse_unique(
                        beta_both_receptors
                    ),

                beta_strong_receptors =
                    collapse_unique(
                        beta_strong_receptors
                    ),

                n_alpha_primary_receptors =
                    length(
                        unique(alpha_any_receptors)
                    ),

                n_alpha_both_conditions_receptors =
                    length(
                        unique(alpha_both_receptors)
                    ),

                alpha_primary_receptors =
                    collapse_unique(
                        alpha_any_receptors
                    ),

                alpha_both_conditions_receptors =
                    collapse_unique(
                        alpha_both_receptors
                    ),

                beta_compatible =
                    any(
                        z$beta_primary_any
                    ),

                beta_stable_across_conditions =
                    any(
                        z$beta_primary_both
                    ),

                alpha_compatible =
                    any(
                        z$alpha_primary_any
                    ),

                alpha_stable_across_conditions =
                    any(
                        z$alpha_primary_both
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(ligand_summary) <- NULL

ligand_summary$tier_rank <- match(
    ligand_summary$candidate_tier,
    c(
        "Tier_A_strict",
        "Tier_B_LODO_robust",
        "Tier_C_directionally_robust"
    )
)

ligand_summary$tier_rank[
    is.na(ligand_summary$tier_rank)
] <- 99

ligand_summary <- ligand_summary[
    order(
        !ligand_summary$beta_stable_across_conditions,
        -ligand_summary$n_beta_strong_receptors,
        -ligand_summary$n_beta_primary_receptors,
        !ligand_summary$alpha_stable_across_conditions,
        ligand_summary$tier_rank,
        ligand_summary$FDR
    ),
    ,
    drop = FALSE
]

ligand_summary$tier_rank <- NULL

write.csv(
    ligand_summary,
    file.path(
        islet_dir,
        paste0(
            "Table_P3_DRG_Ligand_",
            "Human_Islet_Compatibility_1492cells.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 14. Beta-focused shortlist
# ============================================================

beta_shortlist <- ligand_summary[
    ligand_summary$beta_compatible,
    ,
    drop = FALSE
]

write.csv(
    beta_shortlist,
    file.path(
        islet_dir,
        paste0(
            "Table_P3_BetaCell_Compatible_",
            "DRG_Ligands_1492cells.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 15. Console summary
# ============================================================

cat("\n============================================\n")
cat("RECEPTOR COMPATIBILITY SUMMARY\n")
cat("============================================\n")

cat(
    "Candidate receptors present:",
    length(receptors_present),
    "/",
    length(candidate_receptors),
    "\n"
)

cat(
    "LR pairs with primary beta compatibility:",
    sum(pair_support$beta_primary_any),
    "\n"
)

cat(
    "LR pairs with beta compatibility in BOTH conditions:",
    sum(pair_support$beta_primary_both),
    "\n"
)

cat(
    "LR pairs with strong beta compatibility:",
    sum(pair_support$beta_strong_any),
    "\n"
)

cat(
    "LR pairs with primary alpha compatibility:",
    sum(pair_support$alpha_primary_any),
    "\n"
)

cat(
    "LR pairs with alpha compatibility in BOTH conditions:",
    sum(pair_support$alpha_primary_both),
    "\n"
)

cat(
    "LR pairs with any alpha/beta compatibility:",
    sum(pair_support$islet_primary_any),
    "\n"
)

cat("\n============================================\n")
cat("LIGAND-LEVEL SUMMARY\n")
cat("============================================\n")

cat(
    "Prioritized DRG ligands tested:",
    nrow(ligand_summary),
    "\n"
)

cat(
    "Ligands with beta-compatible receptor(s):",
    sum(ligand_summary$beta_compatible),
    "\n"
)

cat(
    "Ligands with beta receptor(s) supported in BOTH conditions:",
    sum(
        ligand_summary$beta_stable_across_conditions
    ),
    "\n"
)

cat(
    "Ligands with alpha-compatible receptor(s):",
    sum(ligand_summary$alpha_compatible),
    "\n"
)

cat(
    "Ligands with alpha receptor(s) supported in BOTH conditions:",
    sum(
        ligand_summary$alpha_stable_across_conditions
    ),
    "\n"
)

cat("\n============================================\n")
cat("TOP BETA-COMPATIBLE DRG LIGANDS\n")
cat("============================================\n")

show_cols <- c(
    "ligand",
    "candidate_tier",
    "logFC",
    "FDR",
    "enriched_region",
    "n_beta_primary_receptors",
    "n_beta_both_conditions_receptors",
    "n_beta_strong_receptors",
    "beta_primary_receptors",
    "beta_both_conditions_receptors",
    "highlighted_original_GSE295206_study"
)

show_cols <- intersect(
    show_cols,
    names(beta_shortlist)
)

print(
    head(
        beta_shortlist[
            ,
            show_cols,
            drop = FALSE
        ],
        50
    ),
    row.names = FALSE
)

# ============================================================
# 16. Selected neuronal-side candidates
# ============================================================

neuronal_candidates <- c(
    "CALCB",
    "CALCA",
    "TAC1",
    "NRG1",
    "NMB",
    "UCN",
    "ADCYAP1",
    "GAL",
    "PENK"
)

neuronal_check <- ligand_summary[
    ligand_summary$ligand %in%
        neuronal_candidates,
    ,
    drop = FALSE
]

neuronal_check <- neuronal_check[
    match(
        intersect(
            neuronal_candidates,
            neuronal_check$ligand
        ),
        neuronal_check$ligand
    ),
    ,
    drop = FALSE
]

cat("\n============================================\n")
cat("NEURONAL-SIDE CANDIDATE CHECK\n")
cat("============================================\n")

print(
    neuronal_check[
        ,
        intersect(
            c(
                "ligand",
                "candidate_tier",
                "logFC",
                "FDR",
                "n_beta_primary_receptors",
                "beta_primary_receptors",
                "beta_both_conditions_receptors",
                "n_alpha_primary_receptors",
                "alpha_primary_receptors"
            ),
            names(neuronal_check)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)

write.csv(
    neuronal_check,
    file.path(
        islet_dir,
        paste0(
            "Table_P3_NeuronalCandidate_",
            "IsletCompatibility_1492cells.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# 17. Analysis metadata
# ============================================================

analysis_metadata <- data.frame(

    dataset =
        "GSE81608_Xin_MuSiC_1492",

    assay =
        "read_counts",

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
            "Receptor compatibility is descriptive",
            "and does not demonstrate signaling or",
            "disease-dependent regulation."
        ),

    stringsAsFactors = FALSE
)

write.csv(
    analysis_metadata,
    file.path(
        islet_dir,
        "GSE81608_1492_ReceptorCompatibility_Method.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 18. Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "04d_DRG_islet_receptor_compatibility_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nThresholds:\n")
print(analysis_metadata)

cat("\nEvaluable donors:\n")
print(eval_table)

cat("\nInput dimensions:\n")
print(dim(expr))

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04d COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

