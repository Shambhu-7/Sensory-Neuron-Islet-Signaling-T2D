
# ============================================================
# PAPER 3
# Script 02b
# Leave-one-donor-out robustness analysis
#
# Tests whether gene-level and pathway-level conclusions
# depend strongly on any single donor.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

project_dir <- normalizePath(".", mustWork = TRUE)

processed_dir <- file.path(project_dir, "data", "processed")
redox_dir <- file.path(project_dir, "results", "redox")
repro_dir <- file.path(project_dir, "reproducibility")

# ------------------------------------------------------------
# Load existing results
# ------------------------------------------------------------

pb_counts <- readRDS(
    file.path(
        processed_dir,
        "GSE295206_donor_class_pseudobulk_counts.rds"
    )
)

sample_info <- read.csv(
    file.path(
        redox_dir,
        "GSE295206_donor_pseudobulk_sample_info.csv"
    ),
    stringsAsFactors = FALSE
)

full_de <- read.csv(
    file.path(
        redox_dir,
        "Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv"
    ),
    stringsAsFactors = FALSE
)

membership <- read.csv(
    file.path(
        redox_dir,
        "GSE295206_focused_MSigDB_membership_used.csv"
    ),
    stringsAsFactors = FALSE
)

full_camera_fixed <- read.csv(
    file.path(
        redox_dir,
        "Table_P3_CAMERA_Focused_Redox_fixedCor001.csv"
    ),
    stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

stopifnot(
    all(sample_info$group_key %in% colnames(pb_counts))
)

pb_counts <- pb_counts[
    ,
    sample_info$group_key,
    drop = FALSE
]

donors <- unique(sample_info$donor_id)

cat("\nDonors:\n")
print(donors)

cat(
    "\nRunning",
    length(donors),
    "leave-one-donor-out models.\n"
)

# ------------------------------------------------------------
# Focused gene sets
# ------------------------------------------------------------

membership <- unique(
    membership[
        ,
        c(
            "gs_name",
            "gene_symbol"
        )
    ]
)

focused_sets <- split(
    membership$gene_symbol,
    membership$gs_name
)

# ------------------------------------------------------------
# Function: run one edgeR + CAMERA model
# ------------------------------------------------------------

run_lodo <- function(excluded_donor) {

    cat("\n============================================\n")
    cat("EXCLUDING DONOR:", excluded_donor, "\n")
    cat("============================================\n")

    si <- sample_info[
        sample_info$donor_id != excluded_donor,
        ,
        drop = FALSE
    ]

    counts <- pb_counts[
        ,
        si$group_key,
        drop = FALSE
    ]

    si$donor_id <- factor(si$donor_id)

    si$annotation_class <- factor(
        si$annotation_class,
        levels = c(
            "Nearby_neuron",
            "Nageotte"
        )
    )

    design <- model.matrix(
        ~ donor_id + annotation_class,
        data = si
    )

    rownames(design) <- si$group_key

    coef_name <- "annotation_classNageotte"

    coef_index <- which(
        colnames(design) == coef_name
    )

    if (length(coef_index) != 1) {
        stop("Contrast coefficient not found.")
    }

    # --------------------------------------------------------
    # edgeR
    # --------------------------------------------------------

    y <- DGEList(
        counts = counts
    )

    keep <- filterByExpr(
        y,
        design = design
    )

    y <- y[
        keep,
        ,
        keep.lib.sizes = FALSE
    ]

    y <- calcNormFactors(
        y,
        method = "TMM"
    )

    y <- estimateDisp(
        y,
        design,
        robust = TRUE
    )

    fit <- glmQLFit(
        y,
        design,
        robust = TRUE
    )

    test <- glmQLFTest(
        fit,
        coef = coef_index
    )

    de <- topTags(
        test,
        n = Inf,
        sort.by = "PValue"
    )$table

    de$gene <- rownames(de)

    de$excluded_donor <- excluded_donor

    # --------------------------------------------------------
    # CAMERA
    # --------------------------------------------------------

    v <- voom(
        y,
        design,
        plot = FALSE
    )

    pathway_index <- lapply(
        focused_sets,
        function(gs) {

            which(
                rownames(v) %in% gs
            )
        }
    )

    pathway_index <- pathway_index[
        lengths(pathway_index) >= 5
    ]

    cam_fixed <- camera(
        v,
        index = pathway_index,
        design = design,
        contrast = coef_index,
        inter.gene.cor = 0.01
    )

    cam_fixed$gene_set <- rownames(cam_fixed)
    rownames(cam_fixed) <- NULL

    cam_fixed$excluded_donor <- excluded_donor
    cam_fixed$correlation_model <- "fixed_0.01"

    cam_est <- camera(
        v,
        index = pathway_index,
        design = design,
        contrast = coef_index,
        inter.gene.cor = NA
    )

    cam_est$gene_set <- rownames(cam_est)
    rownames(cam_est) <- NULL

    cam_est$excluded_donor <- excluded_donor
    cam_est$correlation_model <- "estimated"

    list(
        de = de,
        camera_fixed = cam_fixed,
        camera_estimated = cam_est
    )
}

# ------------------------------------------------------------
# Run all six LODO analyses
# ------------------------------------------------------------

lodo_results <- lapply(
    donors,
    run_lodo
)

names(lodo_results) <- donors

# ------------------------------------------------------------
# Combine gene-level results
# ------------------------------------------------------------

de_all <- do.call(
    rbind,
    lapply(
        lodo_results,
        `[[`,
        "de"
    )
)

rownames(de_all) <- NULL

write.csv(
    de_all,
    file.path(
        redox_dir,
        "Table_P3_LODO_All_Gene_Results.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Gene-level robustness relative to full model
# ------------------------------------------------------------

gene_robustness <- data.frame(
    gene = full_de$gene,
    full_logFC = full_de$logFC,
    full_FDR = full_de$FDR,
    stringsAsFactors = FALSE
)

for (d in donors) {

    tmp <- lodo_results[[d]]$de

    m <- match(
        gene_robustness$gene,
        tmp$gene
    )

    gene_robustness[[paste0("logFC_without_", d)]] <- tmp$logFC[m]

    gene_robustness[[paste0("FDR_without_", d)]] <- tmp$FDR[m]
}

logfc_cols <- grep(
    "^logFC_without_",
    names(gene_robustness),
    value = TRUE
)

fdr_cols <- grep(
    "^FDR_without_",
    names(gene_robustness),
    value = TRUE
)

gene_robustness$n_LODO_tested <- rowSums(
    !is.na(
        gene_robustness[
            ,
            logfc_cols,
            drop = FALSE
        ]
    )
)

gene_robustness$n_LODO_same_direction <- apply(
    gene_robustness,
    1,
    function(z) {

        full_effect <- as.numeric(
            z["full_logFC"]
        )

        vals <- as.numeric(
            z[logfc_cols]
        )

        sum(
            sign(vals) == sign(full_effect),
            na.rm = TRUE
        )
    }
)

gene_robustness$n_LODO_FDR05 <- apply(
    gene_robustness[
        ,
        fdr_cols,
        drop = FALSE
    ],
    1,
    function(x) {
        sum(
            as.numeric(x) < 0.05,
            na.rm = TRUE
        )
    }
)

gene_robustness$all_LODO_direction_concordant <-
    gene_robustness$n_LODO_same_direction ==
    gene_robustness$n_LODO_tested

gene_robustness$robust_FDR05_all_LODO <-
    gene_robustness$n_LODO_FDR05 ==
    length(donors)

write.csv(
    gene_robustness,
    file.path(
        redox_dir,
        "Table_P3_LODO_Gene_Robustness.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Per-donor descriptive effect sizes
#
# These are NOT inferential tests.
# They show whether each donor individually has the same
# direction of expression difference.
# ------------------------------------------------------------

y_all <- DGEList(
    counts = pb_counts
)

y_all <- calcNormFactors(
    y_all,
    method = "TMM"
)

logCPM <- cpm(
    y_all,
    log = TRUE,
    prior.count = 2
)

donor_effects <- matrix(
    NA_real_,
    nrow = nrow(logCPM),
    ncol = length(donors),
    dimnames = list(
        rownames(logCPM),
        donors
    )
)

for (d in donors) {

    nageotte_key <- paste(
        d,
        "Nageotte",
        sep = "__"
    )

    neuron_key <- paste(
        d,
        "Nearby_neuron",
        sep = "__"
    )

    donor_effects[, d] <-
        logCPM[, nageotte_key] -
        logCPM[, neuron_key]
}

donor_effects_df <- data.frame(
    gene = rownames(donor_effects),
    donor_effects,
    check.names = FALSE
)

write.csv(
    donor_effects_df,
    file.path(
        redox_dir,
        "Table_P3_PerDonor_LogCPM_Differences.csv"
    ),
    row.names = FALSE
)

# Add donor concordance to main robustness table

m <- match(
    gene_robustness$gene,
    rownames(donor_effects)
)

donor_matrix_for_de <- donor_effects[
    m,
    ,
    drop = FALSE
]

gene_robustness$n_donors_same_direction <- vapply(
    seq_len(nrow(gene_robustness)),
    function(i) {

        sum(
            sign(
                donor_matrix_for_de[i, ]
            ) ==
            sign(
                gene_robustness$full_logFC[i]
            ),
            na.rm = TRUE
        )
    },
    integer(1)
)

gene_robustness$all_6_donors_same_direction <-
    gene_robustness$n_donors_same_direction == 6

write.csv(
    gene_robustness,
    file.path(
        redox_dir,
        "Table_P3_LODO_Gene_Robustness.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# CAMERA LODO results
# ------------------------------------------------------------

camera_fixed_all <- do.call(
    rbind,
    lapply(
        lodo_results,
        `[[`,
        "camera_fixed"
    )
)

camera_est_all <- do.call(
    rbind,
    lapply(
        lodo_results,
        `[[`,
        "camera_estimated"
    )
)

rownames(camera_fixed_all) <- NULL
rownames(camera_est_all) <- NULL

write.csv(
    camera_fixed_all,
    file.path(
        redox_dir,
        "Table_P3_LODO_CAMERA_Fixed001_All.csv"
    ),
    row.names = FALSE
)

write.csv(
    camera_est_all,
    file.path(
        redox_dir,
        "Table_P3_LODO_CAMERA_Estimated_All.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Summarize pathway robustness
# ------------------------------------------------------------

pathways <- sort(
    unique(
        camera_fixed_all$gene_set
    )
)

pathway_summary <- lapply(
    pathways,
    function(gs) {

        fixed <- camera_fixed_all[
            camera_fixed_all$gene_set == gs,
            ,
            drop = FALSE
        ]

        estimated <- camera_est_all[
            camera_est_all$gene_set == gs,
            ,
            drop = FALSE
        ]

        full_row <- full_camera_fixed[
            full_camera_fixed$gene_set == gs,
            ,
            drop = FALSE
        ]

        full_direction <- if (
            nrow(full_row) == 1
        ) {
            full_row$Direction
        } else {
            NA_character_
        }

        data.frame(

            gene_set = gs,

            full_direction =
                full_direction,

            full_FDR_fixed001 =
                if (nrow(full_row) == 1)
                    full_row$FDR
                else NA_real_,

            fixed_same_direction =
                sum(
                    fixed$Direction ==
                        full_direction,
                    na.rm = TRUE
                ),

            fixed_FDR05 =
                sum(
                    fixed$FDR < 0.05,
                    na.rm = TRUE
                ),

            fixed_P05 =
                sum(
                    fixed$PValue < 0.05,
                    na.rm = TRUE
                ),

            estimated_same_direction =
                sum(
                    estimated$Direction ==
                        full_direction,
                    na.rm = TRUE
                ),

            estimated_FDR05 =
                sum(
                    estimated$FDR < 0.05,
                    na.rm = TRUE
                ),

            estimated_P05 =
                sum(
                    estimated$PValue < 0.05,
                    na.rm = TRUE
                ),

            stringsAsFactors = FALSE
        )
    }
)

pathway_summary <- do.call(
    rbind,
    pathway_summary
)

write.csv(
    pathway_summary,
    file.path(
        redox_dir,
        "Table_P3_LODO_CAMERA_Pathway_Robustness.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Specifically inspect anchor genes
# ------------------------------------------------------------

anchor <- read.csv(
    file.path(
        redox_dir,
        "Table_P3_Redox_Stress_Anchor_Genes.csv"
    ),
    stringsAsFactors = FALSE
)

anchor_robustness <- merge(
    anchor[
        ,
        c(
            "gene",
            "logFC",
            "FDR",
            "direction"
        )
    ],
    gene_robustness[
        ,
        c(
            "gene",
            "n_LODO_tested",
            "n_LODO_same_direction",
            "n_LODO_FDR05",
            "all_LODO_direction_concordant",
            "n_donors_same_direction",
            "all_6_donors_same_direction"
        )
    ],
    by = "gene",
    all.x = TRUE
)

write.csv(
    anchor_robustness,
    file.path(
        redox_dir,
        "Table_P3_Redox_Anchor_LODO_Robustness.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LODO PATHWAY ROBUSTNESS\n")
cat("============================================\n")

print(
    pathway_summary,
    row.names = FALSE
)

cat("\n============================================\n")
cat("REDOX ANCHOR ROBUSTNESS\n")
cat("============================================\n")

print(
    anchor_robustness[
        anchor_robustness$FDR < 0.05,
        ,
        drop = FALSE
    ],
    row.names = FALSE
)

cat("\n============================================\n")
cat("GLOBAL GENE ROBUSTNESS\n")
cat("============================================\n")

cat(
    "Full-model FDR < 0.05 genes:",
    sum(
        gene_robustness$full_FDR < 0.05
    ),
    "\n"
)

cat(
    "FDR < 0.05 genes with same direction in all 6 LODO models:",
    sum(
        gene_robustness$full_FDR < 0.05 &
        gene_robustness$all_LODO_direction_concordant,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "FDR < 0.05 genes with same direction in all 6 individual donors:",
    sum(
        gene_robustness$full_FDR < 0.05 &
        gene_robustness$all_6_donors_same_direction,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "FDR < 0.05 in full model AND every LODO model:",
    sum(
        gene_robustness$full_FDR < 0.05 &
        gene_robustness$robust_FDR05_all_LODO,
        na.rm = TRUE
    ),
    "\n"
)

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

sink(
    file.path(
        repro_dir,
        "02b_LODO_sessionInfo.txt"
    )
)

sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 02b COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

