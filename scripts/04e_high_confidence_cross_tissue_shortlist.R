
options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", mustWork = TRUE)

islet_dir <- file.path(project_dir, "results", "islet")
repro_dir <- file.path(project_dir, "reproducibility")

dir.create(islet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 04e: HIGH-CONFIDENCE CROSS-TISSUE SHORTLIST\n")
cat("============================================\n")

# ============================================================
# 1. Load 04d pair-level compatibility table
# ============================================================

input_file <- file.path(
    islet_dir,
    "Table_P3_DRG_Islet_LR_Compatibility_1492cells.csv"
)

if (!file.exists(input_file)) {
    stop("04d compatibility table not found.")
}

x <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat(
    "Input LR pairs:",
    nrow(x),
    "\n"
)

# ============================================================
# 2. Validate required columns
# ============================================================

required <- c(
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

missing <- setdiff(
    required,
    names(x)
)

if (length(missing)) {
    stop(
        "Missing columns:\n",
        paste(missing, collapse = "\n")
    )
}

to_bool <- function(z) {

    if (is.logical(z)) {
        return(z)
    }

    tolower(
        trimws(
            as.character(z)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
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
    x[[cc]] <- to_bool(x[[cc]])
}

# ============================================================
# 3. High-confidence LR edge subset
# ============================================================

hc <- x[
    x$high_confidence_sensitivity,
    ,
    drop = FALSE
]

cat(
    "High-confidence LR pairs:",
    nrow(hc),
    "\n"
)

cat(
    "High-confidence ligands:",
    length(unique(hc$ligand)),
    "\n"
)

cat(
    "High-confidence receptors:",
    length(unique(hc$receptor)),
    "\n"
)

# ============================================================
# 4. Compatibility definitions
# ============================================================

hc$beta_primary_any <-
    hc$beta_ND_primary_compatible |
    hc$beta_T2D_primary_compatible

hc$beta_primary_both <-
    hc$beta_ND_primary_compatible &
    hc$beta_T2D_primary_compatible

hc$beta_strong_any <-
    hc$beta_ND_strong_compatible |
    hc$beta_T2D_strong_compatible

hc$beta_strong_both <-
    hc$beta_ND_strong_compatible &
    hc$beta_T2D_strong_compatible

hc$alpha_primary_any <-
    hc$alpha_ND_primary_compatible |
    hc$alpha_T2D_primary_compatible

hc$alpha_primary_both <-
    hc$alpha_ND_primary_compatible &
    hc$alpha_T2D_primary_compatible

# ============================================================
# 5. Edge-level interpretation
#
# Highest-confidence cross-tissue edge:
# - passes frozen OmniPath HC filter
# - receptor supported in BOTH ND and T2D beta cells
#
# Very stringent:
# - above +
# - strong receptor support in BOTH conditions
# ============================================================

hc$beta_cross_tissue_supported <-
    hc$beta_primary_both

hc$beta_cross_tissue_strong <-
    hc$beta_strong_both

hc$alpha_cross_tissue_supported <-
    hc$alpha_primary_both

# ============================================================
# 6. Source plausibility
#
# Nearby-neuron enrichment is treated as greater sensory-neuron
# source plausibility than Nageotte-associated enrichment.
#
# This is still NOT proof of neuronal origin.
# ============================================================

hc$source_class <- ifelse(
    hc$enriched_region ==
        "Nearby_neuron_region",
    "Nearby_neuron_enriched",
    "Nageotte_associated"
)

hc$sensory_source_priority <-
    hc$source_class ==
        "Nearby_neuron_enriched" &
    hc$beta_cross_tissue_supported

hc$sensory_source_strict <-
    hc$source_class ==
        "Nearby_neuron_enriched" &
    hc$beta_cross_tissue_strong

# ============================================================
# 7. Collapse to ligand level
# ============================================================

collapse_unique <- function(z) {

    z <- sort(
        unique(
            as.character(z[
                !is.na(z) &
                nzchar(as.character(z))
            ])
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

ligand_groups <- split(
    seq_len(nrow(hc)),
    hc$ligand
)

ligand_summary <- do.call(
    rbind,
    lapply(
        names(ligand_groups),
        function(g) {

            z <- hc[
                ligand_groups[[g]],
                ,
                drop = FALSE
            ]

            beta_both <- z$receptor[
                z$beta_cross_tissue_supported
            ]

            beta_strong <- z$receptor[
                z$beta_cross_tissue_strong
            ]

            beta_any <- z$receptor[
                z$beta_primary_any
            ]

            alpha_both <- z$receptor[
                z$alpha_cross_tissue_supported
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

                n_HC_LR_edges =
                    nrow(z),

                n_HC_beta_any_receptors =
                    length(
                        unique(beta_any)
                    ),

                n_HC_beta_both_receptors =
                    length(
                        unique(beta_both)
                    ),

                n_HC_beta_strong_both_receptors =
                    length(
                        unique(beta_strong)
                    ),

                HC_beta_any_receptors =
                    collapse_unique(beta_any),

                HC_beta_both_receptors =
                    collapse_unique(beta_both),

                HC_beta_strong_both_receptors =
                    collapse_unique(beta_strong),

                n_HC_alpha_both_receptors =
                    length(
                        unique(alpha_both)
                    ),

                HC_alpha_both_receptors =
                    collapse_unique(alpha_both),

                beta_HC_supported =
                    length(beta_both) > 0,

                beta_HC_strong =
                    length(beta_strong) > 0,

                alpha_HC_supported =
                    length(alpha_both) > 0,

                sensory_source_priority =
                    z$enriched_region[1] ==
                        "Nearby_neuron_region" &&
                    length(beta_both) > 0,

                sensory_source_strict =
                    z$enriched_region[1] ==
                        "Nearby_neuron_region" &&
                    length(beta_strong) > 0,

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(ligand_summary) <- NULL

# ============================================================
# 8. Rank without inventing a statistical score
# ============================================================

tier_rank <- match(
    ligand_summary$candidate_tier,
    c(
        "Tier_A_strict",
        "Tier_B_LODO_robust",
        "Tier_C_directionally_robust"
    )
)

tier_rank[is.na(tier_rank)] <- 99

ligand_summary <- ligand_summary[
    order(
        !ligand_summary$sensory_source_strict,
        !ligand_summary$sensory_source_priority,
        !ligand_summary$beta_HC_strong,
        !ligand_summary$beta_HC_supported,
        tier_rank,
        ligand_summary$FDR
    ),
    ,
    drop = FALSE
]

# ============================================================
# 9. Main cross-tissue tables
# ============================================================

write.csv(
    hc,
    file.path(
        islet_dir,
        "Table_P3_HighConfidence_LR_Islet_Compatibility.csv"
    ),
    row.names = FALSE
)

write.csv(
    ligand_summary,
    file.path(
        islet_dir,
        "Table_P3_HighConfidence_Ligand_Islet_Summary.csv"
    ),
    row.names = FALSE
)

sensory_priority <- ligand_summary[
    ligand_summary$sensory_source_priority,
    ,
    drop = FALSE
]

write.csv(
    sensory_priority,
    file.path(
        islet_dir,
        "Table_P3_SensorySource_BetaCell_HC_Shortlist.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 10. Tier-A sensory-source subset
# ============================================================

tierA_sensory <- ligand_summary[
    ligand_summary$candidate_tier ==
        "Tier_A_strict" &
    ligand_summary$sensory_source_priority,
    ,
    drop = FALSE
]

write.csv(
    tierA_sensory,
    file.path(
        islet_dir,
        "Table_P3_TierA_SensorySource_BetaCell_HC_Shortlist.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 11. Nageotte-associated compartment kept separately
# ============================================================

nageotte_beta <- ligand_summary[
    ligand_summary$enriched_region ==
        "Nageotte_associated_region" &
    ligand_summary$beta_HC_supported,
    ,
    drop = FALSE
]

write.csv(
    nageotte_beta,
    file.path(
        islet_dir,
        "Table_P3_NageotteAssociated_BetaCell_HC_Ligands.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 12. Specific neuronal candidates
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

neuronal_edges <- x[
    x$ligand %in%
        neuronal_candidates,
    ,
    drop = FALSE
]

neuronal_edges <- neuronal_edges[
    order(
        match(
            neuronal_edges$ligand,
            neuronal_candidates
        ),
        !neuronal_edges$high_confidence_sensitivity,
        neuronal_edges$receptor
    ),
    ,
    drop = FALSE
]

neuronal_edge_cols <- intersect(
    c(
        "ligand",
        "receptor",
        "candidate_tier",
        "logFC",
        "FDR",
        "high_confidence_sensitivity",
        "n_query_resources",
        "max_curation_effort",
        "ligand_consensus_score",
        "receptor_consensus_score",
        "beta_ND_primary_compatible",
        "beta_T2D_primary_compatible",
        "beta_ND_strong_compatible",
        "beta_T2D_strong_compatible",
        "alpha_ND_primary_compatible",
        "alpha_T2D_primary_compatible"
    ),
    names(neuronal_edges)
)

neuronal_edges_out <- neuronal_edges[
    ,
    neuronal_edge_cols,
    drop = FALSE
]

write.csv(
    neuronal_edges_out,
    file.path(
        islet_dir,
        "Table_P3_NeuronalCandidates_HC_EdgeAudit.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 13. Console summary
# ============================================================

cat("\n============================================\n")
cat("HIGH-CONFIDENCE NETWORK SUMMARY\n")
cat("============================================\n")

cat(
    "High-confidence LR pairs:",
    nrow(hc),
    "\n"
)

cat(
    "High-confidence ligands:",
    length(unique(hc$ligand)),
    "\n"
)

cat(
    "HC LR pairs beta-compatible BOTH conditions:",
    sum(hc$beta_cross_tissue_supported),
    "\n"
)

cat(
    "HC LR pairs STRONGLY beta-compatible BOTH conditions:",
    sum(hc$beta_cross_tissue_strong),
    "\n"
)

cat(
    "HC ligands beta-compatible BOTH conditions:",
    sum(ligand_summary$beta_HC_supported),
    "\n"
)

cat(
    "HC ligands STRONGLY beta-compatible BOTH conditions:",
    sum(ligand_summary$beta_HC_strong),
    "\n"
)

cat(
    "Nearby-neuron ligands with HC beta support:",
    nrow(sensory_priority),
    "\n"
)

cat(
    "Tier-A nearby-neuron ligands with HC beta support:",
    nrow(tierA_sensory),
    "\n"
)

cat(
    "Nageotte-associated ligands with HC beta support:",
    nrow(nageotte_beta),
    "\n"
)

# ============================================================
# 14. Tier-A sensory-neuron-priority shortlist
# ============================================================

cat("\n============================================\n")
cat("TIER-A SENSORY-SOURCE HIGH-CONFIDENCE SHORTLIST\n")
cat("============================================\n")

display_cols <- c(
    "ligand",
    "logFC",
    "FDR",
    "n_HC_beta_both_receptors",
    "n_HC_beta_strong_both_receptors",
    "HC_beta_both_receptors",
    "HC_beta_strong_both_receptors",
    "n_HC_alpha_both_receptors",
    "HC_alpha_both_receptors"
)

print(
    tierA_sensory[
        ,
        intersect(
            display_cols,
            names(tierA_sensory)
        ),
        drop = FALSE
    ],
    row.names = FALSE
)

# ============================================================
# 15. Exact neuronal edge audit
# ============================================================

cat("\n============================================\n")
cat("NEURONAL CANDIDATE HIGH-CONFIDENCE EDGE AUDIT\n")
cat("============================================\n")

print(
    neuronal_edges_out,
    row.names = FALSE
)

# ============================================================
# 16. Nageotte-associated top candidates
# ============================================================

cat("\n============================================\n")
cat("NAGEOTTE-ASSOCIATED HC BETA-COMPATIBLE LIGANDS\n")
cat("============================================\n")

print(
    head(
        nageotte_beta[
            ,
            intersect(
                c(
                    "ligand",
                    "candidate_tier",
                    "logFC",
                    "FDR",
                    "n_HC_beta_both_receptors",
                    "HC_beta_both_receptors"
                ),
                names(nageotte_beta)
            ),
            drop = FALSE
        ],
        30
    ),
    row.names = FALSE
)

# ============================================================
# 17. Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "04e_high_confidence_cross_tissue_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nInput file:\n")
print(input_file)

cat("\nHigh-confidence LR pairs:\n")
print(nrow(hc))

cat("\nTier-A sensory-source shortlist:\n")
print(tierA_sensory)

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 04e COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

