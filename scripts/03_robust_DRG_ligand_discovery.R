
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(edgeR)
})

project_dir <- normalizePath(".", mustWork = TRUE)

processed_dir <- file.path(project_dir, "data", "processed")
redox_dir <- file.path(project_dir, "results", "redox")
ligand_dir <- file.path(project_dir, "results", "ligands")
repro_dir <- file.path(project_dir, "reproducibility")

dir.create(ligand_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("SCRIPT 03: ROBUST DRG LIGAND DISCOVERY\n")
cat("============================================\n")

# ============================================================
# Required files
# ============================================================

required_files <- c(
    file.path(
        processed_dir,
        "GSE295206_donor_class_pseudobulk_counts.rds"
    ),
    file.path(
        redox_dir,
        "GSE295206_donor_pseudobulk_sample_info.csv"
    ),
    file.path(
        redox_dir,
        "Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv"
    ),
    file.path(
        redox_dir,
        "Table_P3_LODO_Gene_Robustness.csv"
    ),
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_curated_interactions_combined.tsv"
    ),
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_ligand_annotations_RAW.tsv"
    ),
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_receptor_annotations_RAW.tsv"
    ),
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_resource_status.csv"
    ),
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_retrieval_metadata.csv"
    )
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {
    stop(
        "Missing required files:\n",
        paste(missing_files, collapse = "\n")
    )
}

# ============================================================
# Load validated upstream data
# ============================================================

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
    stringsAsFactors = FALSE,
    check.names = FALSE
)

de <- read.csv(
    file.path(
        redox_dir,
        "Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

robustness <- read.csv(
    file.path(
        redox_dir,
        "Table_P3_LODO_Gene_Robustness.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# ============================================================
# Load frozen OmniPath snapshots
# ============================================================

interactions <- read.delim(
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_curated_interactions_combined.tsv"
    ),
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

lig_ann <- read.delim(
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_ligand_annotations_RAW.tsv"
    ),
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rec_ann <- read.delim(
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_receptor_annotations_RAW.tsv"
    ),
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

resource_status <- read.csv(
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_resource_status.csv"
    ),
    stringsAsFactors = FALSE
)

retrieval_meta <- read.csv(
    file.path(
        ligand_dir,
        "OmniPath_DIRECT_retrieval_metadata.csv"
    ),
    stringsAsFactors = FALSE
)

if (!all(resource_status$status == "SUCCESS")) {
    stop(
        "Not all frozen OmniPath resources were retrieved successfully."
    )
}

# ============================================================
# Validate input columns
# ============================================================

needed_int <- c(
    "source_genesymbol",
    "target_genesymbol",
    "query_resource"
)

if (!all(needed_int %in% names(interactions))) {

    stop(
        "Frozen interaction table lacks: ",
        paste(
            setdiff(
                needed_int,
                names(interactions)
            ),
            collapse = ", "
        )
    )
}

needed_ann <- c(
    "genesymbol",
    "entity_type",
    "consensus_score",
    "transmitter",
    "receiver",
    "secreted",
    "plasma_membrane_transmembrane",
    "plasma_membrane_peripheral"
)

if (
    !all(needed_ann %in% names(lig_ann)) ||
    !all(needed_ann %in% names(rec_ann))
) {
    stop(
        "Frozen intercell annotation tables do not contain expected columns."
    )
}

if (
    !all(
        c(
            "gene",
            "logFC",
            "FDR"
        ) %in% names(de)
    )
) {
    stop(
        "DE table does not contain gene/logFC/FDR columns."
    )
}

needed_rob <- c(
    "gene",
    "n_LODO_tested",
    "n_LODO_same_direction",
    "n_LODO_FDR05",
    "all_LODO_direction_concordant",
    "n_donors_same_direction",
    "all_6_donors_same_direction"
)

if (!all(needed_rob %in% names(robustness))) {

    stop(
        "LODO table lacks: ",
        paste(
            setdiff(
                needed_rob,
                names(robustness)
            ),
            collapse = ", "
        )
    )
}

# ============================================================
# Align pseudobulk metadata
# ============================================================

m <- match(
    colnames(pb_counts),
    sample_info$group_key
)

if (anyNA(m)) {
    stop(
        "Could not align pseudobulk counts to sample metadata."
    )
}

sample_info <- sample_info[
    m,
    ,
    drop = FALSE
]

if (!all(
    sample_info$group_key ==
        colnames(pb_counts)
)) {
    stop(
        "Pseudobulk metadata ordering failure."
    )
}

donors <- sort(
    unique(
        sample_info$donor_id
    )
)

if (length(donors) != 6) {
    stop(
        "Expected six donors; found ",
        length(donors),
        "."
    )
}

cat(
    "Frozen OmniPath interaction records:",
    nrow(interactions),
    "\n"
)

cat(
    "Frozen ligand annotation records:",
    nrow(lig_ann),
    "\n"
)

cat(
    "Frozen receptor annotation records:",
    nrow(rec_ann),
    "\n"
)

cat(
    "Successful predefined resources:",
    sum(
        resource_status$status ==
            "SUCCESS"
    ),
    "\n"
)

# ============================================================
# Helpers
# ============================================================

to_bool <- function(x) {

    if (is.logical(x)) {
        return(x)
    }

    z <- tolower(
        trimws(
            as.character(x)
        )
    )

    z %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
    )
}

safe_num <- function(x) {

    suppressWarnings(
        as.numeric(
            as.character(x)
        )
    )
}

collapse_unique <- function(x) {

    x <- trimws(
        as.character(x)
    )

    x <- x[
        !is.na(x) &
        nzchar(x)
    ]

    if (!length(x)) {
        return(NA_character_)
    }

    paste(
        sort(
            unique(x)
        ),
        collapse = " | "
    )
}

safe_max <- function(x) {

    x <- safe_num(x)

    if (
        !length(x) ||
        all(is.na(x))
    ) {
        return(NA_real_)
    }

    max(
        x,
        na.rm = TRUE
    )
}

# ============================================================
# Donor-level expression support
# ============================================================

y <- DGEList(
    counts = pb_counts
)

y <- calcNormFactors(
    y,
    method = "TMM"
)

cpm_mat <- cpm(
    y,
    log = FALSE,
    normalized.lib.sizes = TRUE
)

get_donor_cols <- function(class_name) {

    vapply(
        donors,
        function(d) {

            z <- sample_info$group_key[
                sample_info$donor_id == d &
                sample_info$annotation_class ==
                    class_name
            ]

            if (length(z) != 1) {

                stop(
                    "Expected one pseudobulk for ",
                    d,
                    " / ",
                    class_name
                )
            }

            z
        },
        character(1)
    )
}

nageotte_cols <- get_donor_cols(
    "Nageotte"
)

neuron_cols <- get_donor_cols(
    "Nearby_neuron"
)

expr_summary <- data.frame(

    gene = rownames(cpm_mat),

    mean_CPM_Nageotte =
        rowMeans(
            cpm_mat[
                ,
                nageotte_cols,
                drop = FALSE
            ]
        ),

    mean_CPM_NearbyNeuron =
        rowMeans(
            cpm_mat[
                ,
                neuron_cols,
                drop = FALSE
            ]
        ),

    median_CPM_Nageotte =
        apply(
            cpm_mat[
                ,
                nageotte_cols,
                drop = FALSE
            ],
            1,
            median
        ),

    median_CPM_NearbyNeuron =
        apply(
            cpm_mat[
                ,
                neuron_cols,
                drop = FALSE
            ],
            1,
            median
        ),

    n_donors_CPM1_Nageotte =
        rowSums(
            cpm_mat[
                ,
                nageotte_cols,
                drop = FALSE
            ] >= 1
        ),

    n_donors_CPM1_NearbyNeuron =
        rowSums(
            cpm_mat[
                ,
                neuron_cols,
                drop = FALSE
            ] >= 1
        ),

    stringsAsFactors = FALSE
)

write.csv(
    expr_summary,
    file.path(
        ligand_dir,
        "GSE295206_AllGenes_DonorExpressionSupport.csv"
    ),
    row.names = FALSE
)

# ============================================================
# Prepare OmniPath ligand/receptor annotations
#
# Restrict to simple proteins.
# Complexes are not used here because the DRG transcriptomic
# analysis cannot directly establish complete complex assembly.
# ============================================================

prepare_annotation <- function(
    x,
    role = c(
        "ligand",
        "receptor"
    )
) {

    role <- match.arg(role)

    x$genesymbol <- trimws(
        as.character(
            x$genesymbol
        )
    )

    x$entity_type <- tolower(
        trimws(
            as.character(
                x$entity_type
            )
        )
    )

    x$consensus_score <-
        safe_num(
            x$consensus_score
        )

    x$transmitter_bool <-
        to_bool(
            x$transmitter
        )

    x$receiver_bool <-
        to_bool(
            x$receiver
        )

    x$secreted_bool <-
        to_bool(
            x$secreted
        )

    x$pmtm_bool <-
        to_bool(
            x$plasma_membrane_transmembrane
        )

    x$pmp_bool <-
        to_bool(
            x$plasma_membrane_peripheral
        )

    keep <-
        !is.na(x$genesymbol) &
        nzchar(x$genesymbol) &
        x$entity_type == "protein"

    if (role == "ligand") {

        keep <-
            keep &
            x$transmitter_bool
    }

    if (role == "receptor") {

        keep <-
            keep &
            x$receiver_bool
    }

    x <- x[
        keep,
        ,
        drop = FALSE
    ]

    spl <- split(
        seq_len(
            nrow(x)
        ),
        x$genesymbol
    )

    out <- do.call(
        rbind,
        lapply(
            names(spl),
            function(g) {

                ii <- spl[[g]]

                data.frame(

                    gene = g,

                    consensus_score =
                        safe_max(
                            x$consensus_score[ii]
                        ),

                    secreted =
                        any(
                            x$secreted_bool[ii],
                            na.rm = TRUE
                        ),

                    plasma_membrane_transmembrane =
                        any(
                            x$pmtm_bool[ii],
                            na.rm = TRUE
                        ),

                    plasma_membrane_peripheral =
                        any(
                            x$pmp_bool[ii],
                            na.rm = TRUE
                        ),

                    stringsAsFactors = FALSE
                )
            }
        )
    )

    rownames(out) <- NULL

    out
}

lig_info <- prepare_annotation(
    lig_ann,
    "ligand"
)

rec_info <- prepare_annotation(
    rec_ann,
    "receptor"
)

names(lig_info)[-1] <-
    paste0(
        "ligand_",
        names(lig_info)[-1]
    )

names(rec_info)[-1] <-
    paste0(
        "receptor_",
        names(rec_info)[-1]
    )

cat(
    "Simple-protein ligand genes:",
    nrow(lig_info),
    "\n"
)

cat(
    "Simple-protein receptor genes:",
    nrow(rec_info),
    "\n"
)

# ============================================================
# Construct frozen ligand -> receptor network
#
# IMPORTANT:
# Keep source -> target orientation only.
#
# Source must be independently annotated as a ligand/transmitter.
# Target must be independently annotated as a receptor/receiver.
# ============================================================

interactions$source_genesymbol <-
    trimws(
        as.character(
            interactions$source_genesymbol
        )
    )

interactions$target_genesymbol <-
    trimws(
        as.character(
            interactions$target_genesymbol
        )
    )

lr0 <- interactions[
    interactions$source_genesymbol %in%
        lig_info$gene &
    interactions$target_genesymbol %in%
        rec_info$gene,
    ,
    drop = FALSE
]

if (!nrow(lr0)) {

    stop(
        paste(
            "No ligand-to-receptor records remained",
            "after annotation intersection."
        )
    )
}

lr0$ligand <-
    lr0$source_genesymbol

lr0$receptor <-
    lr0$target_genesymbol

# ============================================================
# Collapse duplicate evidence across the 10 resources
# ============================================================

pair_key <- paste(
    lr0$ligand,
    lr0$receptor,
    sep = "|||PAIR|||"
)

pair_groups <- split(
    seq_len(
        nrow(lr0)
    ),
    pair_key
)

lr_pairs <- do.call(
    rbind,
    lapply(
        pair_groups,
        function(ii) {

            data.frame(

                ligand =
                    lr0$ligand[ii[1]],

                receptor =
                    lr0$receptor[ii[1]],

                n_query_resources =
                    length(
                        unique(
                            lr0$query_resource[ii]
                        )
                    ),

                query_resources =
                    collapse_unique(
                        lr0$query_resource[ii]
                    ),

                max_curation_effort =
                    if (
                        "curation_effort" %in%
                        names(lr0)
                    ) {
                        safe_max(
                            lr0$curation_effort[ii]
                        )
                    } else {
                        NA_real_
                    },

                sources =
                    if (
                        "sources" %in%
                        names(lr0)
                    ) {
                        collapse_unique(
                            lr0$sources[ii]
                        )
                    } else {
                        NA_character_
                    },

                references =
                    if (
                        "references" %in%
                        names(lr0)
                    ) {
                        collapse_unique(
                            lr0$references[ii]
                        )
                    } else {
                        NA_character_
                    },

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(lr_pairs) <- NULL

# ============================================================
# Add ligand/receptor annotation evidence
# ============================================================

lr_pairs <- merge(
    lr_pairs,
    lig_info,
    by.x = "ligand",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
)

lr_pairs <- merge(
    lr_pairs,
    rec_info,
    by.x = "receptor",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
)

lr_pairs$ligand_surface_or_secreted <-
    lr_pairs$ligand_secreted |
    lr_pairs$ligand_plasma_membrane_transmembrane |
    lr_pairs$ligand_plasma_membrane_peripheral

# ============================================================
# Conservative sensitivity network
#
# This is NOT claimed to exactly reproduce
# OmnipathR high_confidence=TRUE.
#
# It applies analogous principles to our frozen snapshots:
# - ligand is secreted or membrane exposed
# - receptor is plasma-membrane transmembrane
# - curation effort >= 2
# - ligand/receptor annotation consensus >= category median
# ============================================================

lig_q50 <- as.numeric(
    quantile(
        lig_info$ligand_consensus_score,
        0.50,
        na.rm = TRUE,
        names = FALSE
    )
)

rec_q50 <- as.numeric(
    quantile(
        rec_info$receptor_consensus_score,
        0.50,
        na.rm = TRUE,
        names = FALSE
    )
)

lr_pairs$high_confidence_sensitivity <-
    lr_pairs$ligand_surface_or_secreted &
    lr_pairs$receptor_plasma_membrane_transmembrane &
    !is.na(
        lr_pairs$max_curation_effort
    ) &
    lr_pairs$max_curation_effort >= 2 &
    !is.na(
        lr_pairs$ligand_consensus_score
    ) &
    lr_pairs$ligand_consensus_score >= lig_q50 &
    !is.na(
        lr_pairs$receptor_consensus_score
    ) &
    lr_pairs$receptor_consensus_score >= rec_q50

lr_pairs <- lr_pairs[
    order(
        lr_pairs$ligand,
        lr_pairs$receptor
    ),
    ,
    drop = FALSE
]

write.csv(
    lr_pairs,
    file.path(
        ligand_dir,
        "OmniPath_FROZEN_LR_network_primary.csv"
    ),
    row.names = FALSE
)

write.csv(
    lr_pairs[
        lr_pairs$high_confidence_sensitivity,
        ,
        drop = FALSE
    ],
    file.path(
        ligand_dir,
        paste0(
            "OmniPath_FROZEN_LR_network_",
            "high_confidence_sensitivity.csv"
        )
    ),
    row.names = FALSE
)

cat(
    "Primary frozen LR pairs:",
    nrow(lr_pairs),
    "\n"
)

cat(
    "Primary frozen LR ligands:",
    length(
        unique(
            lr_pairs$ligand
        )
    ),
    "\n"
)

cat(
    "Primary frozen LR receptors:",
    length(
        unique(
            lr_pairs$receptor
        )
    ),
    "\n"
)

cat(
    "High-confidence sensitivity LR pairs:",
    sum(
        lr_pairs$high_confidence_sensitivity
    ),
    "\n"
)

cat(
    "Ligand annotation q50 consensus score:",
    lig_q50,
    "\n"
)

cat(
    "Receptor annotation q50 consensus score:",
    rec_q50,
    "\n"
)

# ============================================================
# Ligand-level network evidence
# ============================================================

ligand_groups <- split(
    seq_len(
        nrow(lr_pairs)
    ),
    lr_pairs$ligand
)

ligand_evidence <- do.call(
    rbind,
    lapply(
        names(ligand_groups),
        function(g) {

            ii <- ligand_groups[[g]]

            resource_string <- paste(
                lr_pairs$query_resources[ii],
                collapse = " | "
            )

            resource_vector <- unlist(
                strsplit(
                    resource_string,
                    " \\| "
                )
            )

            resource_vector <-
                resource_vector[
                    !is.na(resource_vector) &
                    nzchar(resource_vector)
                ]

            data.frame(

                gene = g,

                n_curated_receptors =
                    length(
                        unique(
                            lr_pairs$receptor[ii]
                        )
                    ),

                n_high_confidence_receptors =
                    length(
                        unique(
                            lr_pairs$receptor[ii][
                                lr_pairs$high_confidence_sensitivity[ii]
                            ]
                        )
                    ),

                n_supporting_query_resources =
                    length(
                        unique(
                            resource_vector
                        )
                    ),

                ligand_consensus_score =
                    safe_max(
                        lr_pairs$ligand_consensus_score[ii]
                    ),

                ligand_secreted =
                    any(
                        lr_pairs$ligand_secreted[ii],
                        na.rm = TRUE
                    ),

                ligand_plasma_membrane_transmembrane =
                    any(
                        lr_pairs$ligand_plasma_membrane_transmembrane[ii],
                        na.rm = TRUE
                    ),

                ligand_plasma_membrane_peripheral =
                    any(
                        lr_pairs$ligand_plasma_membrane_peripheral[ii],
                        na.rm = TRUE
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)

rownames(ligand_evidence) <- NULL

# ============================================================
# Merge DE + LODO + expression + LR evidence
# ============================================================

ligand_de <- de[
    de$gene %in%
        ligand_evidence$gene,
    ,
    drop = FALSE
]

ligand_de <- merge(
    ligand_de,
    robustness[
        ,
        needed_rob,
        drop = FALSE
    ],
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

ligand_de <- merge(
    ligand_de,
    expr_summary,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

ligand_de <- merge(
    ligand_de,
    ligand_evidence,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

# ============================================================
# Spatial interpretation
# ============================================================

ligand_de$enriched_region <-
    ifelse(
        ligand_de$logFC > 0,
        "Nageotte_associated_region",
        "Nearby_neuron_region"
    )

ligand_de$n_donors_CPM1_enriched_region <-
    ifelse(
        ligand_de$logFC > 0,
        ligand_de$n_donors_CPM1_Nageotte,
        ligand_de$n_donors_CPM1_NearbyNeuron
    )

ligand_de$mean_CPM_enriched_region <-
    ifelse(
        ligand_de$logFC > 0,
        ligand_de$mean_CPM_Nageotte,
        ligand_de$mean_CPM_NearbyNeuron
    )

# ============================================================
# Flag signals already highlighted in original GSE295206 work
#
# Keep them as biological positive controls, but do not call
# them novel findings from our reanalysis.
# ============================================================

original_highlighted <- c(
    "SPP1",
    "CLU",
    "CLSTN1",
    "NEGR1",
    "CALCA",
    "GAL",
    "PENK",
    "ADCYAP1"
)

ligand_de$highlighted_original_GSE295206_study <-
    ligand_de$gene %in%
    original_highlighted

# ============================================================
# Candidate tiers
#
# Common:
# FDR < 0.05
# |logFC| >= 0.25
# same direction in all six LODO analyses
# CPM >= 1 in >=4/6 donors in enriched region
#
# Tier A:
# additionally FDR <0.05 in >=4/6 LODO models,
# same direction in all six individual donors,
# >=1 high-confidence sensitivity receptor.
#
# Tier B:
# >=4/6 LODO significance but does not meet every Tier-A
# requirement.
#
# Tier C:
# directionally robust but <4/6 LODO FDR significance.
# ============================================================

base_candidate <-
    ligand_de$FDR < 0.05 &
    abs(
        ligand_de$logFC
    ) >= 0.25 &
    !is.na(
        ligand_de$all_LODO_direction_concordant
    ) &
    ligand_de$all_LODO_direction_concordant &
    ligand_de$n_donors_CPM1_enriched_region >= 4

tier_A <-
    base_candidate &
    ligand_de$n_LODO_FDR05 >= 4 &
    !is.na(
        ligand_de$all_6_donors_same_direction
    ) &
    ligand_de$all_6_donors_same_direction &
    ligand_de$n_high_confidence_receptors >= 1

tier_B <-
    base_candidate &
    ligand_de$n_LODO_FDR05 >= 4 &
    !tier_A

tier_C <-
    base_candidate &
    ligand_de$n_LODO_FDR05 < 4

ligand_de$candidate_tier <-
    "Not_prioritized"

ligand_de$candidate_tier[tier_C] <-
    "Tier_C_directionally_robust"

ligand_de$candidate_tier[tier_B] <-
    "Tier_B_LODO_robust"

ligand_de$candidate_tier[tier_A] <-
    "Tier_A_strict"

ligand_de$novel_cross_tissue_candidate <-
    ligand_de$candidate_tier !=
        "Not_prioritized" &
    !ligand_de$highlighted_original_GSE295206_study

# ============================================================
# Ranking score
#
# Used ONLY for ordering candidates.
# It is not a statistical test.
# ============================================================

ligand_de$ranking_score <-
    pmin(
        -log10(
            pmax(
                ligand_de$FDR,
                1e-300
            )
        ),
        50
    ) +
    abs(
        ligand_de$logFC
    ) +
    ligand_de$n_LODO_FDR05 / 6 +
    ligand_de$n_donors_same_direction / 6 +
    pmin(
        ligand_de$n_high_confidence_receptors,
        5
    ) / 5

ligand_de$ranking_score[
    is.na(
        ligand_de$ranking_score
    )
] <- -Inf

tier_order <- c(
    "Tier_A_strict",
    "Tier_B_LODO_robust",
    "Tier_C_directionally_robust",
    "Not_prioritized"
)

ligand_de$tier_rank <-
    match(
        ligand_de$candidate_tier,
        tier_order
    )

ligand_de <- ligand_de[
    order(
        ligand_de$tier_rank,
        -ligand_de$ranking_score,
        ligand_de$FDR
    ),
    ,
    drop = FALSE
]

ligand_de$tier_rank <- NULL

# ============================================================
# Prioritized candidates
# ============================================================

prioritized <- ligand_de[
    ligand_de$candidate_tier !=
        "Not_prioritized",
    ,
    drop = FALSE
]

novel_prioritized <- prioritized[
    !prioritized$highlighted_original_GSE295206_study,
    ,
    drop = FALSE
]

write.csv(
    ligand_de,
    file.path(
        ligand_dir,
        paste0(
            "Table_P3_All_Curated_Ligands_",
            "with_DE_Robustness.csv"
        )
    ),
    row.names = FALSE
)

write.csv(
    prioritized,
    file.path(
        ligand_dir,
        "Table_P3_Prioritized_Robust_DRG_Ligands.csv"
    ),
    row.names = FALSE
)

write.csv(
    novel_prioritized,
    file.path(
        ligand_dir,
        paste0(
            "Table_P3_Prioritized_",
            "Novel_CrossTissue_Ligands.csv"
        )
    ),
    row.names = FALSE
)

# ============================================================
# LR pairs for prioritized ligands
# ============================================================

priority_pairs <- lr_pairs[
    lr_pairs$ligand %in%
        prioritized$gene,
    ,
    drop = FALSE
]

priority_pairs <- merge(
    priority_pairs,
    prioritized[
        ,
        c(
            "gene",
            "candidate_tier",
            "logFC",
            "FDR",
            "enriched_region",
            "highlighted_original_GSE295206_study"
        ),
        drop = FALSE
    ],
    by.x = "ligand",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
)

priority_pairs <- priority_pairs[
    order(
        match(
            priority_pairs$ligand,
            prioritized$gene
        ),
        priority_pairs$receptor
    ),
    ,
    drop = FALSE
]

write.csv(
    priority_pairs,
    file.path(
        ligand_dir,
        "Table_P3_Prioritized_Ligand_Receptor_Pairs.csv"
    ),
    row.names = FALSE
)

# ============================================================
# Marker sanity check
# ============================================================

sanity_genes <- c(
    "SOX10",
    "FABP7",
    "TYROBP",
    "TREM2",
    "SPP1",
    "APOE",
    "CLU",
    "NEFL",
    "NEFM",
    "NEFH",
    "PRPH",
    "SNAP25",
    "CALCA"
)

sanity <- data.frame(
    gene = sanity_genes,
    stringsAsFactors = FALSE
)

san_de_cols <- intersect(
    c(
        "gene",
        "logFC",
        "FDR",
        "direction"
    ),
    names(de)
)

sanity <- merge(
    sanity,
    de[
        ,
        san_de_cols,
        drop = FALSE
    ],
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

sanity <- merge(
    sanity,
    expr_summary,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
)

sanity <- sanity[
    match(
        sanity_genes,
        sanity$gene
    ),
    ,
    drop = FALSE
]

write.csv(
    sanity,
    file.path(
        ligand_dir,
        "GSE295206_marker_expression_sanity_check.csv"
    ),
    row.names = FALSE
)

# ============================================================
# Frozen network metadata
# ============================================================

network_meta <- data.frame(

    analysis_time =
        as.character(
            Sys.time()
        ),

    frozen_api_retrieval_time =
        if (
            "retrieval_time" %in%
            names(retrieval_meta)
        ) {
            as.character(
                retrieval_meta$retrieval_time[1]
            )
        } else {
            NA_character_
        },

    predefined_resources =
        paste(
            resource_status$resource,
            collapse = ";"
        ),

    n_predefined_resources =
        nrow(
            resource_status
        ),

    raw_combined_interaction_records =
        nrow(
            interactions
        ),

    simple_protein_ligands =
        nrow(
            lig_info
        ),

    simple_protein_receptors =
        nrow(
            rec_info
        ),

    primary_LR_pairs =
        nrow(
            lr_pairs
        ),

    high_confidence_sensitivity_LR_pairs =
        sum(
            lr_pairs$high_confidence_sensitivity
        ),

    ligand_consensus_q50 =
        lig_q50,

    receptor_consensus_q50 =
        rec_q50,

    stringsAsFactors = FALSE
)

write.csv(
    network_meta,
    file.path(
        ligand_dir,
        "OmniPath_FROZEN_LR_network_metadata.csv"
    ),
    row.names = FALSE
)

# ============================================================
# Console summary
# ============================================================

cat("\n============================================\n")
cat("LIGAND DISCOVERY SUMMARY\n")
cat("============================================\n")

cat(
    "DE-tested genes recognized as LR ligands:",
    nrow(ligand_de),
    "\n"
)

cat(
    "Tier A strict:",
    sum(
        ligand_de$candidate_tier ==
            "Tier_A_strict"
    ),
    "\n"
)

cat(
    "Tier B LODO robust:",
    sum(
        ligand_de$candidate_tier ==
            "Tier_B_LODO_robust"
    ),
    "\n"
)

cat(
    "Tier C directionally robust:",
    sum(
        ligand_de$candidate_tier ==
            "Tier_C_directionally_robust"
    ),
    "\n"
)

cat(
    "Prioritized total:",
    nrow(prioritized),
    "\n"
)

cat(
    paste(
        "Prioritized, not highlighted by",
        "original GSE295206 study:"
    ),
    nrow(novel_prioritized),
    "\n"
)

cat("\n============================================\n")
cat("TOP PRIORITIZED LIGANDS\n")
cat("============================================\n")

display_cols <- c(
    "gene",
    "logFC",
    "FDR",
    "enriched_region",
    "candidate_tier",
    "n_LODO_FDR05",
    "n_donors_same_direction",
    "n_donors_CPM1_enriched_region",
    "mean_CPM_enriched_region",
    "ligand_secreted",
    "n_curated_receptors",
    "n_high_confidence_receptors",
    "highlighted_original_GSE295206_study"
)

display_cols <- intersect(
    display_cols,
    names(prioritized)
)

print(
    head(
        prioritized[
            ,
            display_cols,
            drop = FALSE
        ],
        40
    ),
    row.names = FALSE
)

cat("\n============================================\n")
cat("MARKER SANITY CHECK\n")
cat("============================================\n")

print(
    sanity,
    row.names = FALSE
)

# ============================================================
# Reproducibility
# ============================================================

sink(
    file.path(
        repro_dir,
        "03_ligand_discovery_sessionInfo.txt"
    )
)

cat("Analysis time:\n")
print(Sys.time())

cat("\nFrozen OmniPath retrieval metadata:\n")
print(retrieval_meta)

cat("\nFrozen resource status:\n")
print(resource_status)

cat("\nFrozen LR network metadata:\n")
print(network_meta)

cat("\n")
sessionInfo()

sink()

cat("\n============================================\n")
cat("SCRIPT 03 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

