
# ============================================================
# PAPER 3
# Redox-associated sensory neuron-islet signaling networks
#
# Script 02
# Donor-aware differential expression and redox pathway analysis
#
# PRIMARY CONTRAST:
#   Nageotte-associated spatial transcriptome
#        versus
#   nearby-neuron spatial transcriptome
#
# Positive logFC / CAMERA "Up":
#   enriched in Nageotte-associated spatial transcriptome
#
# Negative logFC / CAMERA "Down":
#   enriched in nearby-neuron spatial transcriptome
#
# IMPORTANT INTERPRETATION:
# Nageotte nodules are predominantly non-neuronal structures.
# Therefore this analysis identifies degeneration-associated
# spatial programs, NOT neuron-intrinsic redox changes.
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Install required packages if missing
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
    )
}

bioc_pkgs <- c(
    "edgeR",
    "limma"
)

for (p in bioc_pkgs) {

    if (!requireNamespace(p, quietly = TRUE)) {

        BiocManager::install(
            p,
            ask = FALSE,
            update = FALSE
        )
    }
}

cran_pkgs <- c(
    "msigdbr"
)

for (p in cran_pkgs) {

    if (!requireNamespace(p, quietly = TRUE)) {

        install.packages(
            p,
            repos = "https://cloud.r-project.org"
        )
    }
}

suppressPackageStartupMessages({
    library(Matrix)
    library(edgeR)
    library(limma)
    library(msigdbr)
})

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

processed_dir <- file.path(
    project_dir,
    "data",
    "processed"
)

redox_dir <- file.path(
    project_dir,
    "results",
    "redox"
)

figure_dir <- file.path(
    project_dir,
    "figures"
)

repro_dir <- file.path(
    project_dir,
    "reproducibility"
)

dir.create(
    redox_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    repro_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# Load Script 01 objects
# ------------------------------------------------------------

counts_list <- readRDS(
    file.path(
        processed_dir,
        "GSE295206_annotated_counts_by_section.rds"
    )
)

metadata <- readRDS(
    file.path(
        processed_dir,
        "GSE295206_master_spot_metadata.rds"
    )
)

cat("\n============================================\n")
cat("INPUT VALIDATION\n")
cat("============================================\n")

cat(
    "Sections in count object:",
    length(counts_list),
    "\n"
)

cat(
    "Metadata rows:",
    nrow(metadata),
    "\n"
)

# ------------------------------------------------------------
# Restrict to genes represented in ALL 16 sections
#
# This is important because the first two sections were aligned
# against a feature reference containing 58,288 features whereas
# the remaining sections contain 36,601.
#
# We do NOT fill unavailable genes with artificial zeros.
# ------------------------------------------------------------

gene_numbers <- vapply(
    counts_list,
    nrow,
    integer(1)
)

cat("\nGenes per section:\n")
print(gene_numbers)

common_genes <- Reduce(
    intersect,
    lapply(
        counts_list,
        rownames
    )
)

cat(
    "\nGenes shared across all 16 sections:",
    length(common_genes),
    "\n"
)

writeLines(
    common_genes,
    file.path(
        redox_dir,
        "GSE295206_common_genes_all16.txt"
    )
)

# ------------------------------------------------------------
# Align matrices to common gene space
# ------------------------------------------------------------

aligned_counts <- lapply(
    counts_list,
    function(x) {

        x[
            common_genes,
            ,
            drop = FALSE
        ]
    }
)

counts_all <- do.call(
    cbind,
    aligned_counts
)

if (anyDuplicated(colnames(counts_all))) {
    stop("Duplicated spot IDs after matrix combination.")
}

cat(
    "Combined annotated matrix:",
    nrow(counts_all),
    "genes x",
    ncol(counts_all),
    "spots\n"
)

# ------------------------------------------------------------
# Primary Nageotte / nearby-neuron spots
# ------------------------------------------------------------

primary_md <- metadata[
    metadata$annotation_class %in%
        c(
            "Nageotte",
            "Nearby_neuron"
        ),
    ,
    drop = FALSE
]

idx <- match(
    colnames(counts_all),
    primary_md$spot_id
)

keep_columns <- !is.na(idx)

primary_counts <- counts_all[
    ,
    keep_columns,
    drop = FALSE
]

primary_md <- primary_md[
    match(
        colnames(primary_counts),
        primary_md$spot_id
    ),
    ,
    drop = FALSE
]

if (!all(
    colnames(primary_counts) ==
        primary_md$spot_id
)) {
    stop("Expression/metadata ordering failure.")
}

if (ncol(primary_counts) != 2181) {
    stop(
        "Expected 2181 primary spots; observed ",
        ncol(primary_counts)
    )
}

cat(
    "Primary spots:",
    ncol(primary_counts),
    "\n"
)

# ------------------------------------------------------------
# Construct DONOR x CLASS pseudobulks
#
# Biological replication = donor.
#
# This prevents donors with hundreds of annotated spots from
# being treated as hundreds of independent biological samples.
# ------------------------------------------------------------

primary_md$group_key <- paste(
    primary_md$donor_id,
    primary_md$annotation_class,
    sep = "__"
)

sample_info <- unique(
    primary_md[
        ,
        c(
            "group_key",
            "donor_id",
            "annotation_class"
        )
    ]
)

sample_info <- sample_info[
    order(
        sample_info$donor_id,
        sample_info$annotation_class
    ),
    ,
    drop = FALSE
]

rownames(sample_info) <- NULL

# Verify each donor contributes both classes

donor_class_table <- table(
    sample_info$donor_id,
    sample_info$annotation_class
)

cat("\nDonor x class pseudobulk design:\n")
print(donor_class_table)

if (!all(donor_class_table == 1)) {

    stop(
        "Each donor must contribute exactly one pseudobulk ",
        "for each anatomical class."
    )
}

# Aggregate raw counts

pb_list <- lapply(
    sample_info$group_key,
    function(g) {

        cols <- which(
            primary_md$group_key == g
        )

        Matrix::rowSums(
            primary_counts[
                ,
                cols,
                drop = FALSE
            ]
        )
    }
)

pb_counts <- do.call(
    cbind,
    pb_list
)

rownames(pb_counts) <- common_genes

colnames(pb_counts) <- sample_info$group_key

# Number of Visium spots contributing to each pseudobulk

sample_info$n_spots <- vapply(
    sample_info$group_key,
    function(g) {
        sum(primary_md$group_key == g)
    },
    integer(1)
)

sample_info$library_size <- as.numeric(
    colSums(pb_counts)
)

write.csv(
    sample_info,
    file.path(
        redox_dir,
        "GSE295206_donor_pseudobulk_sample_info.csv"
    ),
    row.names = FALSE
)

saveRDS(
    pb_counts,
    file.path(
        processed_dir,
        "GSE295206_donor_class_pseudobulk_counts.rds"
    )
)

cat(
    "\nPseudobulk matrix:",
    nrow(pb_counts),
    "genes x",
    ncol(pb_counts),
    "samples\n"
)

# ------------------------------------------------------------
# edgeR model
#
# Nearby_neuron = reference.
# Therefore positive logFC = Nageotte-associated enrichment.
# ------------------------------------------------------------

sample_info$annotation_class <- factor(
    sample_info$annotation_class,
    levels = c(
        "Nearby_neuron",
        "Nageotte"
    )
)

sample_info$donor_id <- factor(
    sample_info$donor_id
)

design <- model.matrix(
    ~ donor_id + annotation_class,
    data = sample_info
)

rownames(design) <- sample_info$group_key

cat("\nDesign matrix:\n")
print(design)

coef_name <- "annotation_classNageotte"

if (!coef_name %in% colnames(design)) {
    stop(
        "Could not locate coefficient: ",
        coef_name
    )
}

coef_index <- which(
    colnames(design) == coef_name
)

write.csv(
    design,
    file.path(
        redox_dir,
        "GSE295206_edgeR_design_matrix.csv"
    )
)

# ------------------------------------------------------------
# edgeR filtering and normalization
# ------------------------------------------------------------

y <- edgeR::DGEList(
    counts = pb_counts
)

keep_gene <- edgeR::filterByExpr(
    y,
    design = design
)

cat(
    "\nGenes before expression filtering:",
    nrow(y),
    "\n"
)

cat(
    "Genes retained by filterByExpr:",
    sum(keep_gene),
    "\n"
)

y <- y[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
]

y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
)

write.csv(
    data.frame(
        gene = common_genes,
        retained_for_DE = keep_gene
    ),
    file.path(
        redox_dir,
        "GSE295206_gene_filter_status.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Estimate dispersion and fit quasi-likelihood model
# ------------------------------------------------------------

y <- edgeR::estimateDisp(
    y,
    design,
    robust = TRUE
)

fit <- edgeR::glmQLFit(
    y,
    design,
    robust = TRUE
)

qlf <- edgeR::glmQLFTest(
    fit,
    coef = coef_index
)

de <- edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
)$table

de$gene <- rownames(de)

de <- de[
    ,
    c(
        "gene",
        setdiff(
            names(de),
            "gene"
        )
    )
]

de$direction <- ifelse(
    de$logFC > 0,
    "Nageotte_enriched",
    "Nearby_neuron_enriched"
)

de$significant_FDR05 <- de$FDR < 0.05

de$significant_FDR05_abslogFC025 <-
    de$FDR < 0.05 &
    abs(de$logFC) >= 0.25

write.csv(
    de,
    file.path(
        redox_dir,
        "Table_P3_DE_Nageotte_vs_NearbyNeuron_edgeR.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

n_sig <- sum(
    de$FDR < 0.05,
    na.rm = TRUE
)

n_up <- sum(
    de$FDR < 0.05 &
        de$logFC > 0,
    na.rm = TRUE
)

n_down <- sum(
    de$FDR < 0.05 &
        de$logFC < 0,
    na.rm = TRUE
)

n_effect <- sum(
    de$FDR < 0.05 &
        abs(de$logFC) >= 0.25,
    na.rm = TRUE
)

cat("\n============================================\n")
cat("PRIMARY DIFFERENTIAL EXPRESSION\n")
cat("============================================\n")

cat(
    "FDR < 0.05:",
    n_sig,
    "\n"
)

cat(
    "Nageotte enriched:",
    n_up,
    "\n"
)

cat(
    "Nearby-neuron enriched:",
    n_down,
    "\n"
)

cat(
    "FDR < 0.05 & |logFC| >= 0.25:",
    n_effect,
    "\n"
)

# ------------------------------------------------------------
# Redox / stress anchor genes
#
# These are for transparent gene-level inspection only.
# They are NOT used as an ad-hoc statistical pathway test.
# ------------------------------------------------------------

anchor_genes <- c(

    # NRF2 / antioxidant
    "NFE2L2",
    "KEAP1",
    "HMOX1",
    "NQO1",
    "GCLC",
    "GCLM",
    "GSR",
    "GSS",
    "SLC7A11",
    "TXN",
    "TXNRD1",
    "PRDX1",
    "PRDX2",
    "SRXN1",
    "SOD1",
    "SOD2",
    "CAT",

    # Glutathione / ferroptosis
    "GPX4",
    "ACSL4",
    "LPCAT3",
    "AIFM2",
    "FTH1",
    "FTL",
    "TFRC",
    "CHAC1",

    # ER stress / integrated stress response
    "HSPA5",
    "ATF4",
    "ATF6",
    "XBP1",
    "ERN1",
    "EIF2AK3",
    "DDIT3",
    "PPP1R15A",

    # Mitochondrial stress / quality control
    "OMA1",
    "DELE1",
    "EIF2AK1",
    "YME1L1",
    "HSPD1",
    "HSPE1",
    "HSPA9",
    "LONP1",
    "CLPP",
    "ATF5",
    "SIRT3",

    # Neuroinflammatory / degeneration-associated
    "SPP1",
    "APOE",
    "TYROBP",
    "CLU",
    "CCL2",
    "CX3CL1",
    "PTGS2",
    "TREM2"
)

anchor_idx <- match(
    anchor_genes,
    de$gene
)

anchor_present <- !is.na(anchor_idx)

anchor_de <- de[
    anchor_idx[
        anchor_present
    ],
    ,
    drop = FALSE
]

anchor_de$requested_order <- match(
    anchor_de$gene,
    anchor_genes
)

anchor_de <- anchor_de[
    order(
        anchor_de$requested_order
    ),
    ,
    drop = FALSE
]

anchor_de$requested_order <- NULL

write.csv(
    anchor_de,
    file.path(
        redox_dir,
        "Table_P3_Redox_Stress_Anchor_Genes.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# limma-voom object for correlation-aware CAMERA
# ------------------------------------------------------------

v <- limma::voom(
    y,
    design,
    plot = FALSE
)

# ------------------------------------------------------------
# MSigDB
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LOADING MSigDB GENE SETS\n")
cat("============================================\n")

hallmark_db <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
)

c2_db <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "C2"
)

c5_db <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "C5"
)

# ------------------------------------------------------------
# Helper to convert MSigDB table -> named gene-set list
# ------------------------------------------------------------

make_gene_sets <- function(db) {

    split(
        unique(
            db[
                ,
                c(
                    "gs_name",
                    "gene_symbol"
                )
            ]
        )$gene_symbol,
        unique(
            db[
                ,
                c(
                    "gs_name",
                    "gene_symbol"
                )
            ]
        )$gs_name
    )
}

# Safer implementation preserving correct grouping

make_gene_sets <- function(db) {

    tmp <- unique(
        db[
            ,
            c(
                "gs_name",
                "gene_symbol"
            )
        ]
    )

    split(
        tmp$gene_symbol,
        tmp$gs_name
    )
}

hallmark_sets <- make_gene_sets(
    hallmark_db
)

# ------------------------------------------------------------
# ALL Hallmark CAMERA
# ------------------------------------------------------------

hallmark_index <- lapply(
    hallmark_sets,
    function(gs) {
        which(
            rownames(v) %in% gs
        )
    }
)

hallmark_index <- hallmark_index[
    lengths(hallmark_index) >= 5
]

camera_hallmark_fixed <- limma::camera(
    v,
    index = hallmark_index,
    design = design,
    contrast = coef_index,
    inter.gene.cor = 0.01
)

camera_hallmark_fixed$gene_set <-
    rownames(camera_hallmark_fixed)

rownames(camera_hallmark_fixed) <- NULL

camera_hallmark_fixed$interpretation <-
    ifelse(
        camera_hallmark_fixed$Direction == "Up",
        "Nageotte_enriched",
        "Nearby_neuron_enriched"
    )

camera_hallmark_fixed <- camera_hallmark_fixed[
    ,
    c(
        "gene_set",
        setdiff(
            names(camera_hallmark_fixed),
            "gene_set"
        )
    )
]

write.csv(
    camera_hallmark_fixed,
    file.path(
        redox_dir,
        "Table_P3_CAMERA_ALL_Hallmark_fixedCor001.csv"
    ),
    row.names = FALSE
)

# Correlation-estimated sensitivity analysis

camera_hallmark_estimated <- limma::camera(
    v,
    index = hallmark_index,
    design = design,
    contrast = coef_index,
    inter.gene.cor = NA
)

camera_hallmark_estimated$gene_set <-
    rownames(camera_hallmark_estimated)

rownames(camera_hallmark_estimated) <- NULL

camera_hallmark_estimated$interpretation <-
    ifelse(
        camera_hallmark_estimated$Direction == "Up",
        "Nageotte_enriched",
        "Nearby_neuron_enriched"
    )

camera_hallmark_estimated <- camera_hallmark_estimated[
    ,
    c(
        "gene_set",
        setdiff(
            names(camera_hallmark_estimated),
            "gene_set"
        )
    )
]

write.csv(
    camera_hallmark_estimated,
    file.path(
        redox_dir,
        "Table_P3_CAMERA_ALL_Hallmark_estimatedCorrelation.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Predefined focused redox/stress pathways
# ------------------------------------------------------------

focused_names <- c(

    # ROS / oxidative stress
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
    "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS",

    # NRF2
    "SINGH_NFE2L2_TARGETS",
    "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2",
    "REACTOME_NFE2L2_REGULATING_ANTI_OXIDANT_DETOXIFICATION_ENZYMES",

    # Glutathione
    "REACTOME_GLUTATHIONE_SYNTHESIS_AND_RECYCLING",
    "GOBP_GLUTATHIONE_METABOLIC_PROCESS",

    # Ferroptosis
    "GOBP_FERROPTOSIS",
    "GOBP_POSITIVE_REGULATION_OF_FERROPTOSIS",
    "GOBP_NEGATIVE_REGULATION_OF_FERROPTOSIS",

    # Mitochondria
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "REACTOME_CELLULAR_RESPONSE_TO_MITOCHONDRIAL_STRESS",
    "REACTOME_MITOCHONDRIAL_UNFOLDED_PROTEIN_RESPONSE_UPRMT",
    "REACTOME_MITOPHAGY",

    # ER stress
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
    "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",

    # Inflammatory / neuroinflammatory
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "GOBP_NEUROINFLAMMATORY_RESPONSE",

    # Cell injury/death context
    "HALLMARK_APOPTOSIS"
)

all_selected_db <- rbind(
    hallmark_db[
        hallmark_db$gs_name %in% focused_names,
        ,
        drop = FALSE
    ],
    c2_db[
        c2_db$gs_name %in% focused_names,
        ,
        drop = FALSE
    ],
    c5_db[
        c5_db$gs_name %in% focused_names,
        ,
        drop = FALSE
    ]
)

all_selected_db <- unique(
    all_selected_db
)

found_sets <- sort(
    unique(
        all_selected_db$gs_name
    )
)

missing_sets <- setdiff(
    focused_names,
    found_sets
)

cat(
    "\nFocused pathways requested:",
    length(focused_names),
    "\n"
)

cat(
    "Focused pathways found:",
    length(found_sets),
    "\n"
)

if (length(missing_sets) > 0) {

    cat("\nPathway names not found in current MSigDB:\n")

    cat(
        paste(
            missing_sets,
            collapse = "\n"
        ),
        "\n"
    )
}

writeLines(
    found_sets,
    file.path(
        redox_dir,
        "GSE295206_focused_pathways_found.txt"
    )
)

writeLines(
    missing_sets,
    file.path(
        redox_dir,
        "GSE295206_focused_pathways_missing.txt"
    )
)

# Freeze exact MSigDB membership used locally

membership_out <- unique(
    all_selected_db[
        ,
        intersect(
            c(
                "gs_name",
                "gene_symbol",
                "gs_collection",
                "gs_subcollection",
                "gs_description"
            ),
            names(all_selected_db)
        )
    ]
)

write.csv(
    membership_out,
    file.path(
        redox_dir,
        "GSE295206_focused_MSigDB_membership_used.csv"
    ),
    row.names = FALSE
)

focused_sets <- make_gene_sets(
    all_selected_db
)

focused_index <- lapply(
    focused_sets,
    function(gs) {
        which(
            rownames(v) %in% gs
        )
    }
)

# Require at least five represented genes

focused_index <- focused_index[
    lengths(focused_index) >= 5
]

# Save represented gene-set sizes

represented_sizes <- data.frame(
    gene_set = names(focused_index),
    represented_genes = lengths(focused_index),
    stringsAsFactors = FALSE
)

write.csv(
    represented_sizes,
    file.path(
        redox_dir,
        "GSE295206_focused_pathway_gene_counts.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Focused CAMERA: fixed correlation 0.01
# ------------------------------------------------------------

camera_focused_fixed <- limma::camera(
    v,
    index = focused_index,
    design = design,
    contrast = coef_index,
    inter.gene.cor = 0.01
)

camera_focused_fixed$gene_set <-
    rownames(camera_focused_fixed)

rownames(camera_focused_fixed) <- NULL

camera_focused_fixed$interpretation <-
    ifelse(
        camera_focused_fixed$Direction == "Up",
        "Nageotte_enriched",
        "Nearby_neuron_enriched"
    )

camera_focused_fixed <- camera_focused_fixed[
    ,
    c(
        "gene_set",
        setdiff(
            names(camera_focused_fixed),
            "gene_set"
        )
    )
]

write.csv(
    camera_focused_fixed,
    file.path(
        redox_dir,
        "Table_P3_CAMERA_Focused_Redox_fixedCor001.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Focused CAMERA sensitivity:
# estimate inter-gene correlation separately for each pathway
# ------------------------------------------------------------

camera_focused_estimated <- limma::camera(
    v,
    index = focused_index,
    design = design,
    contrast = coef_index,
    inter.gene.cor = NA
)

camera_focused_estimated$gene_set <-
    rownames(camera_focused_estimated)

rownames(camera_focused_estimated) <- NULL

camera_focused_estimated$interpretation <-
    ifelse(
        camera_focused_estimated$Direction == "Up",
        "Nageotte_enriched",
        "Nearby_neuron_enriched"
    )

camera_focused_estimated <- camera_focused_estimated[
    ,
    c(
        "gene_set",
        setdiff(
            names(camera_focused_estimated),
            "gene_set"
        )
    )
]

write.csv(
    camera_focused_estimated,
    file.path(
        redox_dir,
        "Table_P3_CAMERA_Focused_Redox_estimatedCorrelation.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Compare CAMERA robustness
# ------------------------------------------------------------

camera_compare <- merge(
    camera_focused_fixed[
        ,
        c(
            "gene_set",
            "Direction",
            "PValue",
            "FDR"
        )
    ],
    camera_focused_estimated[
        ,
        c(
            "gene_set",
            "Direction",
            "PValue",
            "FDR"
        )
    ],
    by = "gene_set",
    suffixes = c(
        "_fixed001",
        "_estimated"
    ),
    all = TRUE
)

camera_compare$direction_concordant <-
    camera_compare$Direction_fixed001 ==
    camera_compare$Direction_estimated

write.csv(
    camera_compare,
    file.path(
        redox_dir,
        "Table_P3_CAMERA_Focused_Robustness.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# MDS diagnostic plot
# ------------------------------------------------------------

pdf(
    file.path(
        figure_dir,
        "Figure_P3_Pseudobulk_MDS.pdf"
    ),
    width = 8,
    height = 7
)

pch_vec <- ifelse(
    sample_info$annotation_class == "Nageotte",
    16,
    1
)

limma::plotMDS(
    y,
    labels = as.character(
        sample_info$donor_id
    ),
    pch = pch_vec,
    main = "GSE295206 donor-level pseudobulk MDS"
)

legend(
    "topright",
    legend = c(
        "Nageotte-associated",
        "Nearby neuron"
    ),
    pch = c(
        16,
        1
    ),
    bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Compact text report
# ------------------------------------------------------------

summary_file <- file.path(
    redox_dir,
    "GSE295206_primary_redox_summary.txt"
)

sink(summary_file)

cat(
    "GSE295206 PAPER 3 — PRIMARY REDOX ANALYSIS\n"
)

cat(
    "==========================================\n\n"
)

cat(
    "Biological donors: ",
    length(
        unique(
            sample_info$donor_id
        )
    ),
    "\n",
    sep = ""
)

cat(
    "Pseudobulk samples: ",
    nrow(sample_info),
    "\n",
    sep = ""
)

cat(
    "Primary Visium spots: ",
    ncol(primary_counts),
    "\n",
    sep = ""
)

cat(
    "Genes shared across all sections: ",
    length(common_genes),
    "\n",
    sep = ""
)

cat(
    "Genes retained for edgeR: ",
    nrow(y),
    "\n\n",
    sep = ""
)

cat(
    "Differential genes FDR < 0.05: ",
    n_sig,
    "\n",
    sep = ""
)

cat(
    "Nageotte enriched: ",
    n_up,
    "\n",
    sep = ""
)

cat(
    "Nearby-neuron enriched: ",
    n_down,
    "\n",
    sep = ""
)

cat(
    "FDR < 0.05 & |logFC| >= 0.25: ",
    n_effect,
    "\n\n",
    sep = ""
)

cat(
    "TOP 20 DIFFERENTIAL GENES\n"
)

cat(
    "-------------------------\n"
)

print(
    head(
        de[
            ,
            c(
                "gene",
                "logFC",
                "logCPM",
                "F",
                "PValue",
                "FDR",
                "direction"
            )
        ],
        20
    ),
    row.names = FALSE
)

cat(
    "\n\nFOCUSED REDOX/STRESS CAMERA — fixed correlation 0.01\n"
)

cat(
    "--------------------------------------------------\n"
)

print(
    camera_focused_fixed,
    row.names = FALSE
)

cat(
    "\n\nFOCUSED REDOX/STRESS CAMERA — estimated correlation\n"
)

cat(
    "---------------------------------------------------\n"
)

print(
    camera_focused_estimated,
    row.names = FALSE
)

cat(
    "\n\nREDOX/STRESS ANCHOR GENES\n"
)

cat(
    "-------------------------\n"
)

print(
    anchor_de[
        ,
        c(
            "gene",
            "logFC",
            "logCPM",
            "PValue",
            "FDR",
            "direction"
        )
    ],
    row.names = FALSE
)

sink()

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

sink(
    file.path(
        repro_dir,
        "02_redox_analysis_sessionInfo.txt"
    )
)

cat(
    "Analysis date:\n"
)

print(
    Sys.time()
)

cat(
    "\nmsigdbr package version:\n"
)

print(
    packageVersion(
        "msigdbr"
    )
)

if ("db_version" %in% names(hallmark_db)) {

    cat(
        "\nMSigDB database version(s):\n"
    )

    print(
        unique(
            hallmark_db$db_version
        )
    )
}

cat("\n")

sessionInfo()

sink()

# ------------------------------------------------------------
# Final console output
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SCRIPT 02 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
    "Common genes:",
    length(common_genes),
    "\n"
)

cat(
    "edgeR-tested genes:",
    nrow(y),
    "\n"
)

cat(
    "FDR < 0.05 genes:",
    n_sig,
    "\n"
)

cat(
    "Nageotte enriched:",
    n_up,
    "\n"
)

cat(
    "Nearby-neuron enriched:",
    n_down,
    "\n"
)

cat(
    "Focused pathways tested:",
    length(focused_index),
    "\n\n"
)

cat(
    "Primary report:\n",
    summary_file,
    "\n\n"
)

cat(
    "IMPORTANT:\n",
    "Interpret pathway differences as degeneration-associated ",
    "spatial programs, not neuron-intrinsic changes.\n",
    sep = ""
)

