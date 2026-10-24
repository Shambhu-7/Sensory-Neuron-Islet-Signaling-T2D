
# ============================================================
# PAPER 3
# Redox-associated sensory neuron-islet signaling networks
#
# Script 01
# Build validated GSE295206 spatial-transcriptomic objects
#
# Outputs:
#   1. Master spot metadata
#   2. Annotated-spot count matrices by section
#   3. Section/donor annotation summaries
#   4. Barcode matching validation
#   5. QC summaries
#
# Primary biological comparison:
#   Nageotte nodule vs nearby neuron
#
# Overlap spots are retained but excluded from primary comparison.
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

required <- c("Seurat", "Matrix")

missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
    stop(
        "Missing required packages: ",
        paste(missing, collapse = ", "),
        "\nInstall them before rerunning."
    )
}

library(Matrix)

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

raw_dir <- file.path(
    project_dir,
    "data",
    "raw",
    "GSE295206"
)

processed_dir <- file.path(
    project_dir,
    "data",
    "processed"
)

qc_dir <- file.path(
    project_dir,
    "results",
    "qc"
)

repro_dir <- file.path(
    project_dir,
    "reproducibility"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Section / donor manifest
#
# Six subjects and seven DRGs.
# 181030DHB contributes two distinct L4 DRGs.
# ------------------------------------------------------------

manifest <- data.frame(

    GSM = paste0("GSM", 8943583:8943598),

    sample_name = c(
        "sam-UTD-DN0003-DRG-4",
        "sam-UTD-DN0003-DRG-8",
        "sam-UTD-180714DHB-DRGL4-1",
        "sam-UTD-180714DHB-DRGL4-5",
        "sam-UTD-181030DHB-1L4-3",
        "sam-UTD-181030DHB-1L4-7",
        "sam-UTD-181030DHB-2L4-2",
        "sam-UTD-181030DHB-2L4-6",
        "sam-UTD-DN0134-DRGL3-S1-1",
        "sam-UTD-DN0134-DRGL3-S1-2",
        "sam-UTD-DN0134-DRGL3-S1-3",
        "sam-UTD-DN0095-DRGL3-S2-1",
        "sam-UTD-DN0095-DRGL3-S2-2",
        "sam-UTD-DN0095-DRGL3-S2-3",
        "sam-UTD-DN0106-DRGL-S5-1",
        "sam-UTD-DN0106-DRGL-S5-2"
    ),

    donor_id = c(
        "DN0003",
        "DN0003",
        "180714DHB",
        "180714DHB",
        "181030DHB",
        "181030DHB",
        "181030DHB",
        "181030DHB",
        "DN0134",
        "DN0134",
        "DN0134",
        "DN0095",
        "DN0095",
        "DN0095",
        "DN0106",
        "DN0106"
    ),

    drg_id = c(
        "DN0003_DRG",
        "DN0003_DRG",
        "180714DHB_L4",
        "180714DHB_L4",
        "181030DHB_1L4",
        "181030DHB_1L4",
        "181030DHB_2L4",
        "181030DHB_2L4",
        "DN0134_L3",
        "DN0134_L3",
        "DN0134_L3",
        "DN0095_L3",
        "DN0095_L3",
        "DN0095_L3",
        "DN0106_DRG",
        "DN0106_DRG"
    ),

    section_id = sprintf(
        "section_%02d",
        1:16
    ),

    stringsAsFactors = FALSE
)

write.csv(
    manifest,
    file.path(qc_dir, "GSE295206_section_manifest.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------
# Helper: normalize annotation labels
# ------------------------------------------------------------

normalize_annotation <- function(x) {

    # Replace non-breaking spaces
    x <- gsub("\u00A0", " ", x)

    # Collapse repeated whitespace
    x <- gsub("[[:space:]]+", " ", x)

    x <- trimws(x)

    low <- tolower(x)

    result <- rep(NA_character_, length(low))

    # Nageotte nodule
    result[
        low %in% c(
            "nageotte nodule",
            "nageotte nodules"
        )
    ] <- "Nageotte"

    # Explicit overlap
    result[
        low == "nageotte also touching neuron"
    ] <- "Overlap"

    # Nearby neuronal barcodes
    result[
        grepl(
            "neuron.*near nageotte",
            low
        )
    ] <- "Nearby_neuron"

    result
}

# ------------------------------------------------------------
# Helper: read annotation CSV
# ------------------------------------------------------------

read_annotation <- function(path) {

    con <- gzfile(path, open = "rt")

    x <- read.csv(
        con,
        header = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        fill = TRUE
    )

    close(con)

    if (ncol(x) < 2) {
        stop("Annotation file has fewer than two columns: ", path)
    }

    x <- x[, 1:2]

    colnames(x) <- c(
        "barcode",
        "annotation_original"
    )

    x$barcode <- trimws(x$barcode)

    x$annotation_original <- trimws(
        x$annotation_original
    )

    x$annotation_class <- normalize_annotation(
        x$annotation_original
    )

    if (any(is.na(x$annotation_class))) {

        bad <- unique(
            x$annotation_original[
                is.na(x$annotation_class)
            ]
        )

        stop(
            "\nUnrecognized annotation label(s):\n",
            paste(bad, collapse = "\n")
        )
    }

    if (anyDuplicated(x$barcode)) {
        stop(
            "Duplicated annotated barcodes in: ",
            path
        )
    }

    x
}

# ------------------------------------------------------------
# Helper: read Visium positions
# ------------------------------------------------------------

read_positions <- function(path) {

    con <- gzfile(path, open = "rt")

    first_line <- readLines(
        con,
        n = 1
    )

    close(con)

    has_header <- grepl(
        "^barcode,",
        first_line,
        ignore.case = TRUE
    )

    con <- gzfile(path, open = "rt")

    x <- read.csv(
        con,
        header = has_header,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    close(con)

    if (!has_header) {

        if (ncol(x) < 6) {
            stop(
                "Position file has fewer than 6 columns: ",
                path
            )
        }

        colnames(x)[1:6] <- c(
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres"
        )
    }

    # Standardize possible Space Ranger naming variants
    names(x) <- sub(
        "^pxl_row_in_fullres$",
        "pxl_row_in_fullres",
        names(x)
    )

    names(x) <- sub(
        "^pxl_col_in_fullres$",
        "pxl_col_in_fullres",
        names(x)
    )

    if (!"barcode" %in% names(x)) {
        names(x)[1] <- "barcode"
    }

    x$barcode <- trimws(x$barcode)

    x
}

# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

metadata_list <- list()

annotated_counts <- list()

validation_list <- list()

section_summary_list <- list()

# ------------------------------------------------------------
# Process all 16 sections
# ------------------------------------------------------------

for (i in seq_len(nrow(manifest))) {

    gsm <- manifest$GSM[i]

    donor <- manifest$donor_id[i]

    drg <- manifest$drg_id[i]

    section <- manifest$section_id[i]

    sample_name <- manifest$sample_name[i]

    cat("\n")
    cat("====================================================\n")
    cat("PROCESSING ", gsm, "\n", sep = "")
    cat("Donor   : ", donor, "\n", sep = "")
    cat("DRG     : ", drg, "\n", sep = "")
    cat("Section : ", section, "\n", sep = "")
    cat("====================================================\n")

    sample_dir <- file.path(
        raw_dir,
        gsm
    )

    # --------------------------------------------------------
    # Locate files
    # --------------------------------------------------------

    h5_file <- list.files(
        sample_dir,
        pattern = "filtered_feature_bc_matrix\\.h5$",
        full.names = TRUE
    )

    ann_file <- list.files(
        sample_dir,
        pattern = "(NgnNeuSelections\\.csv\\.gz|sample4\\.csv\\.gz|sample8\\.csv\\.gz)$",
        full.names = TRUE
    )

    pos_file <- list.files(
        sample_dir,
        pattern = "tissue_positions.*\\.csv\\.gz$",
        full.names = TRUE
    )

    if (length(h5_file) != 1) {
        stop(
            gsm,
            ": expected exactly one filtered H5 file; found ",
            length(h5_file)
        )
    }

    if (length(ann_file) != 1) {
        stop(
            gsm,
            ": expected exactly one annotation file; found ",
            length(ann_file)
        )
    }

    if (length(pos_file) != 1) {
        stop(
            gsm,
            ": expected exactly one position file; found ",
            length(pos_file)
        )
    }

    # --------------------------------------------------------
    # Read expression matrix
    # --------------------------------------------------------

    cat("Reading count matrix...\n")

    counts <- Seurat::Read10X_h5(
        h5_file,
        use.names = TRUE,
        unique.features = TRUE
    )

    # Some 10X files contain multiple feature types
    if (is.list(counts)) {

        if ("Gene Expression" %in% names(counts)) {

            counts <- counts[["Gene Expression"]]

        } else {

            stop(
                gsm,
                ": H5 contains multiple feature types but no ",
                "'Gene Expression' matrix."
            )
        }
    }

    counts <- as(
        counts,
        "dgCMatrix"
    )

    cat(
        "Genes: ",
        nrow(counts),
        " | filtered spots: ",
        ncol(counts),
        "\n",
        sep = ""
    )

    # --------------------------------------------------------
    # Annotation
    # --------------------------------------------------------

    ann <- read_annotation(
        ann_file
    )

    # --------------------------------------------------------
    # Coordinates
    # --------------------------------------------------------

    pos <- read_positions(
        pos_file
    )

    # --------------------------------------------------------
    # Barcode validation
    # --------------------------------------------------------

    annotated_in_matrix <- ann$barcode %in% colnames(counts)

    annotated_in_positions <- ann$barcode %in% pos$barcode

    matrix_in_positions <- colnames(counts) %in% pos$barcode

    validation <- data.frame(

        GSM = gsm,

        donor_id = donor,

        drg_id = drg,

        section_id = section,

        n_expression_spots = ncol(counts),

        n_annotation_spots = nrow(ann),

        n_annotation_in_matrix =
            sum(annotated_in_matrix),

        n_annotation_missing_matrix =
            sum(!annotated_in_matrix),

        n_annotation_in_positions =
            sum(annotated_in_positions),

        n_annotation_missing_positions =
            sum(!annotated_in_positions),

        n_expression_spots_with_positions =
            sum(matrix_in_positions),

        stringsAsFactors = FALSE
    )

    validation_list[[gsm]] <- validation

    if (any(!annotated_in_matrix)) {

        missing_barcodes <- ann$barcode[
            !annotated_in_matrix
        ]

        stop(
            gsm,
            ": annotation barcodes absent from expression matrix:\n",
            paste(
                missing_barcodes,
                collapse = "\n"
            )
        )
    }

    # --------------------------------------------------------
    # QC metrics for every filtered spot
    # --------------------------------------------------------

    nCount <- Matrix::colSums(
        counts
    )

    nFeature <- Matrix::colSums(
        counts > 0
    )

    mt_idx <- grep(
        "^MT-",
        rownames(counts)
    )

    if (length(mt_idx) > 0) {

        percent_mt <- 100 *
            Matrix::colSums(
                counts[mt_idx, , drop = FALSE]
            ) /
            pmax(nCount, 1)

    } else {

        percent_mt <- rep(
            NA_real_,
            ncol(counts)
        )
    }

    # --------------------------------------------------------
    # Build spot metadata
    # --------------------------------------------------------

    barcodes <- colnames(counts)

    annotation_original <- ann$annotation_original[
        match(
            barcodes,
            ann$barcode
        )
    ]

    annotation_class <- ann$annotation_class[
        match(
            barcodes,
            ann$barcode
        )
    ]

    annotation_class[
        is.na(annotation_class)
    ] <- "Other"

    # Coordinate matching
    pos_index <- match(
        barcodes,
        pos$barcode
    )

    md <- data.frame(

        spot_id = paste(
            gsm,
            barcodes,
            sep = "__"
        ),

        barcode = barcodes,

        GSM = gsm,

        sample_name = sample_name,

        donor_id = donor,

        drg_id = drg,

        section_id = section,

        annotation_original =
            annotation_original,

        annotation_class =
            annotation_class,

        nCount_RNA =
            as.numeric(nCount),

        nFeature_RNA =
            as.numeric(nFeature),

        percent_mt =
            as.numeric(percent_mt),

        stringsAsFactors = FALSE
    )

    # Add spatial coordinates
    coordinate_columns <- setdiff(
        names(pos),
        "barcode"
    )

    for (cc in coordinate_columns) {

        md[[cc]] <- pos[[cc]][
            pos_index
        ]
    }

    metadata_list[[gsm]] <- md

    # --------------------------------------------------------
    # Store annotated expression only
    #
    # Keep:
    # Nageotte
    # Nearby_neuron
    # Overlap
    #
    # "Other" spots remain in metadata but their count matrices
    # are not duplicated into the annotated object.
    # --------------------------------------------------------

    annotated_barcodes <- ann$barcode

    m <- counts[
        ,
        annotated_barcodes,
        drop = FALSE
    ]

    colnames(m) <- paste(
        gsm,
        colnames(m),
        sep = "__"
    )

    annotated_counts[[gsm]] <- m

    # --------------------------------------------------------
    # Section annotation summary
    # --------------------------------------------------------

    tab <- as.data.frame(
        table(
            ann$annotation_class
        ),
        stringsAsFactors = FALSE
    )

    names(tab) <- c(
        "annotation_class",
        "n_spots"
    )

    tab$GSM <- gsm
    tab$donor_id <- donor
    tab$drg_id <- drg
    tab$section_id <- section

    section_summary_list[[gsm]] <- tab

    cat("Annotation counts:\n")
    print(
        table(
            ann$annotation_class
        )
    )

    rm(counts, m)

    invisible(gc())
}

# ------------------------------------------------------------
# Combine metadata
# ------------------------------------------------------------

metadata <- do.call(
    rbind,
    metadata_list
)

rownames(metadata) <- NULL

validation <- do.call(
    rbind,
    validation_list
)

rownames(validation) <- NULL

section_summary <- do.call(
    rbind,
    section_summary_list
)

rownames(section_summary) <- NULL

# ------------------------------------------------------------
# Global annotation validation
# ------------------------------------------------------------

annotated_md <- metadata[
    metadata$annotation_class != "Other",
    ,
    drop = FALSE
]

global_counts <- table(
    annotated_md$annotation_class
)

cat("\n")
cat("====================================================\n")
cat("GLOBAL ANNOTATION COUNTS\n")
cat("====================================================\n")

print(global_counts)

expected <- c(
    Nageotte = 1094L,
    Nearby_neuron = 1087L,
    Overlap = 116L
)

observed <- c(
    Nageotte =
        unname(global_counts["Nageotte"]),

    Nearby_neuron =
        unname(global_counts["Nearby_neuron"]),

    Overlap =
        unname(global_counts["Overlap"])
)

annotation_validation <- data.frame(

    annotation_class =
        names(expected),

    expected =
        as.integer(expected),

    observed =
        as.integer(observed),

    difference =
        as.integer(observed - expected),

    stringsAsFactors = FALSE
)

print(annotation_validation)

if (!all(observed == expected)) {

    stop(
        "\nAnnotation totals do not match expected ",
        "published totals."
    )
}

# ------------------------------------------------------------
# Primary analysis flag
#
# Only Nageotte and Nearby_neuron are used for primary
# comparisons. Overlap is retained for sensitivity/context.
# ------------------------------------------------------------

metadata$primary_analysis <- metadata$annotation_class %in%
    c(
        "Nageotte",
        "Nearby_neuron"
    )

# ------------------------------------------------------------
# Donor-level balance
# ------------------------------------------------------------

primary_md <- metadata[
    metadata$primary_analysis,
    ,
    drop = FALSE
]

donor_balance <- as.data.frame.matrix(
    table(
        primary_md$donor_id,
        primary_md$annotation_class
    )
)

donor_balance$donor_id <- rownames(
    donor_balance
)

rownames(donor_balance) <- NULL

donor_balance <- donor_balance[
    ,
    c(
        "donor_id",
        setdiff(
            names(donor_balance),
            "donor_id"
        )
    ),
    drop = FALSE
]

cat("\n")
cat("====================================================\n")
cat("PRIMARY SPOTS BY DONOR\n")
cat("====================================================\n")

print(donor_balance)

# ------------------------------------------------------------
# DRG-level balance
# ------------------------------------------------------------

drg_balance <- as.data.frame.matrix(
    table(
        primary_md$drg_id,
        primary_md$annotation_class
    )
)

drg_balance$drg_id <- rownames(
    drg_balance
)

rownames(drg_balance) <- NULL

drg_balance <- drg_balance[
    ,
    c(
        "drg_id",
        setdiff(
            names(drg_balance),
            "drg_id"
        )
    ),
    drop = FALSE
]

cat("\n")
cat("====================================================\n")
cat("PRIMARY SPOTS BY DRG\n")
cat("====================================================\n")

print(drg_balance)

# ------------------------------------------------------------
# Section-level QC summary
# ------------------------------------------------------------

qc_summary <- aggregate(

    cbind(
        nCount_RNA,
        nFeature_RNA,
        percent_mt
    ) ~ GSM +
        donor_id +
        drg_id +
        section_id,

    data = metadata,

    FUN = function(x) {
        median(
            x,
            na.rm = TRUE
        )
    }
)

names(qc_summary)[
    names(qc_summary) == "nCount_RNA"
] <- "median_nCount_RNA"

names(qc_summary)[
    names(qc_summary) == "nFeature_RNA"
] <- "median_nFeature_RNA"

names(qc_summary)[
    names(qc_summary) == "percent_mt"
] <- "median_percent_mt"

# ------------------------------------------------------------
# Save processed objects
# ------------------------------------------------------------

saveRDS(
    annotated_counts,
    file.path(
        processed_dir,
        "GSE295206_annotated_counts_by_section.rds"
    )
)

saveRDS(
    metadata,
    file.path(
        processed_dir,
        "GSE295206_master_spot_metadata.rds"
    )
)

write.csv(
    metadata,
    file.path(
        processed_dir,
        "GSE295206_master_spot_metadata.csv"
    ),
    row.names = FALSE
)

write.csv(
    annotated_md,
    file.path(
        processed_dir,
        "GSE295206_annotated_spot_metadata.csv"
    ),
    row.names = FALSE
)

write.csv(
    primary_md,
    file.path(
        processed_dir,
        "GSE295206_primary_Nageotte_neuron_metadata.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Save QC / validation tables
# ------------------------------------------------------------

write.csv(
    validation,
    file.path(
        qc_dir,
        "GSE295206_barcode_validation.csv"
    ),
    row.names = FALSE
)

write.csv(
    section_summary,
    file.path(
        qc_dir,
        "GSE295206_annotation_counts_by_section.csv"
    ),
    row.names = FALSE
)

write.csv(
    annotation_validation,
    file.path(
        qc_dir,
        "GSE295206_annotation_validation.csv"
    ),
    row.names = FALSE
)

write.csv(
    donor_balance,
    file.path(
        qc_dir,
        "GSE295206_primary_spots_by_donor.csv"
    ),
    row.names = FALSE
)

write.csv(
    drg_balance,
    file.path(
        qc_dir,
        "GSE295206_primary_spots_by_DRG.csv"
    ),
    row.names = FALSE
)

write.csv(
    qc_summary,
    file.path(
        qc_dir,
        "GSE295206_section_QC_summary.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Save exact original annotation-label frequencies
# ------------------------------------------------------------

original_labels <- as.data.frame(
    table(
        annotated_md$annotation_original
    ),
    stringsAsFactors = FALSE
)

names(original_labels) <- c(
    "original_annotation",
    "n_spots"
)

original_labels <- original_labels[
    order(
        original_labels$n_spots,
        decreasing = TRUE
    ),
]

write.csv(
    original_labels,
    file.path(
        qc_dir,
        "GSE295206_original_annotation_labels.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

sink(
    file.path(
        repro_dir,
        "01_build_GSE295206_sessionInfo.txt"
    )
)

sessionInfo()

sink()

# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SCRIPT 01 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n\n")

cat(
    "Total filtered spatial spots: ",
    nrow(metadata),
    "\n",
    sep = ""
)

cat(
    "Total annotated spots: ",
    nrow(annotated_md),
    "\n",
    sep = ""
)

cat(
    "Primary-analysis spots: ",
    nrow(primary_md),
    "\n",
    sep = ""
)

cat(
    "Nageotte spots: ",
    sum(
        primary_md$annotation_class ==
            "Nageotte"
    ),
    "\n",
    sep = ""
)

cat(
    "Nearby-neuron spots: ",
    sum(
        primary_md$annotation_class ==
            "Nearby_neuron"
    ),
    "\n",
    sep = ""
)

cat(
    "Overlap spots retained: ",
    sum(
        metadata$annotation_class ==
            "Overlap"
    ),
    "\n",
    sep = ""
)

cat(
    "Donors: ",
    length(
        unique(
            metadata$donor_id
        )
    ),
    "\n",
    sep = ""
)

cat(
    "DRGs: ",
    length(
        unique(
            metadata$drg_id
        )
    ),
    "\n",
    sep = ""
)

cat(
    "Sections: ",
    length(
        unique(
            metadata$GSM
        )
    ),
    "\n\n",
    sep = ""
)

cat(
    "Processed data saved to:\n",
    processed_dir,
    "\n\n",
    sep = ""
)

cat(
    "QC tables saved to:\n",
    qc_dir,
    "\n",
    sep = ""
)

