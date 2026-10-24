
cat("====================================================\n")
cat("SCRIPT 07e2: TRACE-N-SEQ COUNT/METADATA MAPPING\n")
cat("====================================================\n\n")

indir <- "data/raw/E-MTAB-12940"

m1 <- readRDS(
    file.path(indir, "smartseq_metadata_batch1_parsed.rds")
)

m2 <- readRDS(
    file.path(indir, "smartseq_metadata_batch2_parsed.rds")
)

files <- c(
    batch1 = file.path(indir, "smartseq_raw_counts.txt"),
    batch2 = file.path(indir, "smartseq_raw_counts2.txt")
)

stopifnot(all(file.exists(files)))

## ============================================================
## READ HEADER WITHOUT LOADING THE MATRICES
## ============================================================

get_header <- function(f) {

    x <- readLines(f, n = 1, warn = FALSE)

    nt <- length(strsplit(x, "\t", fixed = TRUE)[[1]])
    nc <- length(strsplit(x, ",", fixed = TRUE)[[1]])

    if (nt >= nc) {
        sep <- "\t"
        type <- "TAB"
    } else {
        sep <- ","
        type <- "COMMA"
    }

    fields <- strsplit(x, sep, fixed = TRUE)[[1]]

    fields <- gsub('^"|"$', "", fields)

    list(
        separator = type,
        fields = fields
    )
}

## ============================================================
## AUDIT ONE MATRIX
## ============================================================

audit <- function(f, meta, label) {

    cat("\n====================================================\n")
    cat(label, "\n")
    cat("====================================================\n")

    h <- get_header(f)

    cat("File:", basename(f), "\n")
    cat("Detected separator:", h$separator, "\n")
    cat("Header fields:", length(h$fields), "\n")

    cat("\nFirst 10 header fields:\n")
    print(head(h$fields, 10))

    ## Assume first field is gene identifier.
    cells <- h$fields[-1]

    cat("\nPutative expression columns:", length(cells), "\n")

    candidates <- list(
        CellName = as.character(meta$CellName),
        CellString = as.character(meta$CellString),
        MetadataRowname = rownames(meta)
    )

    cat("\n===== CELL-ID OVERLAP =====\n")

    overlap_table <- data.frame(
        MetadataField = character(),
        N_metadata = integer(),
        N_count_columns = integer(),
        N_overlap = integer(),
        Fraction_metadata = numeric(),
        Fraction_counts = numeric()
    )

    for (nm in names(candidates)) {

        ids <- candidates[[nm]]

        ov <- intersect(ids, cells)

        tmp <- data.frame(
            MetadataField = nm,
            N_metadata = length(unique(ids)),
            N_count_columns = length(unique(cells)),
            N_overlap = length(ov),
            Fraction_metadata =
                length(ov) / length(unique(ids)),
            Fraction_counts =
                length(ov) / length(unique(cells))
        )

        overlap_table <- rbind(
            overlap_table,
            tmp
        )
    }

    print(overlap_table)

    best <- overlap_table[
        which.max(overlap_table$N_overlap),
        ,
        drop = FALSE
    ]

    cat("\nBest metadata identifier:\n")
    print(best)

    best_ids <- candidates[[best$MetadataField]]

    missing_meta <- setdiff(best_ids, cells)
    extra_counts <- setdiff(cells, best_ids)

    cat("\nMetadata IDs absent from count matrix:",
        length(missing_meta), "\n")

    if (length(missing_meta) > 0) {
        print(head(missing_meta, 20))
    }

    cat("\nCount columns absent from metadata:",
        length(extra_counts), "\n")

    if (length(extra_counts) > 0) {
        print(head(extra_counts, 20))
    }

    cat("\nDuplicate count column IDs:",
        sum(duplicated(cells)), "\n")

    cat("Duplicate metadata IDs using best field:",
        sum(duplicated(best_ids)), "\n")

    invisible(
        list(
            cells = cells,
            best = best,
            overlap = overlap_table
        )
    )
}

a1 <- audit(
    files["batch1"],
    m1,
    "COUNT MATRIX 1 vs METADATA BATCH 1"
)

a2 <- audit(
    files["batch2"],
    m2,
    "COUNT MATRIX 2 vs METADATA BATCH 2"
)

## ============================================================
## CROSS-MAPPING CONTROL
## ============================================================

cat("\n====================================================\n")
cat("CROSS-MAPPING CONTROL\n")
cat("====================================================\n")

cat(
    "Matrix1 CellName overlap with Batch2:",
    length(intersect(a1$cells, m2$CellName)),
    "\n"
)

cat(
    "Matrix2 CellName overlap with Batch1:",
    length(intersect(a2$cells, m1$CellName)),
    "\n"
)

## ============================================================
## CANDIDATE GENE SEARCH WITHOUT LOADING MATRIX
## ============================================================

mouse_candidates <- c(
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
)

cat("\n====================================================\n")
cat("TARGET MOUSE LIGANDS\n")
cat("====================================================\n")

print(mouse_candidates)

writeLines(
    mouse_candidates,
    "results/validation/Paper3_mouse_candidate_ligands.txt"
)

saveRDS(
    list(
        matrix1 = a1,
        matrix2 = a2,
        candidates = mouse_candidates
    ),
    "data/processed/E-MTAB-12940_count_mapping_audit.rds"
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07e2_EMTAB12940_sessionInfo.txt"
)

cat("\n====================================================\n")
cat("SCRIPT 07e2 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

