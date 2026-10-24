
cat("====================================================\n")
cat("SCRIPT 07f: TRACE-N-SEQ CANDIDATE EXPRESSION\n")
cat("====================================================\n\n")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table first: install.packages('data.table')")
}
library(data.table)

dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## INPUT
## ------------------------------------------------------------

indir <- "data/raw/E-MTAB-12940"

qcobj <- readRDS(
    "data/processed/E-MTAB-12940_correct_QC_metadata.rds"
)

if (!is.list(qcobj) ||
    !all(c("batch1", "batch2") %in% names(qcobj))) {
    stop("Correct-QC RDS does not contain batch1 and batch2.")
}

m1 <- qcobj$batch1
m2 <- qcobj$batch2

f1 <- file.path(indir, "smartseq_raw_counts.txt")
f2 <- file.path(indir, "smartseq_raw_counts2.txt")

stopifnot(file.exists(f1), file.exists(f2))

## ------------------------------------------------------------
## TARGET GENES
## ------------------------------------------------------------

candidates <- c(
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

controls <- c(
    "Prph",
    "Tubb3",
    "Snap25",
    "Nefh",
    "Nefm",
    "Calca",
    "Calcb",
    "Tac1",
    "Mrgprd",
    "Socs2",
    "Piezo2"
)

target_genes <- unique(c(candidates, controls))

cat("Candidate genes:", length(candidates), "\n")
cat("Control genes:", length(controls), "\n")
cat("Total requested genes:", length(target_genes), "\n")

## ------------------------------------------------------------
## READ + VALIDATE ONE MATRIX
## ------------------------------------------------------------

extract_targets <- function(file, meta, sep, label) {

    cat("\n====================================================\n")
    cat("READING ", label, "\n", sep = "")
    cat("====================================================\n")

    x <- fread(
        file,
        sep = sep,
        header = TRUE,
        check.names = FALSE,
        showProgress = TRUE
    )

    cat("Dimensions:", nrow(x), "genes x",
        ncol(x) - 1, "cells\n")

    genes <- as.character(x[[1]])
    cell_ids <- names(x)[-1]

    cat("Unique genes:", length(unique(genes)), "\n")
    cat("Duplicated gene IDs:", sum(duplicated(genes)), "\n")

    stopifnot(
        length(cell_ids) == nrow(meta),
        setequal(cell_ids, as.character(meta$CellName))
    )

    idx <- match(cell_ids, as.character(meta$CellName))
    stopifnot(!anyNA(idx))

    meta <- meta[idx, , drop = FALSE]

    stopifnot(
        identical(cell_ids, as.character(meta$CellName))
    )

    cat("Cell alignment: PASS\n")

    ## Library-size validation
    cat("Calculating matrix library sizes...\n")

    lib_matrix <- vapply(
        x[, -1, with = FALSE],
        function(z) sum(as.numeric(z), na.rm = TRUE),
        numeric(1)
    )

    lib_meta <- as.numeric(meta$sum)

    cor_lib <- suppressWarnings(
        cor(lib_matrix, lib_meta)
    )

    absdiff <- abs(lib_matrix - lib_meta)
    reldiff <- absdiff /
        pmax(abs(lib_matrix), abs(lib_meta), 1)

    cat("\nLibrary-size validation:\n")
    cat("Correlation:", cor_lib, "\n")
    cat("Median absolute difference:",
        median(absdiff), "\n")
    cat("Maximum absolute difference:",
        max(absdiff), "\n")
    cat("Maximum relative difference:",
        max(reldiff), "\n")

    ## Gene presence
    present <- target_genes[target_genes %in% genes]

    cat("\nTarget genes present:",
        length(present), "/", length(target_genes), "\n")

    cat("Candidate genes present:",
        sum(candidates %in% genes), "/",
        length(candidates), "\n")

    cat("\nPresent candidate genes:\n")
    print(candidates[candidates %in% genes])

    cat("\nAbsent candidate genes:\n")
    print(setdiff(candidates, genes))

    ## Extract target genes only
    small <- x[genes %in% target_genes]

    small_genes <- as.character(small[[1]])

    counts <- as.matrix(
        small[, -1, with = FALSE]
    )

    storage.mode(counts) <- "numeric"
    rownames(counts) <- small_genes

    if (anyDuplicated(rownames(counts))) {
        counts <- rowsum(
            counts,
            group = rownames(counts),
            reorder = FALSE
        )
    }

    missing <- setdiff(
        target_genes,
        rownames(counts)
    )

    if (length(missing) > 0) {

        zero <- matrix(
            0,
            nrow = length(missing),
            ncol = ncol(counts),
            dimnames = list(
                missing,
                colnames(counts)
            )
        )

        counts <- rbind(counts, zero)
    }

    counts <- counts[
        target_genes,
        ,
        drop = FALSE
    ]

    libsize <- as.numeric(meta$sum)
    names(libsize) <- cell_ids

    stopifnot(all(libsize > 0))

    cpm <- sweep(
        counts,
        2,
        libsize,
        "/"
    ) * 1e6

    logcpm <- log1p(cpm)

    presence <- data.frame(
        Gene = target_genes,
        Candidate = target_genes %in% candidates,
        Present = target_genes %in% present,
        stringsAsFactors = FALSE
    )

    rm(x, small)
    gc()

    list(
        counts = counts,
        cpm = cpm,
        logcpm = logcpm,
        meta = meta,
        presence = presence
    )
}

## ------------------------------------------------------------
## PROCESS BATCH 1
## ------------------------------------------------------------

b1 <- extract_targets(
    f1,
    m1,
    "\t",
    "BATCH 1"
)

## ------------------------------------------------------------
## PROCESS BATCH 2
## ------------------------------------------------------------

b2 <- extract_targets(
    f2,
    m2,
    ",",
    "BATCH 2"
)

## ------------------------------------------------------------
## GENE PRESENCE
## ------------------------------------------------------------

presence <- merge(
    b1$presence,
    b2$presence,
    by = c("Gene", "Candidate"),
    suffixes = c("_Batch1", "_Batch2"),
    all = TRUE
)

write.csv(
    presence,
    "results/validation/Table_P3_TraceNSeq_GenePresence.csv",
    row.names = FALSE
)

cat("\n====================================================\n")
cat("CANDIDATE GENE PRESENCE\n")
cat("====================================================\n")

print(
    presence[presence$Candidate, ],
    row.names = FALSE
)

## ------------------------------------------------------------
## SENSORY TYPES
## ------------------------------------------------------------

sensory_types <- c(
    "Peripheral sensory non-peptidergic neurons",
    "Peripheral sensory peptidergic neurons",
    "Peripheral sensory neurofilament neurons"
)

## ------------------------------------------------------------
## STRICT HEALTHY PANCREAS COHORT
## ------------------------------------------------------------

keep1 <- (
    b1$meta$CorrectAuthorQC &
    b1$meta$Ganglion == "DRG" &
    b1$meta$InnervatingTissue == "Pancreas" &
    b1$meta$Tumor == "Healthy" &
    b1$meta$AnnotatedCellType %in% sensory_types
)

cs2 <- as.character(
    b2$meta$CellString
)

healthy_pancreas_string <- (
    grepl(
        "pancreas|pankreas",
        cs2,
        ignore.case = TRUE
    ) &
    grepl(
        "ctrl|healthy",
        cs2,
        ignore.case = TRUE
    ) &
    !grepl(
        paste(
            c(
                "pankreatitis",
                "pancreatitis",
                "PACO",
                "KPC",
                "melanom",
                "xenograft",
                "tumor"
            ),
            collapse = "|"
        ),
        cs2,
        ignore.case = TRUE
    )
)

keep2 <- (
    b2$meta$CorrectAuthorQC &
    healthy_pancreas_string &
    b2$meta$Ganglion %in% c("DRG", "NG") &
    b2$meta$AnnotatedCellType %in% sensory_types
)

cat("\n====================================================\n")
cat("STRICT HEALTHY PANCREAS SENSORY COHORT\n")
cat("====================================================\n")

cat("\nBatch 1:\n")
print(
    table(
        b1$meta$Ganglion[keep1],
        b1$meta$AnnotatedCellType[keep1]
    )
)

cat("\nBatch 2:\n")
print(
    table(
        b2$meta$Ganglion[keep2],
        b2$meta$AnnotatedCellType[keep2]
    )
)

cat("\nBatch 1 DRG:", sum(keep1), "\n")

cat(
    "Batch 2 DRG:",
    sum(
        keep2 &
        b2$meta$Ganglion == "DRG"
    ),
    "\n"
)

cat(
    "Batch 2 NG:",
    sum(
        keep2 &
        b2$meta$Ganglion == "NG"
    ),
    "\n"
)

## ------------------------------------------------------------
## CELL-LEVEL LONG TABLE
## ------------------------------------------------------------

make_long <- function(obj, keep, batch) {

    ids <- colnames(obj$counts)[keep]
    md <- obj$meta[keep, , drop = FALSE]

    out <- vector(
        "list",
        length(target_genes)
    )

    for (i in seq_along(target_genes)) {

        g <- target_genes[i]

        out[[i]] <- data.frame(
            Batch = batch,
            CellName = ids,
            Ganglion =
                as.character(md$Ganglion),
            CellType =
                as.character(md$AnnotatedCellType),
            Gene = g,
            Candidate =
                g %in% candidates,
            Count =
                as.numeric(obj$counts[g, ids]),
            CPM =
                as.numeric(obj$cpm[g, ids]),
            log1p_CPM =
                as.numeric(obj$logcpm[g, ids]),
            stringsAsFactors = FALSE
        )
    }

    do.call(rbind, out)
}

long1 <- make_long(
    b1,
    keep1,
    "Batch1"
)

long2 <- make_long(
    b2,
    keep2,
    "Batch2"
)

long <- rbind(
    long1,
    long2
)

write.csv(
    long,
    "results/validation/Table_P3_TraceNSeq_CellLevel_TargetExpression.csv",
    row.names = FALSE
)

## ------------------------------------------------------------
## SUMMARY FUNCTION
## ------------------------------------------------------------

summarize_expression <- function(dat, groups) {

    dt <- as.data.table(dat)

    bycols <- c(
        groups,
        "Gene",
        "Candidate"
    )

    ans <- dt[
        ,
        .(
            N_cells = .N,
            N_detected = sum(Count > 0),
            Pct_detected =
                100 * mean(Count > 0),
            Mean_count =
                mean(Count),
            Median_count =
                median(Count),
            Mean_CPM =
                mean(CPM),
            Median_CPM =
                median(CPM),
            Mean_log1p_CPM =
                mean(log1p_CPM),
            Median_log1p_CPM =
                median(log1p_CPM)
        ),
        by = bycols
    ]

    as.data.frame(ans)
}

## ------------------------------------------------------------
## BATCH / GANGLION SUMMARY
## ------------------------------------------------------------

summary_bg <- summarize_expression(
    long,
    c("Batch", "Ganglion")
)

write.csv(
    summary_bg,
    "results/validation/Table_P3_TraceNSeq_CandidateExpression_ByBatchGanglion.csv",
    row.names = FALSE
)

## ------------------------------------------------------------
## SUBTYPE SUMMARY
## ------------------------------------------------------------

summary_subtype <- summarize_expression(
    long,
    c(
        "Batch",
        "Ganglion",
        "CellType"
    )
)

write.csv(
    summary_subtype,
    "results/validation/Table_P3_TraceNSeq_CandidateExpression_BySubtype.csv",
    row.names = FALSE
)

## ------------------------------------------------------------
## COMBINED DRG
## ------------------------------------------------------------

drg <- long[
    long$Ganglion == "DRG",
    ,
    drop = FALSE
]

summary_drg <- summarize_expression(
    drg,
    character(0)
)

write.csv(
    summary_drg,
    "results/validation/Table_P3_TraceNSeq_CandidateExpression_DRG_Combined.csv",
    row.names = FALSE
)

## ------------------------------------------------------------
## NG ARM
## ------------------------------------------------------------

ng <- long[
    long$Ganglion == "NG",
    ,
    drop = FALSE
]

if (nrow(ng) > 0) {

    summary_ng <- summarize_expression(
        ng,
        character(0)
    )

    write.csv(
        summary_ng,
        "results/validation/Table_P3_TraceNSeq_CandidateExpression_NG.csv",
        row.names = FALSE
    )

} else {

    summary_ng <- data.frame()
}

## ------------------------------------------------------------
## INTEGRATED 13-LIGAND TABLE
## ------------------------------------------------------------

get_metric <- function(
    tab,
    gene,
    field
) {

    z <- tab[
        tab$Gene == gene,
        field,
        drop = TRUE
    ]

    if (length(z) == 0) {
        return(NA_real_)
    }

    as.numeric(z[1])
}

b1_drg <- summary_bg[
    summary_bg$Batch == "Batch1" &
    summary_bg$Ganglion == "DRG",
    ,
    drop = FALSE
]

b2_drg <- summary_bg[
    summary_bg$Batch == "Batch2" &
    summary_bg$Ganglion == "DRG",
    ,
    drop = FALSE
]

integrated <- data.frame(
    Gene = candidates,
    stringsAsFactors = FALSE
)

integrated$DRG_Batch1_PctDetected <- vapply(
    candidates,
    function(g)
        get_metric(
            b1_drg,
            g,
            "Pct_detected"
        ),
    numeric(1)
)

integrated$DRG_Batch2_PctDetected <- vapply(
    candidates,
    function(g)
        get_metric(
            b2_drg,
            g,
            "Pct_detected"
        ),
    numeric(1)
)

integrated$DRG_Combined_PctDetected <- vapply(
    candidates,
    function(g)
        get_metric(
            summary_drg,
            g,
            "Pct_detected"
        ),
    numeric(1)
)

integrated$NG_PctDetected <- vapply(
    candidates,
    function(g) {

        if (nrow(summary_ng) == 0) {
            return(NA_real_)
        }

        get_metric(
            summary_ng,
            g,
            "Pct_detected"
        )
    },
    numeric(1)
)

integrated$DRG_Batch1_MeanLogCPM <- vapply(
    candidates,
    function(g)
        get_metric(
            b1_drg,
            g,
            "Mean_log1p_CPM"
        ),
    numeric(1)
)

integrated$DRG_Batch2_MeanLogCPM <- vapply(
    candidates,
    function(g)
        get_metric(
            b2_drg,
            g,
            "Mean_log1p_CPM"
        ),
    numeric(1)
)

integrated$NG_MeanLogCPM <- vapply(
    candidates,
    function(g) {

        if (nrow(summary_ng) == 0) {
            return(NA_real_)
        }

        get_metric(
            summary_ng,
            g,
            "Mean_log1p_CPM"
        )
    },
    numeric(1)
)

integrated$Detected_Both_DRG_Batches <-
    integrated$DRG_Batch1_PctDetected > 0 &
    integrated$DRG_Batch2_PctDetected > 0

integrated$Detected_DRG_and_NG <-
    integrated$DRG_Combined_PctDetected > 0 &
    integrated$NG_PctDetected > 0

integrated <- integrated[
    order(
        -integrated$DRG_Combined_PctDetected,
        -integrated$NG_PctDetected
    ),
    ,
    drop = FALSE
]

write.csv(
    integrated,
    "results/validation/Table_P3_TraceNSeq_13Ligand_AnatomicalValidation.csv",
    row.names = FALSE
)

## ------------------------------------------------------------
## PRINT MAIN RESULT
## ------------------------------------------------------------

cat("\n====================================================\n")
cat("13-LIGAND ANATOMICAL VALIDATION\n")
cat("====================================================\n")

print(
    integrated,
    row.names = FALSE,
    digits = 4
)

cat(
    "\nCandidates detected in BOTH DRG batches:",
    sum(
        integrated$Detected_Both_DRG_Batches,
        na.rm = TRUE
    ),
    "/",
    nrow(integrated),
    "\n"
)

cat(
    "Candidates detected in DRG AND NG:",
    sum(
        integrated$Detected_DRG_and_NG,
        na.rm = TRUE
    ),
    "/",
    nrow(integrated),
    "\n"
)

## ------------------------------------------------------------
## SUBTYPE RESULTS
## ------------------------------------------------------------

cand_subtype <- summary_subtype[
    summary_subtype$Candidate,
    ,
    drop = FALSE
]

cand_subtype <- cand_subtype[
    order(
        cand_subtype$Batch,
        cand_subtype$Ganglion,
        cand_subtype$CellType,
        -cand_subtype$Pct_detected
    ),
    ,
    drop = FALSE
]

cat("\n====================================================\n")
cat("CANDIDATE EXPRESSION BY SENSORY SUBTYPE\n")
cat("====================================================\n")

print(
    cand_subtype[
        ,
        c(
            "Batch",
            "Ganglion",
            "CellType",
            "Gene",
            "N_cells",
            "Pct_detected",
            "Mean_log1p_CPM"
        )
    ],
    row.names = FALSE,
    digits = 4
)

## ------------------------------------------------------------
## POSITIVE CONTROLS
## ------------------------------------------------------------

controls_drg <- summary_drg[
    !summary_drg$Candidate,
    c(
        "Gene",
        "N_cells",
        "Pct_detected",
        "Mean_log1p_CPM"
    ),
    drop = FALSE
]

cat("\n====================================================\n")
cat("DRG POSITIVE-CONTROL / MARKER CHECK\n")
cat("====================================================\n")

print(
    controls_drg[
        order(-controls_drg$Pct_detected),
        ,
        drop = FALSE
    ],
    row.names = FALSE,
    digits = 4
)

## ------------------------------------------------------------
## SAVE
## ------------------------------------------------------------

saveRDS(
    list(
        candidate_genes = candidates,
        controls = controls,
        batch1_target_counts = b1$counts,
        batch2_target_counts = b2$counts,
        batch1_metadata = b1$meta,
        batch2_metadata = b2$meta,
        keep_batch1 = keep1,
        keep_batch2 = keep2,
        cell_level_expression = long,
        integrated_validation = integrated,
        subtype_summary = cand_subtype
    ),
    "data/processed/E-MTAB-12940_Paper3_candidate_expression.rds"
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07f_EMTAB12940_sessionInfo.txt"
)

cat("\n====================================================\n")
cat("SCRIPT 07f COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

