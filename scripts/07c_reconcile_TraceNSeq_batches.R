
cat("============================================\n")
cat("SCRIPT 07c: RECONCILE TRACE-N-SEQ METADATA\n")
cat("============================================\n\n")

indir <- "data/raw/E-MTAB-12940"

m1 <- readRDS(file.path(indir, "smartseq_metadata_batch1_parsed.rds"))
m2 <- readRDS(file.path(indir, "smartseq_metadata_batch2_parsed.rds"))

cat("Batch 1:", nrow(m1), "cells\n")
cat("Batch 2:", nrow(m2), "cells\n\n")

## ============================================================
## 1. CELL ID OVERLAP
## ============================================================

stopifnot("CellName" %in% names(m1))
stopifnot("CellName" %in% names(m2))

id1 <- as.character(m1$CellName)
id2 <- as.character(m2$CellName)

common <- intersect(id1, id2)

cat("============================================\n")
cat("CELLNAME OVERLAP\n")
cat("============================================\n")

cat("Batch 1 unique CellNames:", length(unique(id1)), "\n")
cat("Batch 2 unique CellNames:", length(unique(id2)), "\n")
cat("Shared CellNames:", length(common), "\n")
cat(
    "Batch 1 fraction shared:",
    round(length(common) / length(unique(id1)), 4),
    "\n"
)
cat(
    "Batch 2 fraction shared:",
    round(length(common) / length(unique(id2)), 4),
    "\n"
)

if (length(common) > 0) {
    cat("\nFirst shared CellNames:\n")
    print(head(common, 20))
}

## ============================================================
## 2. COMPARE SHARED-CELL METADATA
## ============================================================

common_cols <- intersect(names(m1), names(m2))

compare_cols <- intersect(
    c(
        "CellName",
        "CellString",
        "Ganglion",
        "MouseModel",
        "FastBlue",
        "Tumor",
        "CellLine",
        "InnervatingTissue",
        "Replicate",
        "Is_Cell",
        "is_cell",
        "ClusterNew",
        "AnnotatedCellType",
        "AnnotatedSubCellType",
        "AnnotatedCellTypeSharma"
    ),
    common_cols
)

if (length(common) > 0) {

    a <- m1[match(common, m1$CellName), compare_cols, drop = FALSE]
    b <- m2[match(common, m2$CellName), compare_cols, drop = FALSE]

    cat("\n============================================\n")
    cat("SHARED-CELL FIELD AGREEMENT\n")
    cat("============================================\n")

    for (nm in setdiff(compare_cols, "CellName")) {

        x <- as.character(a[[nm]])
        y <- as.character(b[[nm]])

        agree <- (x == y) | (is.na(x) & is.na(y))

        cat(
            nm, ":",
            sum(agree, na.rm = TRUE), "/",
            length(agree),
            "agreement =",
            round(mean(agree, na.rm = TRUE), 4),
            "\n"
        )
    }
}

## ============================================================
## 3. EXACT BATCH-2 CROSS-TABLES
## ============================================================

cat("\n============================================\n")
cat("BATCH 2: GANGLION x TISSUE x QC\n")
cat("============================================\n")

print(
    xtabs(
        ~ Ganglion + InnervatingTissue + Is_Cell,
        data = m2
    )
)

cat("\n============================================\n")
cat("BATCH 2: GANGLION x TUMOR\n")
cat("============================================\n")

print(
    xtabs(
        ~ Ganglion + Tumor,
        data = m2
    )
)

cat("\n============================================\n")
cat("BATCH 2: TISSUE x TUMOR\n")
cat("============================================\n")

print(
    xtabs(
        ~ InnervatingTissue + Tumor,
        data = m2
    )
)

cat("\n============================================\n")
cat("BATCH 2: PANCREAS-LIKE TISSUES x GANGLION x TUMOR\n")
cat("============================================\n")

pan2 <- grepl(
    "pancreas|pankreas",
    m2$InnervatingTissue,
    ignore.case = TRUE
)

print(
    xtabs(
        ~ InnervatingTissue + Ganglion + Tumor,
        data = m2[pan2, , drop = FALSE]
    )
)

## ============================================================
## 4. QC-PASS PANCREAS-LIKE CELLS
## ============================================================

qc1 <- as.logical(m1$Is_Cell)
qc2 <- as.logical(m2$Is_Cell)

pan1 <- grepl(
    "pancreas|pankreas",
    m1$InnervatingTissue,
    ignore.case = TRUE
)

pan2 <- grepl(
    "pancreas|pankreas",
    m2$InnervatingTissue,
    ignore.case = TRUE
)

cat("\n============================================\n")
cat("QC-PASS PANCREAS-LIKE CELLS\n")
cat("============================================\n")

cat("\nBatch 1:\n")
print(
    table(
        Ganglion = m1$Ganglion[pan1 & qc1],
        useNA = "ifany"
    )
)

cat("Total:", sum(pan1 & qc1), "\n")

cat("\nBatch 2:\n")
print(
    table(
        Ganglion = m2$Ganglion[pan2 & qc2],
        useNA = "ifany"
    )
)

cat("Total:", sum(pan2 & qc2), "\n")

## ============================================================
## 5. EXACT NG AUDIT
## ============================================================

ng <- m2$Ganglion == "NG"

cat("\n============================================\n")
cat("BATCH 2 NG AUDIT\n")
cat("============================================\n")

cat("All NG:", sum(ng), "\n")
cat("QC-pass NG:", sum(ng & qc2), "\n")

cat("\nNG tissue distribution:\n")
print(
    table(
        m2$InnervatingTissue[ng],
        useNA = "ifany"
    )
)

cat("\nQC-pass NG tissue distribution:\n")
print(
    table(
        m2$InnervatingTissue[ng & qc2],
        useNA = "ifany"
    )
)

cat("\nNG Tumor field:\n")
print(
    table(
        m2$Tumor[ng],
        useNA = "ifany"
    )
)

cat("\nQC-pass NG Tumor field:\n")
print(
    table(
        m2$Tumor[ng & qc2],
        useNA = "ifany"
    )
)

cat("\nNG tissue x Tumor:\n")
print(
    xtabs(
        ~ InnervatingTissue + Tumor,
        data = m2[ng, , drop = FALSE]
    )
)

## ============================================================
## 6. LOOK FOR PUBLISHED HEALTHY COHORT COUNTS
## ============================================================

cat("\n============================================\n")
cat("TARGET PUBLISHED HEALTHY COHORT\n")
cat("============================================\n")

cat("Published target counts:\n")
cat("DRG = 589\n")
cat("CG  = 689\n")
cat("JNG = 106\n")
cat("Total = 1384\n")

## Candidate descriptors only.
## We are NOT defining the cohort yet.

for (batch_name in c("Batch1", "Batch2")) {

    x <- if (batch_name == "Batch1") m1 else m2
    qc <- as.logical(x$Is_Cell)

    cat("\n", batch_name, "QC-pass by Ganglion:\n", sep = "")

    print(
        table(
            x$Ganglion[qc],
            useNA = "ifany"
        )
    )
}

## ============================================================
## 7. CHECK CELLSTRING PATTERNS
## ============================================================

cat("\n============================================\n")
cat("PANCREAS-RELATED CELLSTRING PATTERNS\n")
cat("============================================\n")

for (batch_name in c("Batch1", "Batch2")) {

    x <- if (batch_name == "Batch1") m1 else m2

    pan <- grepl(
        "pancreas|pankreas",
        paste(
            x$InnervatingTissue,
            x$Tumor,
            x$CellLine,
            x$CellString
        ),
        ignore.case = TRUE
    )

    cat("\n", batch_name, " top CellString patterns:\n", sep = "")

    cs <- as.character(x$CellString[pan])

    ## Remove terminal cell number to expose experiment naming
    cs2 <- gsub("_[0-9]+$", "", cs)

    print(
        head(
            sort(
                table(cs2),
                decreasing = TRUE
            ),
            40
        )
    )
}

## ============================================================
## SAVE AUDIT
## ============================================================

saveRDS(
    list(
        batch1 = m1,
        batch2 = m2,
        shared_CellNames = common
    ),
    file.path(
        indir,
        "E-MTAB-12940_batch_reconciliation.rds"
    )
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07c_EMTAB12940_sessionInfo.txt"
)

cat("\n============================================\n")
cat("SCRIPT 07c COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

