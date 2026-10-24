
cat("====================================================\n")
cat("SCRIPT 07d: RECONSTRUCT HEALTHY PANCREAS COHORT\n")
cat("====================================================\n\n")

indir <- "data/raw/E-MTAB-12940"
outdir <- "results/validation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create("reproducibility", showWarnings = FALSE)

m1 <- readRDS(
    file.path(indir, "smartseq_metadata_batch1_parsed.rds")
)

m2 <- readRDS(
    file.path(indir, "smartseq_metadata_batch2_parsed.rds")
)

m1$Batch <- "Batch1"
m2$Batch <- "Batch2"

## ============================================================
## 1. VERIFY QC METADATA
## ============================================================

qc_audit <- function(x, label) {

    cat("\n====================================================\n")
    cat(label, ": QC METADATA\n")
    cat("====================================================\n")

    for (v in c("sum", "detected", "mito_percent")) {

        cat("\n", v, "\n", sep = "")

        z <- suppressWarnings(as.numeric(x[[v]]))

        print(
            quantile(
                z,
                probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1),
                na.rm = TRUE
            )
        )
    }
}

qc_audit(m1, "BATCH 1")
qc_audit(m2, "BATCH 2")

## ============================================================
## 2. REPRODUCE AUTHORS' EXPRESSION-DERIVED QC
##
## Published workflow:
## nCount_RNA   > 25000
## nFeature_RNA > 2000
## percent.mt   < 5
##
## Metadata analogues:
## sum
## detected
## mito_percent
## ============================================================

apply_author_qc <- function(x) {

    ncount <- suppressWarnings(as.numeric(x$sum))
    nfeat  <- suppressWarnings(as.numeric(x$detected))
    mito   <- suppressWarnings(as.numeric(x$mito_percent))

    pass <- (
        !is.na(ncount) &
        !is.na(nfeat) &
        !is.na(mito) &
        ncount > 25000 &
        nfeat > 2000 &
        mito < 5
    )

    pass
}

m1$AuthorQC <- apply_author_qc(m1)
m2$AuthorQC <- apply_author_qc(m2)

cat("\n====================================================\n")
cat("AUTHOR-QC COUNTS\n")
cat("====================================================\n")

cat("\nBatch 1:\n")
print(table(m1$AuthorQC, useNA = "ifany"))

cat("\nBatch 2:\n")
print(table(m2$AuthorQC, useNA = "ifany"))

cat("\nBatch 1 QC-pass by ganglion:\n")
print(table(m1$Ganglion[m1$AuthorQC]))

cat("\nBatch 2 QC-pass by ganglion:\n")
print(table(m2$Ganglion[m2$AuthorQC]))

## ============================================================
## 3. COMPARE AuthorQC TO Is_Cell
## ============================================================

cat("\n====================================================\n")
cat("AUTHOR-QC vs Is_Cell\n")
cat("====================================================\n")

cat("\nBatch 1:\n")
print(
    table(
        AuthorQC = m1$AuthorQC,
        Is_Cell = m1$Is_Cell,
        useNA = "ifany"
    )
)

cat("\nBatch 2:\n")
print(
    table(
        AuthorQC = m2$AuthorQC,
        Is_Cell = m2$Is_Cell,
        useNA = "ifany"
    )
)

## ============================================================
## 4. DEFINE NEURON CLASSES
## ============================================================

peripheral_types <- c(
    "Sympathetic noradrenergic neurons",
    "Peripheral sensory non-peptidergic neurons",
    "Peripheral sensory peptidergic neurons",
    "Peripheral sensory neurofilament neurons"
)

sensory_types <- c(
    "Peripheral sensory non-peptidergic neurons",
    "Peripheral sensory peptidergic neurons",
    "Peripheral sensory neurofilament neurons"
)

m1$PeripheralNeuron <- m1$AnnotatedCellType %in% peripheral_types
m2$PeripheralNeuron <- m2$AnnotatedCellType %in% peripheral_types

m1$SensoryNeuron <- m1$AnnotatedCellType %in% sensory_types
m2$SensoryNeuron <- m2$AnnotatedCellType %in% sensory_types

## ============================================================
## 5. CONDITION AUDIT AFTER AUTHORS' QC
## ============================================================

condition_audit <- function(x, label) {

    q <- x[x$AuthorQC, , drop = FALSE]

    cat("\n====================================================\n")
    cat(label, ": CONDITION AUDIT AFTER AUTHOR QC\n")
    cat("====================================================\n")

    cat("\nGanglion:\n")
    print(table(q$Ganglion, useNA = "ifany"))

    cat("\nTumor field:\n")
    print(sort(table(q$Tumor, useNA = "ifany"), decreasing = TRUE))

    cat("\nCellLine:\n")
    print(sort(table(q$CellLine, useNA = "ifany"), decreasing = TRUE))

    cat("\nInnervatingTissue:\n")
    print(
        sort(
            table(q$InnervatingTissue, useNA = "ifany"),
            decreasing = TRUE
        )
    )

    cat("\nGanglion x Tumor:\n")
    print(
        table(
            Ganglion = q$Ganglion,
            Tumor = q$Tumor,
            useNA = "ifany"
        )
    )
}

condition_audit(m1, "BATCH 1")
condition_audit(m2, "BATCH 2")

## ============================================================
## 6. SPECIFIC BATCH-2 PANCREAS CONTROL AUDIT
##
## Key observation:
## Tumor == "Pancreas" gives 106 NG cells before filtering,
## exactly matching published JNG n = 106.
## ============================================================

cat("\n====================================================\n")
cat("BATCH 2: Tumor == Pancreas\n")
cat("====================================================\n")

b2_pan <- m2$Tumor == "Pancreas"

cat("All cells:", sum(b2_pan, na.rm = TRUE), "\n")

cat("\nBy ganglion:\n")
print(table(m2$Ganglion[b2_pan]))

cat("\nAuthor-QC pass:\n")
print(
    table(
        Ganglion = m2$Ganglion[b2_pan & m2$AuthorQC]
    )
)

cat("\nAuthor-QC + peripheral neuron:\n")
print(
    table(
        Ganglion = m2$Ganglion[
            b2_pan &
            m2$AuthorQC &
            m2$PeripheralNeuron
        ]
    )
)

cat("\nAnnotated cell type:\n")
print(
    sort(
        table(
            m2$AnnotatedCellType[b2_pan],
            useNA = "ifany"
        ),
        decreasing = TRUE
    )
)

## ============================================================
## 7. BATCH-2 PANKREAS / CTRL AUDIT
## ============================================================

cat("\n====================================================\n")
cat("BATCH 2: HEALTHY-LIKE LABEL AUDIT\n")
cat("====================================================\n")

healthy_like_values <- c(
    "Pancreas",
    "Pankreas",
    "Ctrl",
    "Healthy",
    "None"
)

for (v in healthy_like_values) {

    hit <- (
        m2$Tumor == v |
        m2$CellLine == v
    )

    hit[is.na(hit)] <- FALSE

    cat("\n--------------------------------------------\n")
    cat("VALUE:", v, "\n")
    cat("--------------------------------------------\n")

    cat("All:\n")
    print(table(m2$Ganglion[hit]))

    cat("Author QC:\n")
    print(table(m2$Ganglion[hit & m2$AuthorQC]))

    cat("Author QC + peripheral:\n")
    print(
        table(
            m2$Ganglion[
                hit &
                m2$AuthorQC &
                m2$PeripheralNeuron
            ]
        )
    )
}

## ============================================================
## 8. BATCH-1 HEALTHY-LIKE LABEL AUDIT
## ============================================================

cat("\n====================================================\n")
cat("BATCH 1: HEALTHY-LIKE LABEL AUDIT\n")
cat("====================================================\n")

for (v in healthy_like_values) {

    hit <- (
        m1$Tumor == v |
        m1$CellLine == v
    )

    hit[is.na(hit)] <- FALSE

    cat("\n--------------------------------------------\n")
    cat("VALUE:", v, "\n")
    cat("--------------------------------------------\n")

    cat("All:\n")
    print(table(m1$Ganglion[hit]))

    cat("Author QC:\n")
    print(table(m1$Ganglion[hit & m1$AuthorQC]))

    cat("Author QC + peripheral:\n")
    print(
        table(
            m1$Ganglion[
                hit &
                m1$AuthorQC &
                m1$PeripheralNeuron
            ]
        )
    )
}

## ============================================================
## 9. CONTROL/PANCREAS CELLSTRING PATTERNS
## ============================================================

audit_cellstrings <- function(x, label) {

    cat("\n====================================================\n")
    cat(label, ": CONTROL/PANCREAS CELLSTRING PATTERNS\n")
    cat("====================================================\n")

    cs <- as.character(x$CellString)

    interesting <- grepl(
        "ctrl|control|healthy|pancreas|pankreas",
        cs,
        ignore.case = TRUE
    )

    cs2 <- cs[interesting & x$AuthorQC]

    cs2 <- gsub("_[0-9]+$", "", cs2)

    print(
        head(
            sort(
                table(cs2),
                decreasing = TRUE
            ),
            100
        )
    )
}

audit_cellstrings(m1, "BATCH 1")
audit_cellstrings(m2, "BATCH 2")

## ============================================================
## 10. PUBLISHED TARGET COMPARISON
## ============================================================

cat("\n====================================================\n")
cat("PUBLISHED HEALTHY-PANCREAS TARGET\n")
cat("====================================================\n")

target <- c(
    DRG = 589,
    CG = 689,
    JNG = 106
)

print(target)

cat("Total:", sum(target), "\n")

cat("\nImportant deposited label correspondence:\n")
cat("DRG -> DRG\n")
cat("CG  -> CG\n")
cat("JNG -> deposited NG arm requires explicit validation\n")

## ============================================================
## 11. SAVE AUGMENTED METADATA
## ============================================================

write.csv(
    m1,
    file.path(
        outdir,
        "E-MTAB-12940_Batch1_metadata_with_AuthorQC.csv"
    ),
    row.names = FALSE
)

write.csv(
    m2,
    file.path(
        outdir,
        "E-MTAB-12940_Batch2_metadata_with_AuthorQC.csv"
    ),
    row.names = FALSE
)

saveRDS(
    list(
        batch1 = m1,
        batch2 = m2,
        published_target = target
    ),
    file.path(
        "data/processed",
        "E-MTAB-12940_metadata_AuthorQC_audit.rds"
    )
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07d_EMTAB12940_sessionInfo.txt"
)

cat("\n====================================================\n")
cat("SCRIPT 07d COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

