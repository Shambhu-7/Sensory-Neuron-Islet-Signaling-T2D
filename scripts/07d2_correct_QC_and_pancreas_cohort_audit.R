
cat("====================================================\n")
cat("SCRIPT 07d2: CORRECT QC + PANCREAS COHORT AUDIT\n")
cat("====================================================\n\n")

indir <- "data/raw/E-MTAB-12940"

m1 <- readRDS(
    file.path(indir, "smartseq_metadata_batch1_parsed.rds")
)

m2 <- readRDS(
    file.path(indir, "smartseq_metadata_batch2_parsed.rds")
)

m1$Batch <- "Batch1"
m2$Batch <- "Batch2"

## ============================================================
## 1. CORRECT AUTHORS' QC
##
## Metadata mito_percent is FRACTION.
## Seurat percent.mt = 100 * mito_percent.
##
## Therefore:
## percent.mt < 5  <=>  mito_percent < 0.05
## ============================================================

correct_qc <- function(x) {

    ncount <- as.numeric(x$sum)
    nfeat  <- as.numeric(x$detected)
    mito   <- as.numeric(x$mito_percent)

    (
        !is.na(ncount) &
        !is.na(nfeat) &
        !is.na(mito) &
        ncount > 25000 &
        nfeat > 2000 &
        mito < 0.05
    )
}

m1$CorrectAuthorQC <- correct_qc(m1)
m2$CorrectAuthorQC <- correct_qc(m2)

cat("====================================================\n")
cat("CORRECTED AUTHOR QC\n")
cat("====================================================\n")

for (nm in c("Batch1", "Batch2")) {

    x <- if (nm == "Batch1") m1 else m2

    cat("\n", nm, "\n", sep = "")

    cat("Total:", nrow(x), "\n")
    cat("Pass:", sum(x$CorrectAuthorQC), "\n")
    cat("Fail:", sum(!x$CorrectAuthorQC), "\n")

    cat("\nQC-pass by ganglion:\n")

    print(
        table(
            x$Ganglion[x$CorrectAuthorQC],
            useNA = "ifany"
        )
    )

    cat("\nCells failing mitochondrial threshold only:\n")

    mtfail <- (
        as.numeric(x$sum) > 25000 &
        as.numeric(x$detected) > 2000 &
        as.numeric(x$mito_percent) >= 0.05
    )

    print(
        table(
            x$Ganglion[mtfail],
            useNA = "ifany"
        )
    )
}

## ============================================================
## 2. DEPOSITED PERIPHERAL-NEURON ANNOTATIONS
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

for (obj in c("m1", "m2")) {

    x <- get(obj)

    x$PeripheralNeuron <-
        x$AnnotatedCellType %in% peripheral_types

    x$SensoryNeuron <-
        x$AnnotatedCellType %in% sensory_types

    assign(obj, x)
}

## ============================================================
## 3. FUNCTION TO REPORT A CANDIDATE COHORT
## ============================================================

report_cohort <- function(x, hit, label) {

    hit[is.na(hit)] <- FALSE

    cat("\n----------------------------------------------------\n")
    cat(label, "\n")
    cat("----------------------------------------------------\n")

    cat("Raw rows:", sum(hit), "\n")

    cat("\nRaw by ganglion:\n")
    print(
        table(
            x$Ganglion[hit],
            useNA = "ifany"
        )
    )

    q <- hit & x$CorrectAuthorQC

    cat("\nCorrect-QC by ganglion:\n")
    print(
        table(
            x$Ganglion[q],
            useNA = "ifany"
        )
    )

    p <- q & x$PeripheralNeuron

    cat("\nCorrect-QC + peripheral-neuron annotation:\n")
    print(
        table(
            x$Ganglion[p],
            useNA = "ifany"
        )
    )

    s <- q & x$SensoryNeuron

    cat("\nCorrect-QC + sensory-neuron annotation:\n")
    print(
        table(
            x$Ganglion[s],
            useNA = "ifany"
        )
    )

    invisible(NULL)
}

## ============================================================
## 4. BATCH 1 PANCREAS DEFINITIONS
## ============================================================

cat("\n====================================================\n")
cat("BATCH 1 PANCREAS DEFINITIONS\n")
cat("====================================================\n")

report_cohort(
    m1,
    m1$InnervatingTissue == "Pancreas",
    "B1: InnervatingTissue == Pancreas"
)

report_cohort(
    m1,
    m1$Tumor == "Healthy",
    "B1: Tumor == Healthy"
)

report_cohort(
    m1,
    m1$CellLine == "Healthy",
    "B1: CellLine == Healthy"
)

report_cohort(
    m1,
    m1$InnervatingTissue == "Pancreas" &
        m1$Tumor %in% c("Healthy", "None"),
    "B1: Pancreas AND Tumor Healthy/None"
)

## ============================================================
## 5. BATCH 2 PANCREAS DEFINITIONS
## ============================================================

cat("\n====================================================\n")
cat("BATCH 2 PANCREAS DEFINITIONS\n")
cat("====================================================\n")

report_cohort(
    m2,
    m2$InnervatingTissue == "Pancreas",
    "B2: InnervatingTissue == Pancreas"
)

report_cohort(
    m2,
    m2$InnervatingTissue == "Pankreas",
    "B2: InnervatingTissue == Pankreas"
)

report_cohort(
    m2,
    m2$Tumor == "Pancreas",
    "B2: Tumor == Pancreas"
)

report_cohort(
    m2,
    m2$Tumor == "Pankreas",
    "B2: Tumor == Pankreas"
)

## Healthy/control pancreas by explicit CellString
cs2 <- as.character(m2$CellString)

healthy_pancreas_string <-
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
        "pankreatitis|pancreatitis|PACO|KPC|melanom|xenograft|tumor",
        cs2,
        ignore.case = TRUE
    )

report_cohort(
    m2,
    healthy_pancreas_string,
    "B2: explicit healthy/control pancreas CellString"
)

## ============================================================
## 6. PANCREAS × TUMOR × GANGLION
## ============================================================

cat("\n====================================================\n")
cat("BATCH 1: PANCREAS x TUMOR x GANGLION\n")
cat("====================================================\n")

z1 <- m1[
    m1$InnervatingTissue == "Pancreas",
    ,
    drop = FALSE
]

print(
    xtabs(
        ~ Ganglion + Tumor,
        data = z1
    )
)

cat("\n====================================================\n")
cat("BATCH 2: PANCREAS/PANKREAS x TUMOR x GANGLION\n")
cat("====================================================\n")

z2 <- m2[
    m2$InnervatingTissue %in% c(
        "Pancreas",
        "Pankreas"
    ),
    ,
    drop = FALSE
]

print(
    xtabs(
        ~ Ganglion + InnervatingTissue + Tumor,
        data = z2
    )
)

## ============================================================
## 7. HEALTHY-PANCREAS CELLSTRING PATTERNS AFTER CORRECT QC
## ============================================================

cat("\n====================================================\n")
cat("BATCH 2 HEALTHY PANCREAS CELLSTRING PATTERNS\n")
cat("====================================================\n")

idx <- healthy_pancreas_string & m2$CorrectAuthorQC

cs <- as.character(m2$CellString[idx])

## remove terminal cell number only
cs_base <- gsub("_[0-9]+$", "", cs)

print(
    sort(
        table(cs_base),
        decreasing = TRUE
    )
)

## ============================================================
## 8. PUBLISHED TARGET — REFERENCE ONLY
## ============================================================

cat("\n====================================================\n")
cat("PUBLISHED FINAL HEALTHY PANCREAS COHORT\n")
cat("====================================================\n")

cat("DRG = 589\n")
cat("CG  = 689\n")
cat("JNG = 106\n")
cat("TOTAL = 1384\n")

cat("\nIMPORTANT:\n")
cat(
    "These are final published high-quality neuronal counts.\n",
    "We do NOT force deposited metadata to equal these values.\n",
    "The deposited raw/annotated tables include additional experiments\n",
    "and non-neuronal annotations.\n",
    sep = ""
)

## ============================================================
## 9. SAVE CORRECTED METADATA
## ============================================================

dir.create(
    "data/processed",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    list(
        batch1 = m1,
        batch2 = m2,
        batch2_explicit_healthy_pancreas =
            healthy_pancreas_string
    ),
    "data/processed/E-MTAB-12940_correct_QC_metadata.rds"
)

write.csv(
    m1,
    "results/validation/E-MTAB-12940_Batch1_correct_QC.csv",
    row.names = FALSE
)

write.csv(
    m2,
    "results/validation/E-MTAB-12940_Batch2_correct_QC.csv",
    row.names = FALSE
)

writeLines(
    capture.output(sessionInfo()),
    "reproducibility/07d2_EMTAB12940_sessionInfo.txt"
)

cat("\n====================================================\n")
cat("SCRIPT 07d2 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

